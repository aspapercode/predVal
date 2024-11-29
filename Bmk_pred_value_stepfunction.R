
#### CLEAN WORKING ENVIRONMENT #############################################
rm(list=setdiff(ls(), c()))

#### IMPORT LIBRARIES ######################################################

library(mgcv)
library(splines)
library(stepp)
library(pROC)
library(lmtest)
library(effects)
library(detectseparation)   # to detect infinite estimates of the coeff in the glm model 
# (https://cran.r-project.org/web/packages/brglm2/brglm2.pdf)
library(doParallel)
library(foreach)
library(purrr)
library(doRNG) #for fixing seeds

#Load all functions for the different methods in the folder functions_approaches
source("functions_approaches/Function_AKSA_linear.R")
source("functions_approaches/Function_linearbinary.R")
source("functions_approaches/Function_lrt.R")
source("functions_approaches/Cutoff_search_fnct_mod.R")
source('functions_approaches/FindcutoffOptimizedbinom.R')
source("functions_approaches/Function_GailSimon.R")
source("functions_approaches/function_delong.R")
source("functions_approaches/function_stepp.R")

###################################################

###Array job to run simulations for each scenario

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

set.seed(4219+task_id)

# Build table of scenarios

#Read table with the response rates above and below cutoff for each scenario
rate_unif <- read.csv("rates_unif_step.csv")

names_scenarios <- rate_unif$name_scena

probstep <- rate_unif[, c("rate_below_ctl","rate_above_ctl",
                          "rate_below_exp","rate_above_exp")]

probstep$rate_above_ctl <- ifelse(is.na(probstep$rate_above_ctl), rate_unif$overall_ctl,probstep$rate_above_ctl)
probstep$rate_below_ctl <- ifelse(is.na(probstep$rate_below_ctl), rate_unif$overall_ctl,probstep$rate_below_ctl) 
probstep$rate_above_exp <- ifelse(is.na(probstep$rate_above_exp), rate_unif$overall_exp,probstep$rate_above_exp) 
probstep$rate_below_exp <- ifelse(is.na(probstep$rate_below_exp), rate_unif$overall_exp,probstep$rate_below_exp) 

colnames(probstep) <- c("p1_plb", "p2_plb", 
                        "p1_exp", "p2_exp")

cutoff <- c(1,17, 
            1,17,
            rep(30,3),
            rep(26,3),
            rep( 17,3),
            rep( 1,3),
            rep( 42,3),
            rep( 40,3),
            rep( 35,3),
            rep( 28,3)
)

distributionv <- as.character(c(rep("gamma", length(names_scenarios)),rep("gamma", length(names_scenarios)),
                                rep("gamma", length(names_scenarios)), rep("uniform", length(names_scenarios))))
param_distributionv <- c(rep(0.0827,length(names_scenarios)),
                         rep(0.0498,length(names_scenarios)),
                         rep(0.0688,length(names_scenarios)), 
                         rep(100, length(names_scenarios))) 

#Data frame with all possible scenarios
df <- data.frame(distributionv, param_distributionv,
                 probstep,
                 cutoff,
                 (rep(names_scenarios,4)))
colnames(df) <- c("distributionbmk", "param_distribution", "p1_plb", "p2_plb", 
                  "p1_exp", "p2_exp", "cutoff", "names_scenarios")

#Number of patients on experimental arm 
N.exp = 40
#Number of patients on placebo arm 
N.pbo = 20
#Total number of patients
Ntot = N.exp+N.pbo

#Marked difference
markdiff = 0

#Distribution of the biomarker
distribution = df$distributionbmk[task_id]
#Set parameters of the BMK-distribution
if(distribution == "gamma"){
  rate = df$param_distribution[task_id]
  shape = 650.5434*rate^2
}

if(distribution == "uniform"){
  max = df$param_distribution[task_id]
  min = 0
}

#Take values for specif scenario
cutpoint = df$cutoff[task_id]
p1_exp = df$p1_exp[task_id]
p2_exp = df$p2_exp[task_id]
p1_plb = df$p1_plb[task_id]
p2_plb = df$p2_plb[task_id]
name_scen = df$names_scenarios[task_id]

###################################################
# Function for generating data from Biomarker distr with step function, 1 cutoff
SimDataStep1 = function(N, cutpoint, p2, p1, distribution){
  if(distribution=="gamma"){
    x = matrix(rgamma(N, shape=shape, rate=rate), N, 1)
    p = unlist(lapply(1:length(x), function(j){
      if (x[j]<= cutpoint)  p1
      else  p2
    }) )
    
    y = rbinom(N, 1, p)
    mydata  = cbind.data.frame(x,y,p)
    names(mydata) = c("x","y", "p")
    return(list("data"=mydata))
  }
  
  if(distribution=="uniform"){
    x = matrix(runif(N, min = min, max = max), N, 1)
    p = unlist(lapply(1:length(x), function(j){
      if (x[j]<= cutpoint)  p1
      else  p2
    }) )
    
    y = rbinom(N, 1, p)
    mydata  = cbind.data.frame(x,y,p)
    names(mydata) = c("x","y", "p")
    return(list("data"=mydata))
  }
}


nsim1 = 1
nsim2 = 5000

#Parameters for the different approaches

## for gam model
k=3

#for cutoff identification
diff_thr = 0
p_thr <- 0.8
pp_temp <- 0.2

#for stepp method
r1 = 15 #minimum number of patients overlapping in each subgroup
r2 = 40 #maximum number of patients in each subgroup

#Set the number of cores for the simulations
ncpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
numCores <- makeCluster(ncpus)

registerDoParallel(numCores) #set number of cores

data_saved <- list()

start = Sys.time()

output <- foreach(i_i = nsim1:nsim2, .packages = c("effects", "detectseparation",  "pROC","mgcv", "stepp", "lmtest"), .combine = rbind) %dorng% { #it returns a list
  
  i = i_i - nsim1+1
  
  # Simulate data
  myData.exp = SimDataStep1(N.exp,cutpoint=cutpoint, p2=p2_exp, p1=p1_exp, distribution = distribution)
  myData.pbo = SimDataStep1(N.pbo, cutpoint=cutpoint, p2=p2_plb, p1=p1_plb,distribution = distribution)
  data=cbind(rbind(myData.exp$data, myData.pbo$data), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo$data$x))))
  
  #AKSA
  AKSA_lin <- Linear(data, 0) #remove seed 
  res_AKSA_lin_v <- unlist(AKSA_lin)
  
  #cutoff
  res_onedata <- Cutoff_onedata(diff_thr,p_thr,pp_temp)
  
  #Dichotomous
  res_linbin <- LinearBinary(data) 
  res_linbin_v <- unlist(res_linbin)
  
  #GailSimon
  res_gailsimon <- GailSimon(data) 
  res_gailsimon_v <- unlist(res_gailsimon)
  
  #DeLong
  res_roctest <- Delong(data)
  res_roctest_v <- unlist(res_roctest)
  
  #STEPP
  res_stepp <- stepp_fun(data, r1, r2)
  res_stepp_v <- unlist(res_stepp)
  
  #Likelihood 
  res_lrt <- lrttest(data)
  res_lrt_v <- unlist(res_lrt)
  
  out <- c(res_AKSA_lin_v,
           res_linbin_v, 
           res_gailsimon,
           res_roctest_v,
           res_stepp_v,
           res_lrt_v,
           res_onedata$FoundCutoff,res_onedata$CutoffValue,res_onedata$Prob_diff_p_temp)
  
  out 
  
  
}
end = Sys.time()
end-start 

#give names to output
colnames(output)[which(colnames(output)=="")] <- c("cutoff_found", "cutoffvalue","Prob_diff_p_temp")

data_res <- output

write.csv(data_res, paste0("RawData_StepFun", name_scen, "_param", df$param_distribution[task_id], 
                           "_distribution_",distribution,  
                           "_ntot", Ntot, "_nsim_", nsim2, ".csv"))


table_cutoffraw = output[, c("cutoff_found", "cutoffvalue","Prob_diff_p_temp")]

write.csv(table_cutoffraw, paste0("RawData_StepFun", name_scen,"_diffthr", diff_thr, 
                                  "p_thr_",p_thr, 
                                  "_param", df$param_distribution[task_id], 
                                  "_distribution_",distribution, 
                                  "_ntot", Ntot, "_nsim_", nsim2, ".csv"))

