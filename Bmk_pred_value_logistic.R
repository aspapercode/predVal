
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

names_scenarios <- c("NPNPT","NPNPNT", 
                     "HPNP1_50","HPNP1_40", "MPNP1_30", "LPNP1_20",
                     "HPNP2_50", "HPNP2_40", "MPNP2_30", "LPNP2_20",
                     "NPPT","NPPNT", 
                     "HPMP1_50","HPMP1_40", "MPMP1_30", "LPMP1_20",
                     "HPMP2_50", "HPMP2_40", "MPMP2_30", "LPMP2_20",
                     "HPHP1_50","HPHP1_40", "MPHP1_30", "LPHP1_20",
                     "HPHP2_50", "HPHP2_40", "MPHP2_30", "LPHP2_20")

#coefficients for logistic response-biomarker relationships for each scenario
coefficients_b0123 <- rbind(c(-0.405, 0.811, 0, 0),
                            c(-0.405, 0, 0, 0),
                            
                            c(-0.405, -1.428, 0.049, 0),
                            c(-0.405, -0.868, 0.035,0),
                            c(-0.405, -0.45, 0.026,0),
                            c(-0.405, -0.02, 0.017,0),
                            c(-0.405, -1.863, 0.045, 0),
                            c(-0.405, -1.36, 0.034,0),
                            c(-0.405, -0.853, 0.024,0),
                            c(-0.405, -0.465, 0.017,0),
                            
                            c(-1.655, 0.8, 0, 0.026),
                            c(-1.655, 0, 0, 0.026))
colnames(coefficients_b0123) <- c("b0", "b1", "b2", "b3")

coefficients_b0123 <- rbind(coefficients_b0123, coefficients_b0123[3:10, ])
coefficients_b0123[13:20, "b3"] = 0.2*coefficients_b0123[3:10, "b2"]

coefficients_b0123 <- rbind(coefficients_b0123, coefficients_b0123[3:10, ])
coefficients_b0123[21:28, "b3"] = 0.5*coefficients_b0123[3:10, "b2"]

distributionv <- as.character(c(rep("gamma", length(names_scenarios)),rep("gamma", length(names_scenarios)),
                                rep("gamma", length(names_scenarios)), rep("uniform", length(names_scenarios))))
param_distributionv <- c(rep(0.0827,length(names_scenarios)),
                         rep(0.0498,length(names_scenarios)),
                         rep(0.0688,length(names_scenarios)), 
                         rep(100, length(names_scenarios))) #rate for gamma distribution and shape =  650.5434 *rate^2; max value range for uniform

#Data frame with all possible scenarios
df <- data.frame(distributionv, param_distributionv,
                 rbind(coefficients_b0123,coefficients_b0123,coefficients_b0123,coefficients_b0123),
                 (rep(names_scenarios,4)))
colnames(df) <- c("distributionbmk", "param_distribution", "b0", "b1", "b2", "b3", "names_scenarios")

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

#Coefficients of logistic function for the given scenario
b0 = df$b0[task_id]
b1 = df$b1[task_id]
b2 = df$b2[task_id]
b3 = df$b3[task_id]
name_scen = df$names_scenarios[task_id]


# Function for generating biomarker data, with a logistic relationship with the response
SimDataLogistic = function(N, b0, b1, b2, b3, trt, distribution){ 
  
  if(distribution=="gamma"){
    x = matrix(rgamma(N, shape=shape, rate=rate), N, 1)
    p = exp(b0 + b1*(trt) + b2*x*(trt) + b3*x) / (1 + exp(b0 + b1*(trt) + b2*x*(trt) + b3*x)) 
    y = rbinom(N, 1, p)
    mydata  = cbind.data.frame(x,y,p)
    names(mydata) = c("x","y", "p")
    return(list("data"=mydata))
  }
  
  if(distribution=="uniform"){
    x = matrix(runif(N, min = min, max = max), N, 1)
    p = exp(b0 + b1*(trt) + b2*x*(trt) + b3*x) / (1 + exp(b0 + b1*(trt) + b2*x*(trt) + b3*x)) 
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

start = Sys.time()

output <- foreach(i_i = nsim1:nsim2, .packages = c("effects", "detectseparation", "pROC", "mgcv", "stepp", "lmtest"), .combine = rbind) %dorng% { #it returns a list
  
  i = i_i - nsim1+1
  
  # Simulate data
  myData.exp = SimDataLogistic(N.exp, b0=b0, b1=b1, b2=b2,b3=b3,trt=1, distribution = distribution)
  myData.pbo = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
  data=cbind(rbind(myData.exp$data, myData.pbo$data), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo$data$x))))
  
  #GaelleMK
  res_GaelleMK_lin <- Linear(data, 0) #remove seed 
  res_GaelleMK_lin_v <- unlist(res_GaelleMK_lin)

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
  
  out <- c(res_GaelleMK_lin_v,
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

write.csv(data_res, paste0("RawData", name_scen, "_param", df$param_distribution[task_id], 
                           "_distribution_",distribution,  
                           "_ntot", Ntot, "_nsim_", nsim2, ".csv"))


table_cutoffraw = output[, c("cutoff_found", "cutoffvalue","Prob_diff_p_temp")]

write.csv(table_cutoffraw, paste0("RawData", name_scen,"_diffthr", diff_thr, 
                                  "p_thr_",p_thr, 
                                  "_param", df$param_distribution[task_id], 
                                  "_distribution_",distribution, 
                                  "_ntot", Ntot, "_nsim_", nsim2, ".csv"))

