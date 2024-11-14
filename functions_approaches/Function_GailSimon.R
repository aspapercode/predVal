#install.packages("effects")
GailSimon <- function(data){

  
  gailSimon <- function(thetahat, se){
    stopifnot(is.numeric(thetahat), length(thetahat) > 1,
              !any(is.na(thetahat)), 
              is.numeric(se), length(se) > 1,
              !any(is.na(se)),
              length(thetahat) == length(se))
    
    dname <- paste(deparse1(substitute(thetahat)), "and", 
                   deparse1(substitute(se)))
    
    nSubgroups <- length(thetahat)
    myrange <- 1L:(nSubgroups - 1L)
    Qplus <- sum((thetahat >= 0) * (thetahat/se)^2)
    Qminus <- sum((thetahat < 0) * (thetahat/se)^2)
    minQ <- min(Qplus, Qminus)
    pval <- sum(stats::dbinom(x = myrange, size = nSubgroups - 1L, prob = 0.5) * (1 - stats::pchisq(minQ, df = myrange)))
    return(c("p-value" = pval, "statistic" = minQ))
  }
  
  data$trt=as.factor(data$trt)
  
  subgroupsx <- c(median(data$x))#quantile(data$x, probs = 0.75)#c(median(data$x))
  
  #check if there are sufficient data in each trt group
  
  if(sum(table(data[data$x<=subgroupsx[1],]$trt)>0)==2){
  
    model1 = glm(y ~ trt,  data = data[data$x<=subgroupsx[1],], family = binomial)
    summary(model1)
    effect1 <- effect("trt", model1)
    summary1 <- summary(effect1)
    eff1 <- summary(model1)$coefficients["trt1", "Estimate"]
    eff1_se <- summary(model1)$coefficients["trt1", "Std. Error"]
    
  } else{
    eff1 <- NA
    eff1_se <- NA
  }
  
  if(sum(table(data[data$x>subgroupsx[1],]$trt)>0)==2){
    
    model2 = glm(y ~ trt,  data = data[data$x>subgroupsx[1],], 
                 family = binomial)
    summary(model2)
    effect2 <- effect("trt", model2)
    summary2 <- summary(effect2)
    eff2 <- summary(model2)$coefficients["trt1", "Estimate"]
    eff2_se <- summary(model2)$coefficients["trt1", "Std. Error"]
  } else{
      eff2 <- NA
      eff2_se <- NA
    }
  
  if(length(c(eff1,eff2)[!is.na(c(eff1,eff2))])>1){
    gailsim_pvalue = gailSimon(c(eff1,eff2)[!is.na(c(eff1,eff2))], 
              c(eff1_se, eff2_se)[!is.na(c(eff1,eff2))])["p-value"]
    
  } else{
    
    gailsim_pvalue = NA
  }
  
  return(gailsim_pvalue)

}

