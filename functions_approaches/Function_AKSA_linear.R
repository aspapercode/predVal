Linear = function(data, markdiff,seed=-1) {
  
  if(seed > 0 )
    set.seed(seed)
  
  #check separation and convergence of models with given dataset
  
  separationglm <- glm(y ~ x + trt +  trt*x, data = data, family = binomial, 
                       method = "detect_separation")$outcome
  
  convergence <- FALSE
  
  if(separationglm==0){ 
    
    convergence <- TRUE
  
    ###################################################
    # Fit a linear model
    
    model = glm(y ~ x +  trt +  trt*x, data = data, family = binomial)
    BIC_lin = BIC(model)
    
    # Prob that slope is negative (warning: parameter is for trt - control)
    #P_linear = pnorm(-markdiff, mean=-summary(model)$coefficients["x:trt1", "Estimate"], sd=summary(model)$coefficients["x:trt1", "Std. Error"])
    P_linear_inter = summary(model)$coefficients["x:trt", "Pr(>|z|)"]
    
    # Obtain predicted values and their confidence interval
    predicted  = predict(model, se.fit = TRUE)
    fit =  predicted$fit
    se = predicted$se.fit
    lwr = fit - qt(1 - 0.025, df = model$df.residual) * se
    upr = fit + qt(1 - 0.025, df = model$df.residual) * se
    data_gam=cbind(data, fit, lwr, upr)
    data_gam = data_gam[order(data_gam$x), ]

    x = seq(min(data$x), max(data$x), length.out = 100)
    pred1 = predict(model, newdata = data.frame(x, trt = 0), se.fit = TRUE)
    pred2 = predict(model, newdata = data.frame(x, trt = 1), se.fit = TRUE)

    # Calculate the contrast and its confidence interval
    contrast = pred1$fit - pred2$fit
    se.constrat = sqrt(pred1$se.fit^2 + pred2$se.fit^2)
    lwr.contrast = contrast - qnorm(1 - 0.025) * se.constrat
    upr.contrast = contrast + qnorm(1 - 0.025) * se.constrat

    max_contrast = sample(1:length(x), 10000, replace=T)
    min_contrast = sample(1:length(x), 10000, replace=T)
    diff = contrast[pmax(min_contrast, max_contrast)] - contrast[pmin(min_contrast, max_contrast)]
    se.diff = sqrt(se.constrat[min_contrast]^2 + se.constrat[max_contrast]^2)
    # Prob that the Diff bw the contrasts is negative
    y = rnorm(10000, mean=diff, sd=se.diff)
    P_linear =  mean(y < -markdiff)
    
  }
  if(convergence){
    return(list("P(delta<-markdiff)_lin"=P_linear,
                "P-value_inter" = P_linear_inter,
                "BIC_lin"= BIC_lin
                
    ))
  }else {
    return(list("P(delta<-markdiff)_lin"=NA,
                "P-value_inter" = NA,
                "BIC_lin"= NA
                
    ))
  }
  
  
}