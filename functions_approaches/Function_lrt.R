lrttest = function(data) {
  
  
  library(lmtest)
  
  model1 = glm(y ~ x + trt , data = data, family = binomial)
  model2 = glm(y ~ x + trt + trt*x, data = data, family = binomial)
  
  p_value_lrtest <- lrtest(model1, model2)$`Pr(>Chisq)`[2]
  
  
  return(list("P-value_lrtest"=p_value_lrtest
                
    ))
  
  
  
}