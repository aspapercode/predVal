Delong <- function(data){
  

  data.trt <- data[data$trt == 1,]
  data.plcb <- data[data$trt == 0,]

  # There needs to be at least one observation of both positive and negative response
  # for the ROC method to work
  if (sum(data.plcb$y) %in% c(0, 1, nrow(data.plcb)) | sum(data.trt$y) %in% c(0, 1, nrow(data.trt))){
    
    res <- list("dAUC.DeLong" = NA, 
                "p.value.DeLong" = NA)
    
  } else{

  roc.analysis.trt <- roc(data.trt$y, data.trt$x, quiet = T, levels = c(0, 1), direction = "<")
  roc.analysis.plcb <- roc(data.plcb$y, data.plcb$x, quiet = T, levels = c(0, 1), direction = "<")

  delong.test <- roc.test(roc.analysis.trt, roc.analysis.plcb, method = "delong", alternative = "greater", paired = F)
  
  res <- list("dAUC.DeLong" = delong.test$estimate[1] - delong.test$estimate[2], 
              "p.value.DeLong" = delong.test$p.value)
  
  }
  
  return(res)
                         
  
}
