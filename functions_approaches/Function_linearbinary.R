LinearBinary = function(data, seed=-1) {
  
  if(seed > 0 )
    set.seed(seed)
  
    ###################################################
    # Fit a linear model
    
  data_cutoff <- seq(0, 100, by = 2)
  P_linearbinary_inter <- NULL
  
  for(i in 1:length(data_cutoff)){
    
    model = glm(y ~ I(x>data_cutoff[i]) + trt +  trt*I(x>data_cutoff[i]), data = data, family = binomial)
    
    if("I(x > data_cutoff[i])TRUE:trt" %in% rownames(summary(model)$coefficients)){
      P_linearbinary_inter[i] = summary(model)$coefficients["I(x > data_cutoff[i])TRUE:trt", "Pr(>|z|)"]
      } else{
        P_linearbinary_inter[i] = NA
      }
    
    }
    
    return(list("P-value_inter_linbin_v" = dplyr::if_else(is.finite(min(P_linearbinary_inter, na.rm = TRUE)),
                                                          min(P_linearbinary_inter, na.rm = TRUE), NA)
                
    ))
  
  
  
}