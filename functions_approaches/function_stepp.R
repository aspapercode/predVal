stepp_fun = function(data, r1 = 5, r2 = 12, seed=-1) {

library(stepp)
  if(seed > 0 )
    set.seed(seed)
  
# STEPP analysis
swin <- new("stwin", type = "sliding", r1 = r1, r2 = r2)
subp <- new("stsubpop")
subp <- generate(subp, win = swin, cov = data$x)


res <- new("steppes")
modelGLM <- new("stmodelGLM", coltrt = data$trt,
                colY = data$y, glm = "binomial",
                trts = sort(unique(data$trt)))

res <- estimate(res, subp, modelGLM)

  nperm <- 1000#2500
  res <- test(res, nperm)
  p_value_stepp <- res@testresults$Res[[1]]$pvalue


return(list("P-value_stepp"=p_value_stepp
            
))



}
