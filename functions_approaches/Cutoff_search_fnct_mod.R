############################################################################

# Function returning the step function with the averages on the left and right hand sides of a cutoff

bar_fnct = function(x, y, cp) {
  temp = x
  temp[which(x <= cp) ] = mean(y[which(x <= cp)])
  temp[which(x > cp) ] = mean(y[which(x > cp)])
  temp
}

###################################################
# Function returning the sum of squares
SSE_fnct = function(cp, x1, y1, x2, y2) {
  sum((y1 - bar_fnct(x1, y1, cp))^2+(y2 - bar_fnct(x2, y2, cp))^2) 
}

###################################################
# Function returning the cutoff minimizing the sum of squares
cutoff_minSSE_fnct = function(grid, x1, y1, x2, y2) {
  SSEvals = sapply(grid, SSE_fnct, x1=x1, y1=y1, x2=x2, y2=y2)
  id = which(SSEvals == min(SSEvals))
  optlims = mean(grid[id])
  return(list("optlims"=optlims,"SSEvals"= SSEvals))
}

