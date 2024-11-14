Cutoff_onedata = function(diff_thr,p_thr,pp_temp,seed=-1) {
  
  if(seed>0)
    set.seed(seed)

subgroup_temp <- logical(length=1)
  
X.all <- data$x
Y.all <- data$y

X.exp.all <- data[data$trt==1,]$x
Y.exp.all <- data[data$trt==1,]$y

X.pbo.all <- data[data$trt==0,]$x
Y.pbo.all <- data[data$trt==0,]$y

grid.all = seq(min(X.all),max(X.all),0.05)
cutoff = cutoff_minSSE_fnct(grid.all, x1 = X.exp.all, y1 = Y.exp.all, x2 = X.pbo.all, y2 = Y.pbo.all)
cp_temp.all = cutoff$optlims

p1_temp.exp.all = mean(Y.exp.all[which(X.exp.all <= cp_temp.all)])
p2_temp.exp.all = mean(Y.exp.all[which(X.exp.all > cp_temp.all)])

p1_temp.pbo.all = mean(Y.pbo.all[which(X.pbo.all <= cp_temp.all)])
p2_temp.pbo.all = mean(Y.pbo.all[which(X.pbo.all > cp_temp.all)])

diff_temp = (p2_temp.exp.all-p2_temp.pbo.all)-(p1_temp.exp.all-p1_temp.pbo.all)

# Check if cut-off satisfies the criterion: [P(\hat{p1}>\hat{p0})>diff_thr]> p_thr in at least pp% of the population

n_biompos <- sum(X.all > cp_temp.all)
prop_temp <- n_biompos/length(X.all)
Prob_diff_p_temp <- NA

if(prop_temp>=pp_temp){
  
  
  a1.exp=sum(Y.exp.all[which(X.exp.all <= cp_temp.all)])
  b1.exp=length(Y.exp.all[which(X.exp.all <= cp_temp.all)])-sum(Y.exp.all[which(X.exp.all <= cp_temp.all)])
  a2.exp=sum(Y.exp.all[which(X.exp.all > cp_temp.all)])
  b2.exp=length(Y.exp.all[which(X.exp.all > cp_temp.all)])-sum(Y.exp.all[which(X.exp.all > cp_temp.all)])
  
  a1.pbo=sum(Y.pbo.all[which(X.pbo.all <= cp_temp.all)])
  b1.pbo=length(Y.pbo.all[which(X.pbo.all <= cp_temp.all)])-sum(Y.pbo.all[which(X.pbo.all <= cp_temp.all)])
  a2.pbo=sum(Y.pbo.all[which(X.pbo.all > cp_temp.all)])
  b2.pbo=length(Y.pbo.all[which(X.pbo.all > cp_temp.all)])-sum(Y.pbo.all[which(X.pbo.all > cp_temp.all)])
  
  Prob_diff_p_temp = mean(((rbeta(10000, a2.exp, b2.exp) - rbeta(10000, a2.pbo, b2.pbo))- 
                            (rbeta(10000, a1.exp, b1.exp)-rbeta(10000, a1.pbo, b1.pbo)))> diff_thr)
  
  if(Prob_diff_p_temp> p_thr){
    subgroup_temp <- TRUE
  }
  
}

prob <- rep(NA, length = length(grid.all))

for(i in 1:length(grid.all)){

a1.exp=sum(Y.exp.all[which(X.exp.all <= grid.all[i])])
b1.exp=length(Y.exp.all[which(X.exp.all <= grid.all[i])])-sum(Y.exp.all[which(X.exp.all <= grid.all[i])])
a2.exp=sum(Y.exp.all[which(X.exp.all > grid.all[i])])
b2.exp=length(Y.exp.all[which(X.exp.all > grid.all[i])])-sum(Y.exp.all[which(X.exp.all > grid.all[i])])

a1.pbo=sum(Y.pbo.all[which(X.pbo.all <= grid.all[i])])
b1.pbo=length(Y.pbo.all[which(X.pbo.all <= grid.all[i])])-sum(Y.pbo.all[which(X.pbo.all <= grid.all[i])])
a2.pbo=sum(Y.pbo.all[which(X.pbo.all > grid.all[i])])
b2.pbo=length(Y.pbo.all[which(X.pbo.all > grid.all[i])])-sum(Y.pbo.all[which(X.pbo.all > grid.all[i])])

prob[i] = mean(((rbeta(10000, a2.exp, b2.exp) - rbeta(10000, a2.pbo, b2.pbo))- 
                           (rbeta(10000, a1.exp, b1.exp)-rbeta(10000, a1.pbo, b1.pbo)))> diff_thr)
}

            

return(list("FoundCutoff"=subgroup_temp, 
            "CutoffValue"=cp_temp.all,
            "diff.all" = diff_temp,
            "prop_temp" = prop_temp,
            "Prob_diff_p_temp" = Prob_diff_p_temp,
            "sumsquares" = cutoff$SSEvals,
            "Prob_diff_p_temp_grid" = prob))
}

