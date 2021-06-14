########################################
# RESAMPLING METHOD SIGNIFICANCE LEVEL #
########################################


#### Estimate the adjusted significance level for selection criterion
adjusted_pval_level = function(training_set, promising, nsim, type_var, type_outcome, level_control, D, L=3, S, num_crit, 
M=5, gamma, alpha, ord.bin=10, upper_best=TRUE, M_per_covar=FALSE, seed=42){
    set.seed(seed)
    nb_promising = length(promising[[1]])
    pval_prom = promising[[2]]
    adjusted_pval = NA
    if(nsim > 0 && nb_promising>0){
        nb_pat = nrow(training_set)
        nb_col = ncol(training_set)
        X = as.matrix(training_set[,3:nb_col], ncol=nb_col-3+1, nrow=nb_pat)
        prop_adj_f = numeric(nb_promising)
        for(i in 1:nsim){
            permute = sample(1:nb_pat, nb_pat, replace=FALSE)
            X_perm = X[permute,]
            perm_set = cbind(training_set[,1:2], X_perm)
            if(M_per_covar==TRUE){
                res_promising_childs = subgroup_identification_promising(perm_set, type_var, type_outcome, level_control, D, alpha, 
                                                   L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best)
            } 
            else{
                res_promising_childs = subgroup_identification_promising2(perm_set, type_var, type_outcome, level_control, D, alpha, 
                                                    L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best)
            }
            
            nb_new = length(res_promising_childs[[1]])
            if(nb_new > 0){
                p_min = min(res_promising_childs[[2]])
                # Selection criteria
                for(c in 1:nb_promising){
                    if(p_min <= pval_prom[c]){
                        prop_adj_f[c] = prop_adj_f[c]+1
                    }
                }
            }
        }
        adjusted_pval = prop_adj_f/nsim
    }
    else{
        adjusted_pval = promising[[2]]
    }   
    return(adjusted_pval)
}

