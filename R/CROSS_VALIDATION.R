####################
# CROSS VALIDATION #
####################

#### Function of cross-validation to determine gamma vector
cross_validation = function(training_set, type_var, type_outcome, level_control, D=0, alpha, L=3, S, num_crit, M=5, step, 
nb_sub_cross=5, nsim_cv, ord.bin=10,upper_best=TRUE, M_per_covar=FALSE, seed=42){
    # Grid for gamma values
    grid_gamma = vector("list", L)
    grid_gamma[[1]] = seq(step, 1, by=step)
    if(L > 1){
        for(i in 2:L){
            grid_gamma[[i]] = seq(0, 1, by=step)
        }
    }

    # Divide training set into "nb_sub_cross" sets
    nb_pat = nrow(training_set)
    nb_pat_set = round(nb_pat/nb_sub_cross)
    ind_pat_sets = vector("list", nb_sub_cross)
    ind_pat_remain = 1:nb_pat
    for(i in 1:(nb_sub_cross-1)){
        ind_pat_sets[[i]] = sort(sample(ind_pat_remain, nb_pat_set, replace=FALSE))
        ind_pat_remain = setdiff(ind_pat_remain, ind_pat_sets[[i]])
    }
    ind_pat_sets[[nb_sub_cross]] = sort(ind_pat_remain)

    # Apply subgroup identification
    mat_gamma = expand.grid(grid_gamma)
    nb_gamma= nrow(mat_gamma)
    vec_L = rep(L, nb_gamma)
    for(g in 1:nb_gamma){
        ind_zero = which(mat_gamma[g,]==0)
        if(length(ind_zero)>0){
            first_zero = ind_zero[1]
            vec_L[g] = first_zero-1
            if(first_zero < L){
                mat_gamma[g,(first_zero+1):L] = NA
            }
        }
    }
    doublons = duplicated(mat_gamma,MARGIN=1)
    mat_gamma = mat_gamma[!doublons,]
    nb_gamma= nrow(mat_gamma)
    vec_L = vec_L[!doublons]
    signif_level = numeric(nrow(mat_gamma))
    for(g in 1:nb_gamma){
        pval_g = 0
        for(k in 1:nb_sub_cross){
            ident_g = subgroup_identification_best_cv(training_set[setdiff(1:nb_pat,ind_pat_sets[[k]]),], type_var, 
                 type_outcome, level_control, D, vec_L[g], S, num_crit, M, mat_gamma[g,], alpha, nsim_cv, ord.bin, M_per_covar, seed)
            if(length(ident_g[[2]])>0){
                best_group = ident_g[[1]]
                set_best_group_k = sub_sets_parents(training_set[ind_pat_sets[[k]],], best_group)[[1]]
                new_pval = analyse(set_best_group_k, type_outcome, level_control, D, alpha, upper_best)[2]
                pval_g = pval_g + new_pval
            }
            else{
                pval_g = pval_g + 1
            }
        }
        signif_level[g] =  pval_g/nb_sub_cross 
    }

    # Optimal gamma
    gamma_opt = mat_gamma[which(signif_level==min(signif_level)),]

    return(gamma_opt)
}

