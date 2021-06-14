###########################
# SUBGROUP IDENTIFICATION #
###########################

#### Function that identify promisings subgroups (select the best split per covariate with a maximum of M splits)
subgroup_identification_promising = function(training_set, type_var, type_outcome, level_control, D=0, alpha, 
L=3, S, num_crit, M=5, gamma, bool_best_pval=FALSE, ord.bin, upper_best=TRUE, modified=TRUE){    
    promisings = list()
    pvalues_promisings = c()
    best_pval = +Inf

    # Recursive function on parent
    rec_search = function(subgroup_parent){
        depth_tree = length(subgroup_parent[[1]])
        # stopping rule
        if(depth_tree >= L) {
            return(0)
        }
        # Remaining covariates and subset from parents
        covariates_not_used = setdiff(3:ncol(training_set), subgroup_parent[[1]]) 
        set_subgroup_parent = sub_sets_parents(training_set, subgroup_parent)[[1]]

        # Update all possible leaves of tree from parents and remaining covariates
        subgroups_child = list()
        pvalues_childs = c()
        splitting_crit_childs = c()
        for(covar in covariates_not_used) {
            vec_z1 = c()
            vec_z2 = c() 
            vec_type = c()
            vec_levels = c()

            # Calculate all combinations of two childs from covariate with pvalues and sample sizes of the corresponding sets
            comb_child = comb_child(set_subgroup_parent[,covar], type_var[covar-2], ord.bin)
            if(anyNA(comb_child)){
                break
            }
            set_all_childs = sub_sets_all_childs(training_set, subgroup_parent, covar, comb_child, type_var[covar-2])
      	    nb_splits = length(comb_child)
            z_stats_eff = vector("list", nb_splits)
      	    pvalues = vector("list", nb_splits)
            sizes = vector("list", nb_splits)
      	    for(split in 1:nb_splits){ 
                analyse_temp = analyse(set_all_childs[[split]][[1]], type_outcome, level_control, D, alpha, upper_best)
          	    pvalues[[split]] = c(pvalues[[split]], analyse_temp[2])
                z_stats_eff[[split]] = c(z_stats_eff[[split]], analyse_temp[1])
                sizes[[split]] = c(sizes[[split]], ifelse(is.null(nrow(set_all_childs[[split]][[1]])), 1, nrow(set_all_childs[[split]][[1]])))
                analyse_temp = analyse(set_all_childs[[split]][[2]], type_outcome, level_control, D, alpha, upper_best)
                pvalues[[split]] = c(pvalues[[split]], analyse_temp[2])
                z_stats_eff[[split]] = c(z_stats_eff[[split]], analyse_temp[1])
                sizes[[split]] = c(sizes[[split]], ifelse(is.null(nrow(set_all_childs[[split]][[2]])), 1, nrow(set_all_childs[[split]][[2]])))
      	    }

      	    # If the best of the two child have a sample size < S, split is eliminated
            comb_child_temp = comb_child
            covar_comb_child_temp = rep(covar, nb_splits)
            z_stats_eff_temp = z_stats_eff
            pvalues_temp = pvalues
            sizes_temp = sizes
            nb_rem=0
      	    for(split in 1:nb_splits) {
                if( ( (pvalues[[split]][1] >  pvalues[[split]][2]) && (sizes[[split]][2] < S) )  || 
                    ( (pvalues[[split]][1] <= pvalues[[split]][2]) && (sizes[[split]][1] < S) ) ){
                    z_stats_eff_temp[[split-nb_rem]] = NULL
                    pvalues_temp[[split-nb_rem]] = NULL
                    sizes_temp[[split-nb_rem]] = NULL
                    comb_child_temp[[split-nb_rem]] = NULL
                    covar_comb_child_temp = covar_comb_child_temp[-(split-nb_rem)]
                    nb_rem = nb_rem+1
          	    }  
                else{
                    vec_z1 = c(vec_z1, z_stats_eff[[split]][1])
           		      vec_z2 = c(vec_z2, z_stats_eff[[split]][2]) 
                    vec_type = c(vec_type, type_var[covar-2])
                    if(type_var[covar-2] == "nominal" || type_var[covar-2] == "ordinal"){
                        vec_levels = c(vec_levels, length(unique(set_subgroup_parent[,covar])))
                    }
                    else{
                        vec_levels = c(vec_levels, ord.cut(set_subgroup_parent[,covar],ord.bin)$g+1)
                    }
                }
       	    }
            # Select the best split (per covariate) with a maximum of M
            best_splits = best_child(vec_z1, vec_z2, num_crit, 1, vec_levels, vec_type, modified)
            ind_best_splits = best_splits[[1]]
            splitting_best_splits = best_splits[[2]]
            if(length(ind_best_splits) > 0){
                subgroup_cur = subgroup_parent
                subgroup_cur[[1]] = c(subgroup_cur[[1]], covar_comb_child_temp[ind_best_splits])
                subgroup_cur[[3]] = c(subgroup_cur[[3]], type_var[covar_comb_child_temp[ind_best_splits]-2])
                if(pvalues_temp[[ind_best_splits]][1] > pvalues_temp[[ind_best_splits]][2]){
          	        subgroup_cur[[2]][[length(subgroup_cur[[2]])+1]] = comb_child_temp[[ind_best_splits]][[2]]
                    if(length(subgroups_child) < M){
                      subgroups_child[[length(subgroups_child)+1]] = subgroup_cur
                      pvalues_childs = c(pvalues_childs, pvalues_temp[[ind_best_splits]][2])
                      splitting_crit_childs = c(splitting_crit_childs, splitting_best_splits)
                    }
                    else{
                        ind_worse = which.max(splitting_crit_childs)
                        subgroups_child[[ind_worse]] = subgroup_cur
                        pvalues_childs[ind_worse] = pvalues_temp[[ind_best_splits]][2]
                        splitting_crit_childs[ind_worse] = splitting_best_splits
                    }
                } 
          	    else {
                    subgroup_cur[[2]][[length(subgroup_cur[[2]])+1]] = comb_child_temp[[ind_best_splits]][[1]]
                    if(length(subgroups_child) < M){
                        subgroups_child[[length(subgroups_child)+1]] = subgroup_cur
                        pvalues_childs = c(pvalues_childs, pvalues_temp[[ind_best_splits]][1])
                        splitting_crit_childs = c(splitting_crit_childs, splitting_best_splits)
                    }
                    else{
                        ind_worse = which.max(splitting_crit_childs)
                        subgroups_child[[ind_worse]] = subgroup_cur
                        pvalues_childs[ind_worse] = pvalues_temp[[ind_best_splits]][1]
                        splitting_crit_childs[ind_worse] = splitting_best_splits
                    }
          	    }
            }

        }
        
        if(length(pvalues_childs) > 0){
            # Continuation criterion for a given gamma
            pvalue_parent = analyse(set_subgroup_parent, type_outcome, level_control, D, alpha, upper_best)[2]
            nb_childs = length(subgroups_child)
            promising_cur = list()
            pval_promising = c()
            for(i in 1:nb_childs){
                 promising_bool = continuation_met(pvalue_parent, pvalues_childs[i], gamma[depth_tree+1])
                 if(promising_bool==TRUE){
                     promising_cur = c(promising_cur, list(subgroups_child[[i]]))
                     pval_promising = c(pval_promising, pvalues_childs[i])
                 }
            }
            if(length(pval_promising) > 0){
                # Add the current promisings subgroups to the final promisings
                promisings <<- c(promisings, promising_cur)
                pvalues_promisings <<- c(pvalues_promisings, pval_promising)
            }
            best_pval <<- min(best_pval, pvalues_childs)

            # For each promising child call the function recursively
            for(child in promising_cur) {
                rec_search(child)
            }
        }
    }
    rec_search(list(c(),list(),c()))
    if(bool_best_pval == FALSE){
        return(list(promisings, pvalues_promisings))
    }
    else{
        return(list(promisings, pvalues_promisings, best_pval))
    }
}




#### Function that identify promisings subgroups (select the M best split)
subgroup_identification_promising2 = function(training_set, type_var, type_outcome, level_control, D=0, alpha, 
L=3, S, num_crit, M=5, gamma, bool_best_pval=FALSE, ord.bin, upper_best=TRUE, modified=TRUE){
    promisings = list()
    pvalues_promisings = c()
    best_pval = +Inf

    # Recursive function on parent
    rec_search = function(subgroup_parent){
        depth_tree = length(subgroup_parent[[1]])
        # stopping rule
        if(depth_tree >= L) {
            return(0)
        }
        # Remaining covariates and subset from parents
        covariates_not_used = setdiff(3:ncol(training_set), subgroup_parent[[1]]) 
        set_subgroup_parent = sub_sets_parents(training_set, subgroup_parent)[[1]]

        # Update all possible leaves of tree from parents and remaining covariates
	      vec_z1 = c()
        vec_z2 = c() 
        vec_type = c()
        vec_levels = c()
	      covar_all_comb_child = c()
        all_comb_child = list()
        all_z_stats_eff = list()
        #all_z_stats_tox = list()
        all_pvalues = list()
        all_sizes = list() 
        for(covar in covariates_not_used){
            # Calculate all combinations of two childs from covariate with pvalues and sample sizes of the corresponding sets
           comb_child = comb_child(set_subgroup_parent[,covar], type_var[covar-2], ord.bin)
            if(anyNA(comb_child)){
                break
            }
            set_all_childs = sub_sets_all_childs(training_set, subgroup_parent, covar, comb_child, type_var[covar-2])
      	    nb_splits = length(comb_child)
            z_stats_eff = vector("list", nb_splits)
      	    pvalues = vector("list", nb_splits)
            sizes = vector("list", nb_splits)
      	    for(split in 1:nb_splits){ 
                analyse_temp = analyse(set_all_childs[[split]][[1]], type_outcome, level_control, D, alpha, upper_best)
          	    pvalues[[split]] = c(pvalues[[split]], analyse_temp[2])
                z_stats_eff[[split]] = c(z_stats_eff[[split]], analyse_temp[1])
                sizes[[split]] = c(sizes[[split]], ifelse(is.null(nrow(set_all_childs[[split]][[1]])), 1, nrow(set_all_childs[[split]][[1]])))
                analyse_temp = analyse(set_all_childs[[split]][[2]], type_outcome, level_control, D, alpha, upper_best)
                pvalues[[split]] = c(pvalues[[split]], analyse_temp[2])
                z_stats_eff[[split]] = c(z_stats_eff[[split]], analyse_temp[1])
                sizes[[split]] = c(sizes[[split]], ifelse(is.null(nrow(set_all_childs[[split]][[2]])), 1, nrow(set_all_childs[[split]][[2]])))
      	    }

      	    # If the best of the two child have a sample size < S, split is eliminated
            nb_values = length(all_comb_child) 
            all_comb_child = c(all_comb_child, comb_child)
            covar_all_comb_child = c(covar_all_comb_child, rep(covar, nb_splits))
            all_z_stats_eff = c(all_z_stats_eff, z_stats_eff)
            all_pvalues = c(all_pvalues, pvalues)
            all_sizes = c(all_sizes, sizes)
            nb_rem=0
      	    for(split in 1:nb_splits) {
                if( ( (pvalues[[split]][1] >  pvalues[[split]][2]) && (sizes[[split]][2] < S) )  || 
                    ( (pvalues[[split]][1] <= pvalues[[split]][2]) && (sizes[[split]][1] < S) ) ){
                    all_z_stats_eff[[nb_values+split-nb_rem]] = NULL
                    all_pvalues[[nb_values+split-nb_rem]] = NULL
                    all_sizes[[nb_values+split-nb_rem]] = NULL
                    all_comb_child[[nb_values+split-nb_rem]] = NULL
                    covar_all_comb_child = covar_all_comb_child[-(nb_values+split-nb_rem)]
                    nb_rem = nb_rem+1
          	    }  
                else{
                    vec_z1 = c(vec_z1, z_stats_eff[[split]][1])
           		      vec_z2 = c(vec_z2, z_stats_eff[[split]][2]) 
                    vec_type = c(vec_type, type_var[covar-2])
                    if(type_var[covar-2] == "nominal" || type_var[covar-2] == "ordinal"){
                        vec_levels = c(vec_levels, length(unique(set_subgroup_parent[,covar])))
                    }
                    else{
                        vec_levels = c(vec_levels, ord.cut(set_subgroup_parent[,covar],ord.bin)$g+1)
                    }
                }
       	    }
        }	    
        # Caculate splitting criteria with correction for all childs and select the best M childs
        best_splits = best_child(vec_z1, vec_z2, num_crit, M, vec_levels, vec_type, modified)
        ind_best_splits = best_splits[[1]]
        splitting_best_splits = best_splits[[2]]
        if(length(ind_best_splits) > 0){
            subgroups_child = list()
            pvalues_childs = c()
            for(s in 1:length(ind_best_splits)) {
                subgroup_cur = subgroup_parent
                subgroup_cur[[1]] = c(subgroup_cur[[1]], covar_all_comb_child[ind_best_splits[s]])
                subgroup_cur[[3]] = c(subgroup_cur[[3]], type_var[covar_all_comb_child[ind_best_splits[s]]-2])
                if(all_pvalues[[ind_best_splits[s]]][1] > all_pvalues[[ind_best_splits[s]]][2]){
          	        subgroup_cur[[2]][[length(subgroup_cur[[2]])+1]] = all_comb_child[[ind_best_splits[s]]][[2]]
                    subgroups_child[[length(subgroups_child)+1]] = subgroup_cur
                    pvalues_childs = c(pvalues_childs, all_pvalues[[ind_best_splits[s]]][2])
                } 
          	    else {
                    subgroup_cur[[2]][[length(subgroup_cur[[2]])+1]] = all_comb_child[[ind_best_splits[s]]][[1]]
                    subgroups_child[[length(subgroups_child)+1]] = subgroup_cur
                    pvalues_childs = c(pvalues_childs, all_pvalues[[ind_best_splits[s]]][1])
          	    }
            }
            if(length(pvalues_childs) > 0){
                # Continuation criterion for a given gamma
                pvalue_parent = analyse(set_subgroup_parent, type_outcome, level_control, D, alpha, upper_best)[2]
                nb_childs = length(subgroups_child)
                promising_cur = list()
                pval_promising = c()
                for(i in 1:nb_childs){
                    promising_bool = continuation_met(pvalue_parent, pvalues_childs[i], gamma[depth_tree+1])
                    if(promising_bool==TRUE){
                        promising_cur = c(promising_cur, list(subgroups_child[[i]]))
                        pval_promising = c(pval_promising, pvalues_childs[i])
                    }
                }
                if(length(pval_promising) > 0){
                    # Add the current promisings subgroups to the final promisings
                    promisings <<- c(promisings, promising_cur)
                    pvalues_promisings <<- c(pvalues_promisings, pval_promising)
                }
                best_pval <<- min(best_pval, pvalues_childs)
                # For each promising child call the function recursively
                for(child in promising_cur) {
                    rec_search(child)
                }
            }
        }
    }
    rec_search(list(c(),list(),c()))
    if(bool_best_pval == FALSE){
        return(list(promisings, pvalues_promisings))
    }
    else{
        return(list(promisings, pvalues_promisings, best_pval))
    }
}





#### Function that identify candidates subgroups
subgroup_identification_candidates = function(training_set, type_var, type_outcome, level_control, 
D=0, L=3, S, num_crit, M=5, gamma, alpha, nsim, ord.bin, upper_best=TRUE, M_per_covar=FALSE, seed=42, modified=TRUE){  
    candidates = list()
    pvalues_candidates = c()
    adj_pvalues_candidates = c()

    # Identification of promising subgroups
    if(M_per_covar==TRUE){
        res_prom = subgroup_identification_promising(training_set, type_var, type_outcome, level_control, D, alpha, L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best, modified)
    }
    else{

        res_prom = subgroup_identification_promising2(training_set, type_var, type_outcome, level_control, D, alpha, L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best, modified)
    }
    promising_subgroups = res_prom[[1]]
    pvalues_promising = res_prom[[2]]
    nb_prom = length(pvalues_promising)

    if(nb_prom >= 1){
        # Adjusted pvalue for significance level for selection criterion
        adjusted_pval_level = adjusted_pval_level(training_set, res_prom, nsim, type_var, type_outcome, level_control, D, L, S, num_crit, M, gamma, alpha, ord.bin, upper_best, M_per_covar, seed)
        # Selection criterion
        for(i in 1:nb_prom){
            if(adjusted_pval_level[i] < alpha){
                candidates = c(candidates, list(promising_subgroups[[i]]))
                pvalues_candidates = c(pvalues_candidates, pvalues_promising[i])
                adj_pvalues_candidates = c(adj_pvalues_candidates, adjusted_pval_level[i])
            }
        }
    }
    return(list(candidates, pvalues_candidates, adj_pvalues_candidates))
}



subgroup_identification_best_cv = function(training_set, type_var, type_outcome, level_control, 
D=0, L=3, S, num_crit, M=5, gamma, alpha, nsim, ord.bin, upper_best=TRUE, M_per_covar=FALSE, seed=42, modified=TRUE){
    candidates = list()
    pvalues_candidates = c()
    adj_pvalues_candidates = c()
    best_pval = 1
    # Identification of promising subgroups
    if(M_per_covar==TRUE){
        res_prom = subgroup_identification_promising(training_set, type_var, type_outcome, level_control, D, alpha, L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best, modified)
    }
    else{
        res_prom = subgroup_identification_promising2(training_set, type_var, type_outcome, level_control, D, alpha, L, S, num_crit, M, gamma, FALSE, ord.bin, upper_best, modified)
    }
    promising_subgroups = res_prom[[1]]
    pvalues_promising = res_prom[[2]]
    nb_prom = length(pvalues_promising)

    # Adjusted pvalue for significance level for selection criterion
    adjusted_pval_level = adjusted_pval_level(training_set, res_prom, nsim, type_var, type_outcome, level_control, D, L, S, num_crit, M, gamma, alpha, ord.bin, upper_best, M_per_covar, seed)

    # Selection criterion
    for(i in 1:nb_prom){
        if(adjusted_pval_level[i] < best_pval){
            candidates = promising_subgroups[[i]]
            pvalues_candidates = pvalues_promising[i]
            adj_pvalues_candidates = adjusted_pval_level[i]
            best_pval = adj_pvalues_candidates
        }
    }
    return(list(candidates, pvalues_candidates, adj_pvalues_candidates))
}





