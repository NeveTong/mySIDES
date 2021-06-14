###################
# SIDES ALGORITHM #
###################

catch_entries_commun = function(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                 num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed){
    if(ncol(all_set)<=2){
        stop("Data base (all_set) does not contain any covariate.")
    }
    if(type_outcome != "binary" && type_outcome != "continuous" && type_outcome != "survival" && type_outcome != "count"){
        stop("Outcome (type_outcome) should be either continuous, binary or survival.")
    }
    if(sum((type_var != "continuous") & (type_var != "ordinal") & (type_var != "nominal")) > 0){
        stop("Covariate types (type_var) must be either continuous, ordinal, or nominal.")
    }
    if(as.integer(L) != L){
        L = as.integer(L)
        print(paste("The maximum number of covariates to define subgroups L was transformed into integer: L=", L, sep=""))
    }
    if(L < 1){
        stop("The maximum number of covariates to define subgroups L must be superior or equal to 1.")
    }
    if(as.integer(M) != M){
        M = as.integer(M)
        print(paste("The maximum number of best promising subgroups M was transformed into integer: M=", M, sep=""))
    }
    if(M < 1){
        stop("The maximum number of best promising subgroups M must be superior or equal to 1.")
    }
    if(length(gamma)==1 && is.na(gamma)){
        print("Vector of relative improvment gamma was not supplied and will be chosen by cross-validation. This is time-consuming and not recommended.")
    }
    if(length(gamma) != L){
        stop("Lenght of vector of relative improvment gamma should be equal to L.")
    }
    if(as.integer(H) != H){
        H = as.integer(H)
        print(paste("The number of sets H was transformed into integer: H=", H, sep=""))
    }
    if(H < 1){
        stop("The number of sets H must be superior or equal to 1.")
    }
    if(pct_rand < 0 || pct_rand > 1){
        stop("Percentage of sample size allocated randomly between sets (pct_rand) must be comprised between 0 and 1.")
    }
    if(length(prop_gpe) != H){
        stop("Vector of proportions in each set (prop_gpe) must be equal to the number of sets H.")
    }
    if(sum(prop_gpe<0)>0 || sum(prop_gpe> 1)>0){
        stop("Proportions of sample size in each set (prop_gpe) must be comprised between 0 and 1.")
    }
    if(alloc_high_prob != TRUE && alloc_high_prob != FALSE){
        alloc_high_prob = TRUE
        print("The allocation rule (alloc_high_prob) was misspecified and was thus reset to its default value.")
    }
    if(num_crit != 1 && num_crit != 2 && num_crit != 3){
        num_crit = 1
        print("The number associated to splitting criterion (num_crit) was misspecified and was thus reset to its default value.")
    }
    if(is.na(gamma) && (step < 0 || step > 1)){
        stop("step for cross-validation must be comprised between 0 and 1.")
    }
    if(is.na(gamma) && as.integer(nb_sub_cross) != nb_sub_cross){
        nb_sub_cross = as.integer(nb_sub_cross)
        print(paste("The number of folds for cross-validation (nb_sub_cross) was transformed into integer: nb_sub_cross=", nb_sub_cross, sep=""))
    }
    if(nb_sub_cross < 2){
        stop("The number of folds for cross-validation (nb_sub_cross) must be superior or equal to 2.")
    }
    if(alpha < 0 || alpha > 1){
        stop("Type I error rate (alpha) must be comprised between 0 and 1.")
    }
    if(as.integer(nsim) != nsim){
        nsim = as.integer(nsim)
        print(paste("The number of permutations for resampling-based methods to adjust pvalues (nsim) was transformed into integer: nsim=", nsim, sep=""))
    }
    if(nsim < 0){
        stop("The number of permutations for resampling-based methods to adjust pvalues (nsim) must be superior or equal to 0.")
    }
    if(!is.na(gamma) && as.integer(nsim_cv) != nsim_cv){
        nsim_cv = as.integer(nsim_cv)
        print(paste("The number of permutations for resampling-based methods to adjust pvalues in the cross-validation part (nsim_cv) was transformed into integer: nsim_cv=", nsim_cv, sep=""))
    }
    if(!is.na(gamma) && nsim_cv < 0){
        stop("The number of permutations for resampling-based methods to adjust pvalues in the cross-validation part (nsim_cv) must be superior or equal to 0.")
    }
    if(as.integer(ord.bin) != ord.bin){
        ord.bin = as.integer(ord.bin)
        print(paste("The number of classes to discretize covariates (ord.bin) was transformed into integer: ord.bin=", ord.bin, sep=""))
    }
    if(ord.bin < 2){
        stop("The number of classes to discretize covariates (ord.bin) must be superior or equal to 2.")
    }
    if(M_per_covar != TRUE && M_per_covar != FALSE){
        M_per_covar = FALSE
        print("The selection rule for best promising child subgroups (M_per_covar) was misspecified and was thus reset to its default value.")
    }
    if(upper_best != TRUE && upper_best != FALSE){
        upper_best = TRUE
        print("Boolean indicating if greater values of the outcome mean better responses (upper_best) was misspecified and was thus reset to its default value.")
    }
    if(as.integer(seed) != seed){
        seed = as.integer(seed)
        print(paste("The seed was transformed into integer: seed=", seed, sep=""))
    }
    return(list(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
           num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed))
}

catch_entries1 = function(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                 num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, selec){
    catch = catch_entries_commun(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                         num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed)
    all_set = catch[[1]]
    type_var = catch[[2]]
    type_outcome = catch[[3]]
    level_control = catch[[4]]
    D = catch[[5]]
    L = catch[[6]]
    S = catch[[7]]
    M = catch[[8]]
    gamma = catch[[9]]
    H = catch[[10]]
    pct_rand = catch[[11]]
    prop_gpe = catch[[12]]
    alloc_high_prob = catch[[13]]
    num_crit = catch[[14]]
    step = catch[[15]]
    nb_sub_cross = catch[[16]]
    alpha = catch[[17]]
    nsim = catch[[18]]
    nsim_cv = catch[[19]]
    ord.bin = catch[[20]]
    M_per_covar = catch[[21]]
    upper_best = catch[[22]]
    seed = catch[[23]]
    if(selec != TRUE && selec != FALSE){
        selec = FALSE
        print("Boolean indicating if the function also print subgroups selected and not necessarily validated (selec) was misspecified and was thus reset to its default value.")
    }
    return(list(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
         num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, selec))
}

catch_entries2 = function(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                 num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, nrep){
    catch = catch_entries_commun(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                         num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed)
    all_set = catch[[1]]
    type_var = catch[[2]]
    type_outcome = catch[[3]]
    level_control = catch[[4]]
    D = catch[[5]]
    L = catch[[6]]
    S = catch[[7]]
    M = catch[[8]]
    gamma = catch[[9]]
    H = catch[[10]]
    pct_rand = catch[[11]]
    prop_gpe = catch[[12]]
    alloc_high_prob = catch[[13]]
    num_crit = catch[[14]]
    step = catch[[15]]
    nb_sub_cross = catch[[16]]
    alpha = catch[[17]]
    nsim = catch[[18]]
    nsim_cv = catch[[19]]
    ord.bin = catch[[20]]
    M_per_covar = catch[[21]]
    upper_best = catch[[22]]
    seed = catch[[23]]
    if(as.integer(nrep) != nrep){
        nrep = as.integer(nrep)
        print(paste("The number of simulations (nrep) was transformed into integer: nrep=", nrep, sep=""))
    }
    if(nrep < 1){
        stop("The number of simulations (nrep) must be superior or equal to 1.")
    }
    return(list(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, nrep))
}


#### SIDES algorithm
SIDES = function(all_set, type_var, type_outcome, level_control, D=0, L=3, S, M=5, gamma=rep(1,3), H=1, pct_rand=0.5, prop_gpe=c(1), alloc_high_prob=TRUE, 
                 num_crit=1, step=0.5, nb_sub_cross=5, alpha, nsim=500, nsim_cv=500, ord.bin=10, M_per_covar=FALSE, 
                 upper_best=TRUE, selec=FALSE, seed=42, modified=TRUE){  
                 
    catch = catch_entries1(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                 num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, selec)    
    all_set = catch[[1]]
    type_var = catch[[2]]
    type_outcome = catch[[3]]
    level_control = catch[[4]]
    D = catch[[5]]
    L = catch[[6]]
    S = catch[[7]]
    M = catch[[8]]
    gamma = catch[[9]]
    H = catch[[10]]
    pct_rand = catch[[11]]
    prop_gpe = catch[[12]]
    alloc_high_prob = catch[[13]]
    num_crit = catch[[14]]
    step = catch[[15]]
    nb_sub_cross = catch[[16]]
    alpha = catch[[17]]
    nsim = catch[[18]]
    nsim_cv = catch[[19]]
    ord.bin = catch[[20]]
    M_per_covar = catch[[21]]
    upper_best = catch[[22]]
    seed = catch[[23]]
    selec = catch[[24]]
    X_covariate = all_set[,-c(1,2)]
    # Balanced allocation procedure
    alloc_btw_sets = allocation_procedure(H, pct_rand, X_covariate, type_var, prop_gpe, alloc_high_prob, FALSE, seed)
    base = cbind(alloc_btw_sets, all_set)
    
    # Training set
    training_set = all_set[which(base[,1]==1),]
    # Cross validation to determine gamma
    if(is.na(gamma[1])==TRUE){
        gamma = cross_validation(training_set, type_var, type_outcome, level_control, D, alpha, L, S, num_crit, 
                                  M, step, nb_sub_cross, nsim_cv, ord.bin, upper_best, M_per_covar, seed)
        if(is.null(nrow(gamma))==FALSE){
            gamma = as.numeric(gamma[1,])
        }
    }
    # Candidates subgroups
    res_candidates = subgroup_identification_candidates(training_set, type_var, type_outcome, level_control, D, L, S, num_crit, M, gamma, alpha, nsim, ord.bin, upper_best, M_per_covar, seed, modified)
    candidates = res_candidates[[1]]
    nb_candidates = length(candidates)  
    if(nb_candidates==0){
        print("No subgroup identified")
        res = list("candidates"=list(list(),c()), "confirmed"=list(list(),c()))
    }
    else{
        # Suppress doublons
        if(nb_candidates > 1){
            to_comp = combn(1:nb_candidates,2)
            candidates_temp = candidates
            pval_cand_temp = res_candidates[[3]]
            ind_rem = c()
            for(icol in 1:ncol(to_comp)){
                if(identical_subgroups(candidates[[to_comp[1,icol]]], candidates[[to_comp[2,icol]]])==TRUE){
                    ind_rem = c(ind_rem, to_comp[2,icol])
                }
            }
            ind_rem = sort(unique(ind_rem), decreasing=TRUE)
            for(ir in ind_rem){
                candidates_temp[[ir]] = NULL
            }
            candidates = candidates_temp
            nb_candidates = length(candidates)
            if(length(ind_rem)>0){
                pval_cand = pval_cand_temp[-ind_rem]
            }
            else{
                pval_cand = pval_cand_temp
            }
        }
        else{
          pval_cand = res_candidates[[3]]
        }
        # Validation on other sets
        if(H > 1){
            confirmed = list()
            pval_confirmed = c()
            for(cand in 1:nb_candidates){
                conf_all_set = TRUE
                i=2
                while(i<=H && conf_all_set==TRUE){
                    set_valid_cur = all_set[which(base[,1]==i),]
                    set_subg_cur = sub_sets_parents(set_valid_cur, candidates[[cand]])[[1]]
                    res_analyse = analyse(set_subg_cur, type_outcome, level_control, D, alpha, upper_best)
                    if(res_analyse[3]==FALSE){
                        conf_all_set = FALSE
                    }
                    else{
                        if(i == H){
                            confirmed[[length(confirmed)+1]] = candidates[[cand]]
                            pval_confirmed = c(pval_confirmed,res_analyse[2])
                        }
                    }
                    i=i+1
                }
            }
            if(length(pval_confirmed) > 0){
                if(selec==FALSE){ 
                    res = list("candidates"=list(list(),c()), "confirmed"=list(confirmed,pval_confirmed))
                }
                else{
                    res = list("candidates"=list(candidates,pval_cand),"confirmed"=list(confirmed,pval_confirmed))
                }
            }
            else{
                print("No subgroup confirmed")
                res = list("candidates"=list(list(),c()), "confirmed"=list(list(),c()))
            }
        }
        else{
            res = list("candidates"=list(candidates,pval_cand),"confirmed"=list(list(),c()))
        }  
    }
    res = c(res,"base"=list(all_set),"training"=list(training_set))
    class(res) = "SIDES"
    return(res)
}


#### Simulations on SIDES
simulation_SIDES = function(all_set, type_var, type_outcome, level_control, D=0, L=3, S, M=5, num_crit=1, gamma=rep(1,3), 
                            alpha, nsim=500, ord.bin=10, nrep=100, seed=42, 
                            H=1, pct_rand=0.5, prop_gpe=c(1), alloc_high_prob=TRUE, 
                            step=0.5, nb_sub_cross=5, nsim_cv=500,
                            M_per_covar=FALSE, upper_best=TRUE, nb_cores=NA, ideal=NA, modified=TRUE){
    catch = catch_entries2(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                 num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, seed, nrep)
  
    all_set = catch[[1]]
    type_var = catch[[2]]
    type_outcome = catch[[3]]
    level_control = catch[[4]]
    D = catch[[5]]
    L = catch[[6]]
    S = catch[[7]]
    M = catch[[8]]
    gamma = catch[[9]]
    H = catch[[10]]
    pct_rand = catch[[11]]
    prop_gpe = catch[[12]]
    alloc_high_prob = catch[[13]]
    num_crit = catch[[14]]
    step = catch[[15]]
    nb_sub_cross = catch[[16]]
    alpha = catch[[17]]
    nsim = catch[[18]]
    nsim_cv = catch[[19]]
    ord.bin = catch[[20]]
    M_per_covar = catch[[21]]
    upper_best = catch[[22]]
    seed = catch[[23]]
    nrep = catch[[24]]  
    if(is.na(nb_cores)){
        nb_cores = detectCores()-1
    }   
    if(nb_cores>1){
        cl = makeCluster(nb_cores, outfile="")
        registerDoParallel(cl)
    }
  
    if(H==1){
        n_rep = 1
    }
                 
    list_selected = list()
    list_top = list()
    pct_selected = c()
    pct_top = c()
    pct_no_subgroup = 0
    pct_sous_cov_select1 = 0
    pct_sous_ens_top1 = 0
    pct_sous_cov_select2 = 0
    pct_sous_ens_top2 = 0
    pct_ideal_selected = 0
    pct_ideal_top = 0
    mean_size = 0

    #Simulate nrep replications of analysis
    if(nb_cores>1){
      res_simu = foreach(r=1:nrep, .export=ls(globalenv()), .inorder=FALSE) %dopar% {
          set.seed(seed+r)
print(r)
          res_r = SIDES(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                   num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, selec=FALSE, seed+r, modified)    
          return(res_r)
      }
    }
    else{
      res_simu = list()
      for(r in 1:nrep){
        set.seed(seed+r)
print(r)
        res_r = SIDES(all_set, type_var, type_outcome, level_control, D, L, S, M, gamma, H, pct_rand, prop_gpe, alloc_high_prob, 
                             num_crit, step, nb_sub_cross, alpha, nsim, nsim_cv, ord.bin, M_per_covar, upper_best, selec=FALSE, seed+r, modified)    
        res_simu = c(res_simu, list(res_r))
      }
    }
    #Format results
    for(r in 1:nrep){
        res_r = res_simu[[r]]
        if( (H==1 && length(res_r$candidates[[2]]) > 0) || (H>1 && length(res_r$confirmed[[2]]) > 0) ){
            if(H > 1){
                select_cur = res_r$confirmed[[1]]
                pval_cur = res_r$confirmed[[2]]
            }
            else{
                select_cur = res_r$candidates[[1]]
                pval_cur = res_r$candidates[[2]]
            }
            
            find_sous_cov1 = FALSE
            find_sous_ens1 = FALSE
            find_sous_cov2 = FALSE
            find_sous_ens2 = FALSE
           
            if(length(list_selected)==0){
                list_selected = c(list_selected, select_cur)
                pct_selected = c(pct_selected, rep(1, length(select_cur)))
                for(sg in 1:length(select_cur)){
                    cand_sg = select_cur[[sg]]
                    mean_size = mean_size + nrow(sub_sets_parents(res_r$training, cand_sg)[[1]])/length(select_cur)
                    if(identical_subgroups(ideal, cand_sg)==TRUE){
                        pct_ideal_selected = pct_ideal_selected+1
                        if(pval_cur[sg]==min(pval_cur)){
                            pct_ideal_top = pct_ideal_top+1
                        }
                    }
                    if(find_sous_cov1 == FALSE && included_subgroups(cand_sg, ideal)==TRUE){
                        find_sous_cov1 = TRUE
                        pct_sous_cov_select1 = pct_sous_cov_select1+1
                    }
                    if(find_sous_ens1 == FALSE && included_subgroups(ideal, cand_sg)==TRUE){
                        find_sous_ens1 = TRUE
                        pct_sous_ens_top1 = pct_sous_ens_top1+1
                    }    
                    if(find_sous_cov2 == FALSE && included_subgroups(cand_sg, ideal)==TRUE && identical_subgroups(cand_sg, ideal)==FALSE){
                        find_sous_cov2 = TRUE
                        pct_sous_cov_select2 = pct_sous_cov_select2+1
                    }
                    if(find_sous_ens2 == FALSE && included_subgroups(ideal, cand_sg)==TRUE && identical_subgroups(cand_sg, ideal)==FALSE){
                        find_sous_ens2 = TRUE
                        pct_sous_ens_top2 = pct_sous_ens_top2+1
                    }
                }
            }
            else{ 
                if(length(select_cur)>0){ 
                    for(s in 1:length(select_cur)){
                        cand_s = select_cur[[s]]
                        different = TRUE
                        i=1
                        while(different==TRUE && i <= length(list_selected)){
                            if(identical_subgroups(cand_s, list_selected[[i]])==TRUE){
                                different=FALSE
                                pct_selected[i] = pct_selected[i]+1
                            }
                            i = i+1
                        }
                        if(different == TRUE){
                             list_selected = c(list_selected, list(cand_s))
                             pct_selected = c(pct_selected, 1)
                        }
                        mean_size = mean_size + nrow(sub_sets_parents(res_r$training, cand_s)[[1]])/length(select_cur)
                        if(identical_subgroups(ideal, cand_s)==TRUE){
                            pct_ideal_selected = pct_ideal_selected+1
                            if(pval_cur[s]==min(pval_cur)){
                                pct_ideal_top = pct_ideal_top+1
                            }
                        }
                        if(find_sous_cov1 == FALSE && included_subgroups(cand_s, ideal)==TRUE){
                            find_sous_cov1 = TRUE
                            pct_sous_cov_select1 = pct_sous_cov_select1+1
                        }
                        if(find_sous_ens1 == FALSE && included_subgroups(ideal, cand_s)==TRUE){
                            find_sous_ens1 = TRUE
                            pct_sous_ens_top1 = pct_sous_ens_top1+1
                        }            
                        if(find_sous_cov2 == FALSE && included_subgroups(cand_s, ideal)==TRUE && identical_subgroups(cand_s, ideal)==FALSE){
                            find_sous_cov2 = TRUE
                            pct_sous_cov_select2 = pct_sous_cov_select2+1
                        }
                        if(find_sous_ens2 == FALSE && included_subgroups(ideal, cand_s)==TRUE && identical_subgroups(cand_s, ideal)==FALSE){
                            find_sous_ens2 = TRUE
                            pct_sous_ens_top2 = pct_sous_ens_top2+1
                        }
                    }
                }
            }
        }
        else{
            pct_no_subgroup = pct_no_subgroup+1
        } 
    }  
    mean_size = mean_size/(nrep-pct_no_subgroup)
    pct_selected = pct_selected/nrep*100
    pct_no_subgroup = pct_no_subgroup/nrep*100
    or_pct_selected = order(pct_selected, decreasing=TRUE)
    pct_ideal_selected = pct_ideal_selected/nrep*100
    pct_ideal_top = pct_ideal_top/nrep*100
    pct_sous_cov_select1 = pct_sous_cov_select1/nrep*100
    pct_sous_ens_top1 = pct_sous_ens_top1/nrep*100
    pct_sous_cov_select2 = pct_sous_cov_select2/nrep*100
    pct_sous_ens_top2 = pct_sous_ens_top2/nrep*100
    if(nb_cores>1){
        stopCluster(cl)
    }
    
    res = list( "pct_no_subgroup"=pct_no_subgroup, "mean_size"=mean_size,
    "pct_ideal_selected"=pct_ideal_selected, "pct_ideal_top"=pct_ideal_top,
    "pct_sous_cov_select1"=pct_sous_cov_select1, "pct_sous_ens_top1"=pct_sous_ens_top1, 
    "pct_sous_cov_select2"=pct_sous_cov_select2, "pct_sous_ens_top2"=pct_sous_ens_top2,
    "subgroups"=list_selected[or_pct_selected], "pct_selection"=pct_selected[or_pct_selected],
    "ideal"=ideal )  
    res = c(res,"base"=list(all_set))
    class(res) = "simulation_SIDES"
    return(res)
}


#identical_subgroups = function(g1, g2){
#    res = FALSE
#    if(length(g1[[1]]) == length(g2[[1]])){
#        level_identical = 0
#        for(j in 1:length(g1[[1]])){
#            or_g1 = order(g1[[1]])
#            or_g2 = order(g2[[1]])
#            if(g1[[1]][or_g1][j]==g2[[1]][or_g2][j] && length(g1[[2]][[or_g1[j]]])==length(g2[[2]][[or_g2[j]]])){
#                level_identical_temp = 0
#                for(k in 1:length(g1[[2]][[or_g1[j]]])){
#                    if(g1[[2]][[or_g1[j]]][k]==g2[[2]][[or_g2[j]]][k]){
#                        level_identical_temp = level_identical_temp+1
#                    }
#                }
#                if(level_identical_temp == length(g1[[2]][[or_g1[j]]])){
#                    level_identical = level_identical+1
#                }
#            }
#        }
#        if(level_identical == length(g1[[1]])){
#            res = TRUE
#        }
#    }
#    return(res)
#}

identical_subgroups = function(g1, g2){
    res = FALSE
    if(sum(!is.element(g1[[1]], g2[[1]])) == 0 && sum(!is.element(g2[[1]], g1[[1]])) == 0){
        for(j in 1:length(g1[[1]])){
            ind_j = which(g1[[1]][j]==g2[[1]])
            if(sum(!is.element(g1[[2]][[j]], g2[[2]][[ind_j]])) == 0 && sum(!is.element(g2[[2]][[ind_j]], g1[[2]][[j]])) == 0){
                res = TRUE
            }
        }
    }
    return(res)
}


included_subgroups = function(g1, g2){
    res = FALSE
    if(sum(!is.element(g1[[1]], g2[[1]])) == 0){
        for(j in 1:length(g1[[1]])){
            ind_j = which(g1[[1]][j]==g2[[1]])
            if(sum(!is.element(g1[[2]][[j]], g2[[2]][[ind_j]])) == 0){
                res = TRUE
            }
        }
    }
    return(res)
}



#function to print one subgroup with pvalue
print_gpe = function(subgroup, pval=NA, x, pct=NA){ 
    icov = subgroup[[1]]
    nb_cov = length(icov)
    type_var = subgroup[[3]] 
    levels_icov = subgroup[[2]] 
    txt_sgpe = c()
    for(i in 1:nb_cov){
        levels_theo = sort(unique(x$base[,icov[i]]))
        levels_sgpe = c()
        if(type_var[i]=="ordinal"){
            val_cut = as.numeric(substr(levels_icov[[i]],1,nchar(levels_icov[[i]])-1))
            signe = substr(levels_icov[[i]],nchar(levels_icov[[i]]),nchar(levels_icov[[i]]))
            levels_sgpe = ""
            if(signe == "-"){
                levels_sgpe = levels_theo[which(levels_theo<=val_cut)]
            }
            else{
                levels_sgpe = levels_theo[which(levels_theo>val_cut)]
            }
            tlevels_sgpe = paste(levels_sgpe, collapse=",")
            txt_sgpe = c(txt_sgpe, paste(names(x$base)[icov[i]], " = {", tlevels_sgpe,"}",sep="")) 
        }
        else if(type_var[i]=="nominal"){
            levels_sgpe = levels_icov[[i]]
            tlevels_sgpe = paste(levels_sgpe, collapse=",")
            txt_sgpe = c(txt_sgpe, paste(names(x$base)[icov[i]], " = {", tlevels_sgpe,"}",sep="")) 
        }
        else if(type_var[i]=="continuous"){
            val_cut = as.numeric(substr(levels_icov[[i]],1,nchar(levels_icov[[i]])-1))
            signe = substr(levels_icov[[i]],nchar(levels_icov[[i]]),nchar(levels_icov[[i]]))
            levels_sgpe = ""
            if(signe == "-"){
                signe = "<="
            }
            else{
                signe = ">"
            }
            txt_sgpe = c(txt_sgpe, paste(names(x$base)[icov[i]], " ", signe, " ", val_cut, sep=""))
        }
        if(i < nb_cov){
            txt_sgpe = c(txt_sgpe, " AND ")
        }
        else{
            txt_sgpe = c(txt_sgpe, "\n")
        }
    }
    cat(txt_sgpe)
    if(!is.na(pval)){
        cat("pvalue = ", pval, "\n")
    }
    if(!is.na(pct)){
        cat("Percentage of selection = ", pct, "% \n")
    }
}

#gg=list(c(8,3,13),list(c(0,3),"0-","1.52+"),c("nominal","ordinal","continuous"))
#print_gpe(gg,0.0124)
#gg2=list(c(8,3,13),list("1+","0-",c(1,3)),c("ordinal","ordinal","nominal"))
#print_gpe(gg2,0.00058)
    

print.SIDES = function(x, ...){
    nb_cand = length(x$candidates[[2]])
    nb_conf = length(x$confirmed[[2]])
    if(nb_cand>0){
        cat("Identified candidate subgroups before confirmation phase:\n")
        for(i in 1:nb_cand){
            print_gpe(subgroup=x$candidates[[1]][[i]], pval=x$candidates[[2]][i], x=x)
        }
    }
    else{
        cat("No candidate subgroups identified before confirmation phase:\n")
    }
    if(nb_conf>0){
        cat("Confirmed candidate subgroups:\n")
        for(i in 1:nb_conf){
            print_gpe(subgroup=x$confirmed[[1]][[i]], pval=x$confirmed[[2]][i], x=x)
        }
    }
    else{
        cat("No candidate subgroups confirmed:\n")
    }
}





print.simulation_SIDES = function(x, ...){
    nb_ssgpe = length(x$pct_selection)
    others = FALSE
    cat("No subgroup selected in ", x$pct_no_subgroup, "% \n")
    cat("Average size of the confirmed subgroups in the training data set in ", x$mean_size, "\n")
    if(length(x$ideal) > 1 || (length(x$ideal) == 1 && is.na(x$ideal)==FALSE)){      
        cat("Percentage of simulations where the ideal subgroup is confirmed: ", x$pct_ideal_selected, "% \n")
        cat("Percentage of simulations where the ideal subgroup is the top confirmed subgroup: ", x$pct_ideal_top, "% \n")
        cat("Percentage of simulations where a subgroup containing a subset of the covariates used to define the ideal subgroup is selected (including the ideal): ", x$pct_sous_cov_select1, "% \n")
        cat("Percentage of simulations where a subgroup containing a subset of the covariates used to define the ideal subgroup is selected (excluding the ideal): ", x$pct_sous_cov_select2, "% \n")
        cat("Percentage of simulations where a subset of the ideal subgroup is selected (including the ideal): ", x$pct_sous_ens_top1, "% \n")
        cat("Percentage of simulations where a subset of the ideal subgroup is selected (exluding the ideal): ", x$pct_sous_ens_top2, "% \n")
    }
    if(nb_ssgpe>0){
        cat("Confirmed candidate subgroups:\n")
        for(i in 1:nb_ssgpe){
            if(x$pct_selection[i] >= 10){
                print_gpe(subgroup=x$subgroups[[i]], x=x, pct=x$pct_selection[i])
            }
            else{
                others = TRUE
            }
        }
        if(others == TRUE){
            cat("Others subgroups in less than 10% \n")
        }
    }   
}




