library(nnet)

allocation_procedure = function(H, pct_random, Xcov, type_var, prop_gpe, alloc_hp=TRUE, overall_imb=FALSE, seed=NA){
    if(is.na(seed)==FALSE){
        set.seed(seed)
    }
    nb_patients = nrow(Xcov)
    if(H > 1){
        X = Xcov
        ind_cont = which(type_var=="continuous")
        nb_cont = length(ind_cont)
        if(nb_cont>0){
            for(j in 1:nb_cont){
                quant = quantile(Xcov[,ind_cont[j]], probs = seq(0, 1, 0.33))
                q33 = as.numeric(quant[2])
                q66 = as.numeric(quant[3])
                for(i in 1:nb_patients){
                    if(Xcov[i,ind_cont[j]] < q33){
                        X[i,ind_cont[j]] = 0
                    }
                    else if(Xcov[i,ind_cont[j]] >= q66){
                        X[i,ind_cont[j]] = 2
                    }
                    else{
                        X[i,ind_cont[j]] = 1
                    }
                }
            }
        }
    
        nb_pat_max_gpe = round(prop_gpe*nb_patients)
        J = ncol(X)
        set_alloc = rep(NA,nb_patients)
        nb_pat_random = max(round(nb_patients*pct_random),1)  
        nb_pat_remain = nb_patients-nb_pat_random

        # Allocation of "pct_random" of the sample size randomly between the H sets
        sets = 1:H
        set_full = numeric(0)
        for(s in 1:nb_pat_random){
            set_cur = 1:H
            alloc_cur = table(set_alloc)
            set_allocated = as.numeric(names(alloc_cur))
            pat_allocated = numeric(H)
            pat_allocated[set_allocated] = as.numeric(alloc_cur[paste(set_allocated)])
            prob_alloc = prop_gpe
          	if(length(set_allocated)>0){
                set_full = which(pat_allocated >= nb_pat_max_gpe)
                if(length(set_full)>0){
                    set_cur = sets[-set_full]
                    prob_alloc = prop_gpe[-set_full]
                }
            }
            if(length(set_cur)>1){
                set_alloc[s] = sample(set_cur, 1, replace=TRUE, prob=prob_alloc)
            }
        	  else{
                set_alloc[s] = set_cur
            }
        }

        # Allocation based on imbalanced score
    	  if(nb_pat_remain>0){
            for(s in (nb_pat_random+1):nb_patients){
                covariates_cur = as.numeric(X[s,])
                alloc_cur = table(set_alloc)
                set_allocated = as.numeric(names(alloc_cur))
                set_cur = 1:H
                set_full = numeric(0)
                pat_allocated = numeric(H)
                pat_allocated[set_allocated] = as.numeric(alloc_cur[paste(set_allocated)])
                if(length(set_allocated)>0){
                    set_full = which(pat_allocated >= nb_pat_max_gpe)
                }
                if(length(set_full) == H){
                    set_full = which(pat_allocated >= nb_pat_max_gpe+1)
                }
                if(length(set_full)>0){
                    set_cur = set_cur[-set_full]
                }
                H_cur = length(set_cur)

                if(H_cur > 1){
                    f_hij = array(0, dim=c(H_cur,H_cur,J))
                    for(i in 1:H_cur){
                        n_i_cur = length(which(set_alloc==set_cur[i]))
                        if(n_i_cur > 0){
                            for(j in 1:J){
                                n_ij_cur = length(which(X[,j]==covariates_cur[j] & set_alloc==set_cur[i]))
                                f_hij[,i,j] = (n_ij_cur + (set_cur[i]==set_cur)) / (n_i_cur + (set_cur[i]==set_cur))
                            }
                        }
                        else if(n_i_cur == 0){
                            for(h in 1:H_cur){
                                f_hij[h,h,] = rep(1,J)
                            }
                        }
                    }   
                    d_jh_cur = matrix(NA, nrow=H_cur, ncol=J)
                    for(h in 1:H_cur){
                        for(j in 1:J){
                            d_jh_cur[h,j] = max(f_hij[h,,j])-min(f_hij[h,,j])
                        }
                    }
                    d_h = rowSums(d_jh_cur)
                    d_tot = sum(d_h)
                    f_h = 1/(H_cur-1) * (1 - d_h/d_tot)
                    if(alloc_hp==TRUE){     
                        set_alloc[s] = set_cur[which.is.max(f_h)]
                    }
                    else{
                        set_alloc[s] = sample(set_cur, 1, replace=TRUE, prob=f_h)
                    }
                }
                else{
                   set_alloc[s] = set_cur
                }
            }
        }
        # Overall imbalanced score
        if(overall_imb){
            nb_levels_cov = numeric(J)
            for(j in 1:J){
                nb_levels_cov[j] = length(unique(X[,j]))
            }
            f_end_ijl = array(NA, dim=c(H, J, max(nb_levels_cov)))
            for(i in 1:H){
                for(j in 1:J){
                    levels_cov = unique(X[,j])
                    for(l in 1:nb_levels_cov[j]){
                        f_end_ijl[i,j,l] = length(which(X[,j]==levels_cov[l] & set_alloc==i))/length(which(set_alloc==i))
                    }
                }
            }
            djl = matrix(NA, nrow=J, ncol=max(nb_levels_cov))
            for(j in 1:J){
                for(l in 1:nb_levels_cov[j]){
                    djl[j,l] = max(f_end_ijl[,j,l],na.rm=TRUE)-min(f_end_ijl[,j,l],na.rm=TRUE)
                }
            }
            d_j = numeric(J)
            for(j in 1:J){
                d_j[j] = max(djl[j,],na.rm=TRUE)
            }
            overall_score = 1/J*sum(d_j)
            return(list(set_alloc,overall_score))
        }
        else{
            return(set_alloc)
        }
    }
    else{
        return(rep(1, nb_patients))
    }
}


