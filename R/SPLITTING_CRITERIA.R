######################
# SPLITTING CRITERIA #
######################

#### Criterion 1
criterion1 = function(z1,z2){
    return(2 * (1 - pnorm( abs(z1-z2)/sqrt(2), 0,1)))
}

#### Criterion 2
criterion2 = function(z1,z2){
    return( 2 * min(1-pnorm(z1, 0,1), 1-pnorm(z2, 0,1)) )
}


#### Criterion 3 implemented directly in best_child


#### Calculate splitting criterion
#splitting_crit = function(z1,z2,z3=NA,z4=NA,w=NA,num_crit){
splitting_crit = function(z1,z2,num_crit){
    splitting_criterion = NA
    if(num_crit==1){
        splitting_criterion = criterion1(z1,z2)
    }
    else if(num_crit==2){
        splitting_criterion = criterion2(z1,z2)
    }
    return(splitting_criterion)
}


#### Function that order the pair of child (the splits) from the best to the worst 
#### in terms of splitting criterion, and select the best M pairs
best_child = function(vec_z1, vec_z2, num_crit, M, vec_levels, vec_type, modified){
    nb_splits = length(vec_z1)
    if(nb_splits > 0){
        if(num_crit == 1 || num_crit == 2){
            split_crit = rep(NA, nb_splits)
            for(i in 1:nb_splits){
                split_crit[i] = splitting_crit(vec_z1[i],vec_z2[i],num_crit)
            }
            adj_split_crit = adjusted_pvalue(split_crit, vec_levels, vec_type, num_crit, modified)
        }
        else if(num_crit == 3){
            split_crit = matrix(NA, nrow=nb_splits, ncol=2)
            for(i in 1:nb_splits){
                split_crit[i,1] = splitting_crit(vec_z1[i],vec_z2[i],1)
                split_crit[i,2] = splitting_crit(vec_z1[i],vec_z2[i],2)
            }
            adj_split_crit1 = adjusted_pvalue(split_crit[,1], vec_levels, vec_type, 1, modified) 
            adj_split_crit2 = adjusted_pvalue(split_crit[,2], vec_levels, vec_type, 2, modified) 
            adj_split_crit = pmax(adj_split_crit1, adj_split_crit2)
        }
        ind_best = order(adj_split_crit)
        ind_best_M = ind_best[1:(min(M,nb_splits))]
        return(list(ind_best_M, adj_split_crit[ind_best_M]))
    }
    else{
        return(list(c(),c()))
    }
}




