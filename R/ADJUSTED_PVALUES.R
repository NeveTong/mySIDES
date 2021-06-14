######################
# PVALUES ADJUSTMENT #
######################

library(multicool)
library(memoise)


################
#### Functions to calculate the average correlation r*

#### abs_corr: For splliting criterion 1 (maximum differential effect) only
abs_corr =function(r){
    return( ( 2/pi * ( (1 - r^2)^0.5 + r*asin(r)-1 ) ) / (1 -2/pi) )
}


#### For nominal covariates and splitting criterion 1
ro_nominal_c1 = function(k, kj, ki, kij){
    v = 0.5 * ( kij/(ki*kj)^0.5 + (k-kj-ki+kij)/((k-ki)*(k-kj))^0.5 - (kj-kij)/(kj*(k-ki))^0.5 - (ki-kij)/((k-kj)*ki)^0.5 )
    return(abs_corr(v))
}

average_r_nominal_c1 = function(nb_levels){
    r_star = 0
    if(nb_levels >= 3){
        R = 0
        k = nb_levels - 1
        for(l in 1:k){
            if(l-1 >= max(0,2*l-k)){
                for(h in max(0,2*l-k):(l-1)){
                    R = R + ro_nominal_c1(k+1, l, l, h) * 0.5 * multicool::multinom(c(h,l-h,l-h,k-2*l+h), counts=TRUE) 
                }
            }
            if(k >= l+1){
                for(f in (l+1):k){
                    for(h in max(f+l-k, 0):l){
                        R = R + ro_nominal_c1(k+1, l, f, h) * multicool::multinom(c(h,l-h,f-h,k-l-f+h), counts=TRUE)
                    }
                }
            }
        }
        m = 2^k - 1
        r_star = R / (m * (m - 1) / 2)
    }
    return(r_star)
}


#### For nominal covariates and splitting criterion 2
ro_nominal_c2 = function(k, kj, ki, kij){
    r = 0.5 * ( kij/(ki*kj)^0.5 + (k-kj-ki+kij)/((k-ki)*(k-kj))^0.5 + (kj-kij)/(kj*(k-ki))^0.5 + (ki-kij)/((k-kj)*ki)^0.5 )
    return(r)
}

average_r_nominal_c2 = function(nb_levels){
    r_star = 0
    if(nb_levels >= 3){
        m = 2^nb_levels - 2 
        R = 0
        k = nb_levels - 1
        for(l in 1:k){
            if(l-1 >= max(0,2*l-k)){
                for(h in max(0,2*l-k):(l-1)){   
                    R = R + ro_nominal_c2(k+1, l, l, h) * 0.5 * multicool::multinom(c(h,l-h,l-h,k-2*l+h), counts=TRUE)
                }
            }
            if(k >= l+1){
                for(f in (l+1):k){
                    for(h in max(f+l-k, 0):l){
                        R = R + ro_nominal_c2(k+1, l, f, h) * multicool::multinom(c(h,l-h,f-h,k-l-f+h), counts=TRUE)
                    }
                }
            }
        }
        r_star = R / (m * (m - 1) / 2)
    }
    return(r_star)
}


#### For ordinal covariates and splitting criterion 1
ro_ordinal_c1 = function(k, i, j){
  v = 0.5 * ( i/(i*j)^0.5 + (k-j)/((k-i)*(k-j))^0.5 - (j-i)/((k-i)*j)^0.5 )
  return(abs_corr(v))
}

average_r_ordinal_c1 = function(nb_levels){
    r_star = 0
    if(nb_levels >= 3){
        R = 0
        k=nb_levels
        for(i in 1:(k-1)){
            if(k-1 >= i+1){
                for(j in (i+1):(k-1)){
                    R = R + ro_ordinal_c1(k, i, j)
                }
            }
        }
        r_star = 2 * R / ((k-1)*(k-2))
    }
    return(r_star)
}



#### For ordinal covariates and splitting criterion 2
ro_ordinal_c2 = function(k, i, j){
    r = 0.5 * ( i/(i*j)^0.5 + (k-j)/((k-i)*(k-j))^0.5 + (j-i)/((k-i)*j)^0.5 )
    return(r)
}

average_r_ordinal_c2 = function(nb_levels){
    R = 0
    if(nb_levels >= 3){
        k = nb_levels
        for(i in 1:(k-1)){
            if(k-1 >= i+1){
                for(j in (i+1):(k-1)){
                    R = R + ro_ordinal_c2(k, i, j)
                }
            }
        }
        r_star = R / ((k-1)*(2*k-3))
    }
    return(r_star)
}


#### Function to calculate the effective number of split (G^(1-r*))
effective_nb_splits = function(nb_levels, type_var, num_crit){
    G = NA
    r_star = 0
    if(type_var=="nominal"){
        G = 2^(nb_levels-1)-1
        if(num_crit==1){  
            r_star = average_r_nominal_c1(nb_levels)
        }
        else if(num_crit==2){
            r_star = average_r_nominal_c2(nb_levels)
        }
    }
    else if(type_var=="ordinal" || type_var=="continuous"){
        G = nb_levels-1
        if(num_crit==1){
            r_star = average_r_ordinal_c1(nb_levels)
        }
        else if(num_crit==2){
            r_star = average_r_ordinal_c2(nb_levels)
        }
    }
    return(G^(1-r_star))
}

effective_nb_splits_memo = memoise(effective_nb_splits)


#### Function to calculate the number of split G
number_splits = function(nb_levels, type_var){
  G = NA
  if(type_var=="nominal"){
    G = 2^(nb_levels-1)-1
  }
  else if(type_var=="ordinal" || type_var=="continuous"){
    G = nb_levels-1
  }
  return(G)
}


#### Function to calculate the adjusted p_value
adjusted_pvalue = function(vec_pvalues, vec_levels, vec_type, num_crit, modified){
    nb_pval = length(vec_pvalues)
    adj_pval = rep(NA, nb_pval)
    for(i in 1:nb_pval){
        if(vec_levels[i] >= 3){
            if(modified==TRUE){
                adj_pval[i] = 1-(1-vec_pvalues[i])^(effective_nb_splits_memo(vec_levels[i],vec_type[i],num_crit))
            }
            else{
                adj_pval[i] = 1-(1-vec_pvalues[i])^(number_splits(vec_levels[i],vec_type[i]))
            }
        }
        else{
            adj_pval[i] = vec_pvalues[i]
        }
    }
    return(adj_pval)
}



