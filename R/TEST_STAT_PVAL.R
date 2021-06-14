########################################
# Treatment effect that return pvalues #
########################################

#### Function that substract the subgroup of interest from parents
sub_sets_parents = function(set, parents){
    set_cur = set 
    if(length(parents[[1]]) > 0){
        for(i in 1:length(parents[[1]])){
            ind_pat = c()
            if(parents[[3]][i] == "nominal"){
                for(j in 1:length(parents[[2]][[i]])){
                    ind_pat = c(ind_pat, which(set_cur[,parents[[1]][i]] == parents[[2]][[i]][j]))
                }
            }
            else if(parents[[3]][i] == "ordinal" || parents[[3]][i] == "continuous"){
                val_cut = as.numeric(substr(parents[[2]][[i]],1,nchar(parents[[2]][[i]])-1))
                signe = substr(parents[[2]][[i]],nchar(parents[[2]][[i]]),nchar(parents[[2]][[i]]))
                if(signe == "-"){
                    ind_pat = c(ind_pat, which(set_cur[,parents[[1]][i]] <= val_cut) )
                }
                else{
                    ind_pat = c(ind_pat, which(set_cur[,parents[[1]][i]] > val_cut) )
                }
            }
            set_cur = set_cur[ind_pat,]  
        }
    }
    return(list(set_cur))
}


#### Function that substract the subgroup of interest from parents to the two studied child from the global set
sub_sets_all_childs = function(set, parents, covariate, list_two_child, type_cov){
    set_cur = sub_sets_parents(set, parents)[[1]]
    nb_splits = length(list_two_child)
    set_subgroup_child = vector("list", nb_splits)
    for(split in 1:nb_splits){  
        set_subgroup_child[[split]] = vector("list", 2) 
        set_subgroup_child[[split]][[1]] = sub_sets_parents(set_cur, list(covariate,list(list_two_child[[split]][[1]]),type_cov))[[1]]
        set_subgroup_child[[split]][[2]] = sub_sets_parents(set_cur, list(covariate,list(list_two_child[[split]][[2]]),type_cov))[[1]]
    }
    return(set_subgroup_child)
}




#### Function that return treatment T-stat and p_value from a given data set
analyse = function(set_groupe, type_outcome, level_control, D=0, alpha, upper_best=TRUE){
    if(is.null(dim(set_groupe))){
        set_groupe = t(as.matrix(set_groupe))
    }
    outcome = set_groupe[,1]
    treatment = set_groupe[,2]
    level_trt = unique(treatment)
    ind_pat_trt1= which(treatment!=level_control)
    ind_pat_trt2= which(treatment==level_control)
    n1 = length(ind_pat_trt1)
    n2 = length(ind_pat_trt2)
    outcome_trt1 = outcome[ind_pat_trt1]
    outcome_trt2 = outcome[ind_pat_trt2]
    if(n1 > 1 && n2 > 1){
        if(type_outcome=="continuous"){
            mean1 = mean(outcome_trt1)
            mean2 = mean(outcome_trt2)
            var1 = var(outcome_trt1)
            var2 = var(outcome_trt2)
            var_pool = ( (n1-1)*var1+(n2-1)*var2 ) / (n1+n2-2)
            zeff = (mean1-mean2-D) / sqrt(var_pool*(1/n1+1/n2))
            #temp = summary(lm(outcome~treatment))
            #zeff = temp$coef[2,3]
            #pval_eff = 1-pt(zeff,temp$df[2])
            pval_eff = 1-pt(zeff,n1+n2-2)
            eff_met = (pval_eff<alpha)
        }
        else if(type_outcome=="binary"){
            prop1 = mean(outcome_trt1)
            prop2 = mean(outcome_trt2)
            if(prop1*n1 > 5 && prop2*n2 > 5 && (1-prop1)*n1 > 5 && (1-prop2)*n2 > 5){
                prop_commune = (prop1*n1+prop2*n2)/(n1+n2)
                var_pool = prop_commune*(1-prop_commune)
                zeff = (prop1-prop2-D) / sqrt(var_pool*(1/n1+1/n2)) 
                pval_eff = 1-pnorm(zeff,0,1)
            }
            else{
                mat = matrix(c(round(prop1*n1),round((1-prop1)*n1),round(prop2*n2),round((1-prop2)*n2)),ncol=2)
                test = fisher.test(mat,alternative="greater")
                pval_eff = test$p
                zeff = qnorm(1-pval_eff)      
            }
            eff_met = (pval_eff<alpha)
        }
        else if(type_outcome=="survival"){
            zeff = summary(coxph(Surv(outcome[,1], outcome[,2]) ~ treatment))$coef[4]
            pval_eff = pnorm(zeff,0,1)
            eff_met = (pval_eff<alpha)
        }
        else if(type_outcome=="count"){
            zeff = summary(glm.nb(outcome ~ treatment))$coef[6]
            pval_eff = 1-pnorm(zeff,0,1)
            eff_met = (pval_eff<alpha)
        }
    }
    else{
        zeff = 0
        pval_eff = 1
        eff_met = FALSE
    }
    if(upper_best == FALSE){
        zeff = -zeff
        pval_eff = 1-pval_eff
    }
    if(is.na(pval_eff)){
        zeff = 0
        pval_eff = 1
        eff_met = FALSE
    }
    return(c(zeff,pval_eff,eff_met))
}




