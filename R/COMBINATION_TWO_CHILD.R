###########################################
# ENUMERATE ALL COMBINATIONS OF TWO CHILD #
###########################################

#### Functions to cut continuous variable in ordinal
ord.cut = function(x,ord.bin){
	x=x[!is.na(x)]
	temp=sort(unique(x))
	temp = temp[-length(temp)]
	if(length(temp) < ord.bin)
		return(list(cut=temp,g=length(temp)))

	ordseq=round(seq(0,length(x),length(x)/ord.bin))
	ordseq=ordseq[c(-1,-length(ordseq))]
	temp=unique(sort(x)[ordseq])
	return(list(cut=temp,g=length(temp)))
}


#### Function that enumerate all possible combinations of two child for a given covariate with its type
comb_child = function(covariate, type, ord.bin=NA) {
    if(type == "nominal"){
        levels = sort(unique(covariate))
        nb_levels = length(levels)
        if(nb_levels>1){
            groups = vector("list", 2^(nb_levels-1)-1)
            for(i in 1:(2^(nb_levels-1)-1)) {
                bits = intToBits(i)[1:nb_levels];
                group1 = levels[which(bits==0)]
                group2 = levels[which(bits==1)]
                groups[[i]] = list(group1, group2)
            }
            return(groups)
        }
        else{
            return(NA)
        }
    } 
    if(type == "ordinal"){
        levels = sort(unique(covariate))
        nb_levels = length(levels)
        if(nb_levels>1){
            groups = vector("list", nb_levels-1)
            for(i in 1:(nb_levels-1)) {
                group1 = paste(levels[i],"-",sep="")
                group2 = paste(levels[i],"+",sep="")
                groups[[i]] = list(group1, group2)
            }
            return(groups)
        }
        else{
            return(NA)
        }
    }
    if(type=="continuous"){
        temp = ord.cut(covariate,ord.bin)
        levels = temp$cut
        nb_levels = temp$g+1
        if(nb_levels>1){
            groups = vector("list", nb_levels-1)
            for(i in 1:(nb_levels-1)) {
                group1 = paste(levels[i],"-",sep="")
                group2 = paste(levels[i],"+",sep="")
                groups[[i]] = list(group1, group2)
            }
            return(groups)
        }
        else{
            return(NA)
        }
    }
}



#### Function that excludes splits for which the size of the child subgroup with the largest treatment effect is < S 
exclude_splits = function(comb_child, pvalues, sizes, S){
    nb_splits = length(comb_child)
    update_comb_child = comb_child
    ind_remove = c()
    for(i in 1:nb_splits){
        if( (pvalues[[i-length(ind_remove)]][1] < pvalues[[i-length(ind_remove)]][2] && sizes[[i-length(ind_remove)]][1] < S) ||
            (pvalues[[i-length(ind_remove)]][2] <= pvalues[[i-length(ind_remove)]][1] && sizes[[i-length(ind_remove)]][2] < S) ){
            update_comb_child[[i-length(ind_remove)]] = NULL
            pvalues[[i-length(ind_remove)]] = NULL
            sizes[[i-length(ind_remove)]] = NULL
            ind_remove = c(ind_remove, i)
        }
    }    
    return(list(update_comb_child, pvalues))
}

