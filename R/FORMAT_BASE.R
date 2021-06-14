#### Function to transform ordinal variable into ordinal variable with numbers 
transform_ordinal = function(base,id_cov,order_cov){
    nb_cov_modif = length(id_cov)
    nb_levels = c()
    base_res = data.frame(matrix(NA,nrow=nrow(base),ncol=ncol(base)))
    base_res[,1] = as.numeric(as.character(base[,1]))
    base_res[,2] = as.character(base[,2])
    i_nm = setdiff(3:ncol(base),id_cov)
    for(cov in i_nm){
        base_res[,cov] = as.numeric(as.character(base[,cov]))
    }
    for(icov in 1:nb_cov_modif){
        nb_levels = length(which(!is.na(order_cov[icov,])))
        for(j in 1:nb_levels){
            base_res[which(base[,id_cov[icov]]==order_cov[icov,j]), id_cov[icov]] = j-1
        }
    }
    names(base_res) = names(base)
    return(base_res)
}

