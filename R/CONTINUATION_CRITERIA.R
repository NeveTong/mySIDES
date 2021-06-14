#########################
# CONTINUATION CRITERIA #
#########################

#### Continuation criterion
continuation_met = function(p_val_parents, p_val_child, gamma){
    parent = FALSE
    if(p_val_child <= gamma * p_val_parents){
       parent = TRUE
    }
    return(parent)
}

