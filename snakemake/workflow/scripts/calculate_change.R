library(ape)
library(phytools)
library(RRphylo)

# ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

species_tree_path = args[1]
busco_prot_tree_path = args[2]
feature_matrix_path = args[3]
busco_id = args[4]
out_file = args[5]

# FUNCTIONS
load_busco_matrix = function(busco_tsv){
    busco_df = read.csv(busco_tsv, sep='\t', row.names=1)
    colnames(busco_df) = unlist(lapply(colnames(busco_df), simplify))
    return(busco_df)
}
simplify = function(s){
    v = unlist(strsplit(s, '_'))
    return(paste(v[1], v[2], sep='_'))
}

keep.tip = function(phy, tip){
    all_tips = phy$tip.label
    drop = setdiff(all_tips, tip)
    return(drop.tip(phy,drop))
}

change_on_branch = function(row){
    parent = as.character(row[1])
    parent_val = acr[parent]
    child = row[2]
    # if internal branch
    if (as.character(child) %in% names(acr)){
        child = as.character(child)
        child_val = acr[child]
    }
    # if external branch
    else{
        child_val = species_vals[child]
    }
    abs_change = abs(parent_val - child_val)/parent_val
    return(abs_change)
}

# LOAD DATA
busco_prot_tree = read.tree(busco_prot_tree_path)
if (is.null(busco_prot_tree)){
    res = c(NA,NA,NA)
}else{
    species_list = sapply(strsplit(busco_prot_tree$tip.label, '__'), function(x) simplify(x[1]))
    
    feature_matrix = read.csv(feature_matrix_path, sep='\t', row.names=1)
    colnames(feature_matrix) = unlist(lapply(colnames(feature_matrix), simplify))
    species_vals = as.numeric(feature_matrix[busco_id,species_list])
    names(species_vals) = species_list
    
    full_species_tree = read.tree(species_tree_path)
    full_species_tree$tip.label = sapply(full_species_tree$tip.label, simplify)
    busco_tree = keep.tip(full_species_tree, species_list)
    
    # RRPHYLO
    rrphylo_res = tryCatch({
         RRphylo(tree = busco_tree, y = species_vals)
    }, error=function(err){
        if (conditionMessage(err) == "task 1 failed - \"L-BFGS-B needs finite values of 'fn'\""){
            print('RRphylo optimization error')
            return(NULL)
        }
    })
    if (is.null(rrphylo_res)){
        res = c(NA,NA,NA)
    }else{
    acr = as.numeric(rrphylo_res$aces)
    names(acr) = rownames(rrphylo_res$aces)
    
    # CALCULATE CHANGE
    species_tree_length = sum(busco_tree$edge.length)/1000	# In Byr
    total_prot_change = sum(busco_prot_tree$edge.length)
    total_feature_change = sum(apply(busco_tree$edge, 1, change_on_branch))
    
    res = c(species_tree_length,total_prot_change,total_feature_change)
    }
}

# WRITE OUT
names(res) = c('Total_tree_length','Total_sequence_change','Total_gene_structure_change')
res_df = t(as.data.frame(res))
rownames(res_df) = c(busco_id)
write.table(res_df, sep='\t', quote=FALSE, file=out_file)
