prepare_matrix_B <- function(ST_data, address_candidate, cell_types, lead_genes, tresh_tau, tresh_rho ){
  nTypes <- length(cell_types)
  list_marker_genes_all_types <- vector("list", length = nTypes)
  for(i in 1:nTypes){
    print(i)
    candidate_markers <- unique(unlist(read.table(paste0(address_candidate, cell_types[i],".txt") )))
    list_marker_genes_all_types[[i]] <- find_markers_in_agreement_one_type(ST_data, candidate_markers,  tresh_tau[i], tresh_rho[i], lead_genes[i])
  }
  
  names(list_marker_genes_all_types) <- cell_types
  all_markers <- unique(unlist(list_marker_genes_all_types))
  
  list2env(list_marker_genes_all_types, envir=globalenv())
  B <- as.data.frame.matrix(table(stack(list_marker_genes_all_types)))
  B<- B[match(all_markers,rownames(B) ),]
  
  B
}

find_markers_in_agreement_one_type <- function(ST_data, candidate_markers,  tresh_kendall, tresh_pearson, lead_gene){
  # sprawdz, czy lead gene in candidate_markers
  temp <- c()
  
  for(i in 1:length(candidate_markers)){
    tryCatch({
      helper <- ST_data[candidate_markers[i],]
      if(sum(helper>4)>5) temp <- c(temp, candidate_markers[i])
      
    }, error= function(e) print(paste0("Gene ", candidate_markers[i], " not found in data"  ))  )
    
  }
  markers_one_type <- temp
  
  tmp1 <- ST_data[lead_gene,]
  
  df <- data.frame(matrix(0, ncol=5, nrow=length( markers_one_type )  ))
  rownames(df) <- markers_one_type
  
  for(i in 1:length(markers_one_type )){
    tmp2 <- ST_data[markers_one_type [i],]
    df[i,1] <- cor(tmp1,tmp2, method = "kendall")
    df[i,2] <- cor(tmp1,tmp2, method = "pearson")
    df[i,3] <- cor(tmp1,tmp2, method = "spearman")
  }
  df[,4] <- df[,1]/max(df[-which(rownames(df)==lead_gene),1])
  df[,5] <- df[,2]/max(df[-which(rownames(df)==lead_gene),2])
  
  colnames(df) <- c("kendall","pearson", "spearman", "normalised K", "normalised P")
  which_we_choose <- (df$"normalised K" >= tresh_kendall) & (df$"normalised P" >= tresh_pearson)
  df_yes <- df[ which_we_choose, ]
  
  rownames(df_yes)
}
