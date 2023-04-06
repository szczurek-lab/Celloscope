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


tworz_ramke <- function(h, x, y, type_nr){
  df <- data.frame(x=x,y=y, frequency= h[, type_nr])
  df <- df[order(df$frequency), ]
  head(df)
  
  #kolory <-  c("#000080",  "#ffff99",  "gold", "#cc00cc", "#cc00cc","#cc00cc", "#cc00cc")
  kolory <- c("#000080", "#FFFF99",  "gold","#CC00CC",  "#CC00CC", "#680068",  "#680068")
  
  marginesy <- c(0, 0.01, -0.5, -0.5)
  #if ( type_nr %in% c(4,8,12,16)) marginesy <- c(0, 0.03, -0.5, -0.5)
  
  ggplot(df, aes(x =x, y = -y)) + geom_point(aes(color =frequency), size =1.8) +
    scale_color_gradientn(colors =kolory, breaks=c(0.2,0.4,  0.6, 0.8,  1),  limits=c(0,1), labels = c("0.2","0.4", "0.6", "0.8",  "1")) +
    theme_bw() + labs(x = "", y="") +   ggtitle(typy[type_nr])+
    theme(plot.title = element_text(size=8, hjust = 0, vjust =-2 ), legend.position = "right", legend.key.width=unit(0.4,"cm"),
          legend.background = element_rect(fill='transparent'),
          plot.margin = unit(marginesy, "cm"), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.text=element_text(size=8), legend.title=element_text(size=8), legend.box="horizontal",  legend.key.size = unit(0.4, 'cm'))
}

