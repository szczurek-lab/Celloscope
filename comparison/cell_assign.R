library(tensorflow)
library(SingleCellExperiment)
library(cellassign)
library(scran)





input_dir = '/home/kd/Desktop/cellassign' 
path = input_dir
#dirs = dir(input_dir) 
learn_rate_vec =  c(1:4)


#for (dir in dirs){
  #print(dir)
  #path = paste(input_dir, dir, sep="/")

  C_gs = read.csv(paste(path, "C_gs.csv", sep='/'), row.names = 1)
  gen_names = rownames(C_gs)
  spot_names = colnames(C_gs)
  C_gs = data.matrix(C_gs)
  C_gs = SingleCellExperiment(assays=list(counts=C_gs))


  B = read.csv(paste(path, "matB.csv", sep='/'), row.names = 1)
  rownames(B) = gen_names
  
  #true_h = read.csv(paste(path, "TRUE_h.csv", sep='/'), row.names = 1)
  #colnames(B) = colnames(true_h)

  B = as.matrix(B)
    

  calculated_sum_factors = calculateSumFactors(C_gs)
  names(calculated_sum_factors) = c(colnames(C_gs))
  
  mtx= matrix(nrow=length(spot_names))


  for (i in learn_rate_vec){
    learn_rate = 10^-i
    print(learn_rate)
    fit <- cellassign(exprs_obj = C_gs[rownames(B), ], # selecting only marker genes
                      marker_gene_info = B, # macierz B
                      s = calculated_sum_factors, #  Numeric vector of cell size factors
                      learning_rate = learn_rate, 
                      shrinkage = TRUE, # logical - should the delta parameters have hierarchical shrinkage?
                      verbose = FALSE)
    #write.csv(celltypes(fit), file=paste(input_dir,'test.csv', sep=''))

    mtx <- cbind(as.vector(fit$cell_type), mtx)

  #}
  df = as.data.frame(mtx)
  #df = df[,learn_rate_vec]
  #colnames(df) = learn_rate_vec
  #df$true_max_h = colnames(true_h)[apply(true_h,1,which.max)]
  #results_path = paste(input_dir, paste(dir, '.csv', sep=''), sep = '')
  write.csv(df, paste(paste(paste(input_dir, "mouse_brain_jan", sep="/"), i, sep="_"), ".csv", sep='')) 
}
