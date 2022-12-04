shell_call <- function(command, ...) {
  result <- system(command, intern = TRUE, ...)
  cat(paste0(result, collapse = "\n"))}
shell_call("sudo apt install libgsl-dev")
shell_call("sudo apt-get install libgmp3-dev")
shell_call("pip install igraph")
shell_call("pip install leidenalg")
install.packages('quadprog')
shell_call('add-apt-repository -y ppa:cran/imagemagick')
shell_call('apt-get update')
shell_call("apt-get install libmagick++-dev")
install.packages("magick")
require(remotes)
remotes::install_github("RubD/Giotto") 
library(Giotto)

shell_call("which python3")
python_path = '/usr/bin/python3'
my_working_dir = '/content'
instrs = createGiottoInstructions(python_path = python_path, show_plot = FALSE)
setwd("/content")
dataset_zip = 'dense'
unzip(paste(dataset_zip, 'zip', sep='.'))

main_dir = paste("/content/", dataset_zip, sep="")
setwd(main_dir)
folders <- as.vector(list.dirs())
folders <- grep('setup', folders, value = TRUE)
folders = (folders[!grepl("checkpoints", folders)])
print(folders)

for (sub_dir in folders) {
    print(sub_dir)
    setwd(paste(main_dir, sub_dir, sep="/"))
    
    # simulating scRNAseq, but giving true lambda instead, so one type has one cell
    sc_matrix<-read.csv("TRUE_lambda.csv", row.names = 1, header= TRUE)
    # data frame with col containing vector with single cell labels
    sc_lable = as.data.frame(colnames(read.csv("TRUE_lambda.csv", row.names = 1, header= TRUE)))
    colnames(sc_lable) = 'single_cell_labels'

    sc_cortex <- createGiottoObject(raw_exprs = sc_matrix, instructions = instrs)
    sc_cortex <- normalizeGiotto(gobject = sc_cortex)
    # using single cell labels as leiden cluster like in the https://github.com/rdong08/spatialDWLS_dataset/blob/main/codes/seqFISH_plus_deconvolution.Rmd
    sc_cortex@cell_metadata$leiden_clus<-as.character(sc_lable$single_cell_labels)

    Sig = as.matrix(read.csv("TRUE_lambda.csv", row.names = 1, header= TRUE))

    grid_exp = read.csv("C_gs.csv", row.names = 1, header= TRUE)
    grid_seqFish <- createGiottoObject(raw_exprs = grid_exp, instructions = instrs)
    grid_seqFish <- normalizeGiotto(gobject = grid_seqFish)
    grid_seqFish <- calculateHVG(gobject = grid_seqFish)
    gene_metadata = fDataDT(grid_seqFish)
    gene_metadata$hvg = 'yes'

    featgenes = gene_metadata[hvg == 'yes']$gene_ID
    grid_seqFish <- runPCA(gobject = grid_seqFish, genes_to_use = featgenes, scale_unit = F)
    grid_seqFish <- createNearestNetwork(gobject = grid_seqFish, dimensions_to_use = 1:10, k = length(colnames(sc_matrix)))
    grid_seqFish <- doLeidenCluster(gobject = grid_seqFish, resolution = 0.4, n_iterations = 1000)
    grid_seqFish<-runDWLSDeconv(gobject = grid_seqFish, sign_matrix = Sig)
    
    res_df = grid_seqFish@spatial_enrichment$DWLS
    rownames(res_df) = res_df$cell_ID
    res_df = res_df[,-c(1)]
    write.csv(res_df, file=paste(sub_dir, "csv", sep="."))
}

setwd("/content")
zip('dense.zip', 'dense')