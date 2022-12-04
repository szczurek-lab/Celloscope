
library("devtools")
install_github("Danko-Lab/BayesPrism/BayesPrism")
suppressWarnings(library(BayesPrism))

setwd("/content")
dataset_zip = 'more_than_50_dummy'
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
    c_gs <- t(read.csv("C_gs.csv", row.names = 1, header= TRUE))
    true_lambda <- t(read.csv("TRUE_lambda.csv", row.names = 1, header= TRUE))
    cell.type.labels = rownames(true_lambda)
    cell.state.labels = rownames(true_lambda)

    myPrism <- new.prism(
      reference= true_lambda, #sc.dat.filtered.pc, 
      mixture = c_gs, #bk.dat,
      input.type = "GEP", 
      cell.type.labels = cell.type.labels, 
      cell.state.labels = cell.state.labels,
      key = NULL,
      outlier.fraction=1.)

    bp.res <- run.prism(prism = myPrism, n.cores=40)

    theta <- get.fraction (bp=bp.res,
                which.theta="final",
                state.or.type="type")

    write.table(theta,file=paste(sub_dir, "csv", sep="."))
  
}

setwd("/content")
zip(paste(dataset_zip, 'zip', sep='.'), dataset_zip)