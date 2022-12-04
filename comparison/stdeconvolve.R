require(remotes)
system("apt install libgsl-dev", intern=TRUE)
install.packages("topicmodels")
remotes::install_github('JEFworks-Lab/STdeconvolve')
library(STdeconvolve)

my_working_dir = '/content'
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

    c_gs = as.matrix(read.csv("C_gs.csv", row.names = 1, header= TRUE))
    true_lambda<-read.csv("TRUE_lambda.csv", row.names = 1, header= TRUE)

    cd <- as.matrix(read.csv("C_gs.csv", row.names = 1, header= TRUE)) # df with genes as rownames and spot ids as colnames

    ## choose optimal number of cell-types
    ldas <- fitLDA(t(as.matrix(cd)), Ks = ncol(true_lambda), ncores=parallel::detectCores()) #ks - liczbe typow komorek

    ## get best model results
    optLDA <- optimalModel(models = ldas, opt = "min")
    ## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
    results <- getBetaTheta(optLDA, perc.filt = 0.)
    deconProp <- results$theta
    deconGexp <- results$beta
    write.csv(t(deconGexp), file=paste('marker_gene', "csv", sep="."))
    write.csv(deconProp, file=paste(sub_dir, "csv", sep="."))
}

setwd("/content")
zip('dense.zip', 'dense')