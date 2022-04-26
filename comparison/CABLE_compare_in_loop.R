    library(RCTD)
    library(Matrix)
    
    input_dir = '/home/kd/Desktop/porownanie_cable_magisterka/final/preprocessed' 
    dirs = dir(input_dir) # sciezka do folderu z podfolderami danych
    
    # reference_moje <- dgeToSeurat('/home/kd/Desktop/porownanie/Vignette_SC_wlasne') #SC dowolne, byle sie zgadzamly nazwy - i tak je potem podmieniamy
    

    counts <- read.csv(file.path('/home/kd/Desktop/porownanie_cable_magisterka', "dge.csv")) # load in counts matrix
    rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
    meta_data <- read.csv(file.path('/home/kd/Desktop/porownanie_cable_magisterka', "meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
    cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
    cell_types <- as.factor(cell_types) # convert to factor data type
    nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list
    ### Create the Reference object
    reference <- Reference(counts, cell_types, nUMI)

    
    
    for (dir in dirs){
        print(dir)
        path = paste(input_dir, dir, sep="/")
        #puck_moje <- read.SpatialRNA(path) #trzeba poprzesuwac o dane o kulumne przez tego buga 
        #w barcode colmes chyba zjada 1 kolumne ze spotem
        
        counts <- read.csv(file.path(path,"MappedDGEForR.csv")) # load in counts matrix
        coords <- read.csv(file.path(path,"BeadLocationsForR.csv"))
        rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
        rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
        nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
        ### Create SpatialRNA object
        puck_moje <- SpatialRNA(coords, counts, nUMI)

            
        
        #gene_cutoff i fc_cutoff  i reszta na 0 zeby nie odcinalo genow
        myRCTD_moje <- create.RCTD(puck_moje, reference, max_cores = 7, gene_cutoff = 0, fc_cutoff = 0, gene_cutoff_reg = 0, fc_cutoff_reg = 0)
        TRUE_lambda_path = paste(path, 'TRUE_lambda.csv', sep="/")
        TRUE_lambda = read.csv(TRUE_lambda_path, row.names=1)
        
        #nadpisujemy srednie plikiem z prawdiwymi srednimi (odpowiednik fitbullk)   
        myRCTD_moje@cell_type_info$info[[1]] <- TRUE_lambda
                      
        bulk = fitBulk(myRCTD_moje)
        
        simga = choose_sigma_c(bulk)
              
        moje_wynik <- fitPixels(simga, doublet_mode = 'full')
        
        #normalizujemy wynik  
        norm_weights = sweep(moje_wynik@results$weights, 1, rowSums(moje_wynik@results$weights), '/') 
        wyniki_path = paste(input_dir, dir, sep = "/")
        wyniki_path = paste(wyniki_path, '.csv', sep='.')
        write.csv(as.data.frame(as.matrix(norm_weights)), wyniki_path) 
        

    }
    