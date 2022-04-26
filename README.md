# Celloscope
Celloscope: a probabilistic model for marker-gene-driven cell type deconvolution in spatial transcriptomics data

This is the directory containing source code necessary to reproduce analyses presented in the manuscript:  

> **Celloscope: a probabilistic model for marker-gene-driven cell type deconvolution in spatial transcriptomics data**  
> Agnieszka Geras, Shadi Darvish Shafighi, Kacper Domżał, Igor Filipiuk,  Łukasz Rączkowski, Hosein Toosi, Leszek Kaczmarek,  Łukasz Koperski, Jens Lagergren, Dominika Nowis and Ewa Szczurek

# Basic usage
Celloscope takes as an input the following files provided in `input_data` directory:

* `param.txt` - file containing run setting (values for hiperparameters, number of iterations etc.),
* `matB.csv` - binary matrix with prior knowledge on marker genes,
* `C_gs.csv` - gene expression data,
* `n_cells.csv` - estimates for the number of cells in each spot.

Celloscope can be run from the `code` directory as follows:

```
bash  Celloscope-run.sh 'input_data' 'results' number_of_chains
```

The arguments of the bash script:
* `input_data` - directory containing input data,
* `results` - directory dedicated for results.
* `number_of_chains` - number of independent chains to run

Celloscope gives as an output the estimates for the cell type proportions across as a spots vs. types matrix  in the file `h_est.csv` that can be found in `results` directory.

Examplificatory input files can be found in this resposotory, in the `example` directory. 

Please see the manual for the details on how to prepare data and guidance to visualize results.

# Repository content

* `markers` - files containing marker genes used in the analysis described in the paper (Additional File 1 and Additional File 2) and marker genes for the example from the manual,
* `data` - data used for the analysis presented in the papaer (ST data from mouse brain and human prostate and estimates for the number of cells for each spot),
* `code` - model's implementation,
* `example` - examplificatory input files ready to run on Celloscope,
* `manual` explaining the details how to prepare data, run the model and visualise results,
* `comparison` - code used to perform comparison to preceding aproaches
