# Celloscope

This directory contains the source code necessary to reproduce analyses presented in the manuscript:  

> **Celloscope: a probabilistic model for marker-gene-driven cell type deconvolution in spatial transcriptomics data**  
> Agnieszka Geras, Shadi Darvish Shafighi, Kacper Domżał, Igor Filipiuk,  Łukasz Rączkowski, Hosein Toosi, Leszek Kaczmarek,  Łukasz Koperski, Jens Lagergren, Dominika Nowis and Ewa Szczurek

The preprint version of the manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2022.05.24.493193v1).

Cell counting procedure is [here](https://github.com/szczurek-lab/qupath-spot-utils).

**Questions about the implementation and general questions:**

Agnieszka Geras, A.Geras[at]mini.pw.edu.pl

**Manuscript's corresponding author:**

Ewa Szczurek, szczurek[at]mimuw.edu.pl

# Basic usage
Celloscope takes as input the following files provided in the `input_data` directory:

* `param.txt` - file containing run setting (hyperparameters, number of iterations, etc.),
* `matB.csv` - binary matrix with prior knowledge about marker genes,
* `C_gs.csv` - gene expression data,
* `n_cells.csv` - estimates for the number of cells in each ST spot.

Celloscope can be run from the `code` directory as follows:

```
bash  Celloscope-run.sh 'input_data' 'results' number_of_chains
```

The arguments of the bash script:
* `input_data` - directory containing input data,
* `results` - directory dedicated to results,
* `number_of_chains` - number of independent chains to run.

Exemplary input files can be found in this repository in the `example` directory. 

Please see the manual for the details on preparing data and guidance on visualising results.

# Repository content

* `data` - data used for the analysis presented in the paper (ST data from mouse brain and human prostate and estimates for the number of cells for each spot),
* `code` - model's implementation,
* `example` - exemplary input files ready to run on Celloscope,
* `manual` explaining the details of how to prepare data, run the model, and visualize results,
* `comparison` - code used to compare to preceding approaches.
