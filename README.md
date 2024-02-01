# WFCC
Code for the algorithm Weight-Filtration Comparison Curve (WFCC) to compare weighted graphs with equal labelled vertices and edges, but different weightings. The algorithm takes as inputs the graphs in the form of lists of tuples (start vertex, end vertex, weight) corresponding to edges. The algorithm outputs a curve visualising the number of non-common edges over an edge-filtration. 

## Description of files
The WFCC algorithm and auxiliary code are found in `main_Functions.py`. An example of how to use the algorithm is in `Example_WFCC.ipynb`, there, one can introduce their data files in cell [1] '#open datafiles' and proceed with the comparison.

The synthetic hypercubes we considered, with probability and flux weightings, are obtained in `diff-cubes.R`. The file `main.ipynb` contains the code explained, and it reproduces all the figures and computations in the main text of the paper, except for the MDS. The MDS are produced in `mds.R` and in `analysis-general.R`. 

`analysis-general.R` supports analysis of anti-microbial resistance data from the BV-BRC database https://www.bv-brc.org . First, visit the page corresponding to a bacterium of interest (for example, *Mycobacterium tuberculosis* here https://www.bv-brc.org/view/Taxonomy/1773#view_tab=genomes). Download the contents of the "Genomes" and "AMR Phenotypes" tables as CSV files. Label them `BVBRC_genome_[label].csv` and `BVBRC_genome_amr_[label].csv` respectively. Then set the label in `analysis-general.R` and run the code from the directory in which the data is saved.

The folder `Data_files_hypercubes` contains all the data files needed. `Synthetic_data_to_run_HMM_TraPS` contains the synthetic data mimicking clinical data, to be run on HyperHMM and HyperTraPS, and the corresponding outputs are found in `results_synthetic_data` for HyperHMM and `hypertrapsTransMatrices` for HyperTraPS. `Synthetic_data_all_zeros_to_run_HMM_TraPS` is synthetic data mimicking clinical where all samples are the string of all zeros, to be run on HyperHMM and HyperTraPS, and the corresponding outputs are `Hyperhmm_sample_zero` and `Hypertraps_sample_zero`. Lastly, `Outputs_github` is an empty folder to direct there the plots in png form and the data files in csv form that are generated from the code.

HyperTraPS was run using https://github.com/StochasticBiology/hypertraps-simple. HyperHMM was run using https://github.com/StochasticBiology/hypercube-hmm. 

### Citation
If one uses the WFCC algorithm, it can be cited as

	@article {Garc{\'\i}a Pascual2024.01.29.577802,

	author = {Bel{\'e}n Garc{\'\i}a Pascual and Lars M Salbu and Jessica Renz and Konstantinos Giannakis and Iain Johnston},
 
	title = {Comparing structure and dynamics of transition graphs by the symmetric difference metric over an edge-filtration},
 
	elocation-id = {2024.01.29.577802},
 
	year = {2024},
 
	doi = {10.1101/2024.01.29.577802},
 
	publisher = {Cold Spring Harbor Laboratory},
 
	URL = {https://www.biorxiv.org/content/early/2024/01/31/2024.01.29.577802},
 
	eprint = {https://www.biorxiv.org/content/early/2024/01/31/2024.01.29.577802.full.pdf},
 
	journal = {bioRxiv}
	}
