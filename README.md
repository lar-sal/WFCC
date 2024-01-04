# WFCC
Code for the algorithm Weight-Filtration Comparison Curve (WFCC) to compare weighted graphs with equal labelled vertices and edges, but different weightings. The algorithm takes as inputs the graphs in the form of lists of tuples (start vertex, end vertex, weight) corresponding to edges. The algorithm outputs a curve visualising the number of non-common edges over an edge-filtration. 

## Description of files
The WFCC algorithm and auxiliary code are found in `main_Functions.py`. An example of how to use the algorithm is in `Example_WFCC.ipynb`, there, one can introduce their data files in cell [1] '#open datafiles' and proceed with the comparison.

The synthetic hypercubes we considered, with probability and flux weightings, are obtained in `diff-cubes.R`. The file `main.ipynb` contains the code explained, and it reproduces all the figures and computations in the main text of the paper, except for the MDS. The MDS are produced in `mds.R`. 

The folder `Data_files_hypercubes` contains all the data files needed. `Klebsiella_clinical_data` is the original real data to run HyperHMM. `Klebsiella` is the folder with the HyperHMM outputs of the continents. `Klebsiella_subsampled` contains the HyperHMM outputs from subsampled clinical data. `Synthetic_data_to_run_HMM_TraPS` contains the synthetic data mimicking clinical data, to be run on HyperHMM and HyperTraPS, and the corresponding outputs are found in `results_synthetic_data` for HyperHMM and `hypertrapsTransMatrices` for HyperTraPS. `Synthetic_data_all_zeros_to_run_HMM_TraPS` is synthetic data mimicking clinical where all samples are the string of all zeros, to be run on HyperHMM and HyperTraPS, and the corresponding outputs are `Hyperhmm_sample_zero` and `Hypertraps_sample_zero`. Lastly, `Outputs_github` is an empty folder to direct there the plots in png form and the data files in csv form that are generated from the code.

HyperTraPS was run using https://github.com/StochasticBiology/hypertraps-simple. HyperHMM was run using https://github.com/StochasticBiology/hypercube-hmm. 

### Citation
If one uses the WFCC algorithm, it can be cited as

### References 
The clinical data in `Klebsiella_clinical_data` was obtained from the following publications:

1. Savin, M., Bierbaum, G., Hammerl, J. A., Heinemann, C., Parcina, M., Sib, E., Voigt, A., and Kreyenschmidt, J. (2020). ESKAPE bacteria and extended-spectrum-betalactamase-producing escherichia coli isolated from wastewater and process water from german poultry slaughterhouses. *Applied and Environmental Microbiology*, 86(8):e02748–19.
2. Spadar, A., Phelan, J., Elias, R., Modesto, A., Caneiras, C., Marques, C., Lito, L., Pinto, M., Cavaco-Silva, P., Ferreira, H., Pomba, C., Da Silva, G. J., Saavedra, M. J., Melo-Cristino, J., Duarte, A., Campino, S., Perdigão, J., and Clark, T. G. (2022). Genomic epidemiological analysis of klebsiella pneumoniae from portuguese hospitals reveals insights into circulating antimicrobial resistance. *Scientific Reports*, 12(1):13791.
3. Di Tella, D., Tamburro, M., Guerrizio, G., Fanelli, I., Sammarco, M. L., and Ripabelli, G. (2019). Molecular epidemiological insights into colistin-resistant and carbapenemases-producing clinical klebsiella pneumoniae isolates. *Infection and drug resistance*, 12.
4. Prendergast, D. M., O’Doherty, Á., Burgess, C. M., Howe, N., McMahon, F., Murphy, D., Leonard, F., Morris, D., Harrington, C., Carty, A., Moriarty, J., and Gutierrez, M.
(2022). Critically important antimicrobial resistant enterobacteriaceae in irish farm effluent and their removal in integrated constructed wetlands. *Science of The Total Environment*, 806:151269.
5. Walsh, F., Cooke, N. M., Smith, S. G., Moran, G. P., Cooke, F. J., Ivens, A., Wain, J., and Rogers, T. R. (2010). Comparison of two DNA microarrays for detection of plasmidmediated antimicrobial resistance and virulence factor genes in clinical isolates of enterobacteriaceae and non-enterobacteriaceae. *International Journal of Antimicrobial Agents*, 35(6):593–598.
6. Taitt, C. R., Leski, T. A., Erwin, D. P., Odundo, E. A., Kipkemoi, N. C., Ndonye, J. N., Kirera, R. K., Ombogo, A. N., Walson, J. L., Pavlinac, P. B., Hulseberg, C., and Vora, G. J. (2017). Antimicrobial resistance of klebsiella pneumoniae stool isolates circulating in kenya. *PLOS ONE*, 12(6):e0178880.
7. Henson, S. P., Boinett, C. J., Ellington, M. J., Kagia, N., Mwarumba, S., Nyongesa, S., Mturi, N., Kariuki, S., Scott, J. A. G., Thomson, N. R., and Morpeth, S. C. (2017). Molecular epidemiology of klebsiella pneumoniae invasive infections over a decade at kilifi county hospital in kenya. *International Journal of Medical Microbiology*, 307(7):422–429.
8. Norsigian, C. J., Attia, H., Szubin, R., Yassin, A. S., Palsson, B. Ø., Aziz, R. K., and Monk, J. M. (2019). Comparative genome-scale metabolic modeling of metallo-betalactamase– producing multidrug-resistant klebsiella pneumoniae clinical isolates. *Frontiers in Cellular and Infection Microbiology*, 9.
9. Shawa, M., Furuta, Y., Mulenga, G., Mubanga, M., Mulenga, E., Zorigt, T., Kaile, C., Simbotwe, M., Paudel, A., Hang’ombe, B., and Higashi, H. (2021). Novel chromosomal
insertions of ISEcp1-blaCTX-m-15 and diverse antimicrobial resistance genes in zambian clinical isolates of enterobacter cloacae and escherichia coli. *Antimicrobial Resistance & Infection Control*, 10(1):79.
10. Fadare, F. T. and Okoh, A. I. (2021). Distribution and molecular characterization of ESBL, pAmpC beta-lactamases, and non-beta-lactam encoding genes in enterobacteriaceae
isolated from hospital wastewater in eastern cape province, south africa. *PLOS ONE*, 16(7):e0254753.
11. Selmi, R., Tayh, G., Srairi, S., Mamlouk, A., Ben Chehida, F., Lahmar, S., Bouslama, M., Daaloul-Jedidi, M., and Messadi, L. (2022). Prevalence, risk factors and emergence of extended-spectrum beta-lactamase producing-, carbapenem- and colistin-resistant enterobacterales isolated from wild boar (sus scrofa) in tunisia. *Microbial Pathogenesis*, 163:105385.
12. Carlos, C. C., Masim, M. A. L., Lagrada, M. L., Gayeta, J. M., Macaranas, P. K. V., Sia, S. B., Facun, M. A. M., Palarca, J. F. C., Olorosa, A. M., Cueno, G. A. C., Abrudan, M., Abudahab, K., Argimón, S., Kekre, M., Underwood, A., Stelling, J., Aanensen, D. M., and on Genomic Surveillance of Antimicrobial Resistance, N. G. H. R. U. (2021). Genome Sequencing Identifies Previously Unrecognized Klebsiella pneumoniae Outbreaks in Neonatal Intensive Care Units in the Philippines. *Clinical Infectious Diseases*, 73(Supplement 4):S316–S324.
13. Nguyen, T., Nguyen, P., Le, N., Nguyen, L., Duong, T. B., Ho, N., Nguyen, Q., Pham, T. D., Tran, A. T., The, H. C., Nguyen, H. H., Nguyen, C., Thwaites, G. E., Rabaa,
M. A., and Pham, D. T. (2021). Emerging carbapenem-resistant klebsiella pneumoniae sequence type 16 causing multiple outbreaks in a tertiary hospital in southern vietnam. *Microbial genomics*, 7,3.
14. Lomonaco, S., Crawford, M. A., Lascols, C., Timme, R. E., Anderson, K., Hodge, D. R., Fisher, D. J., Pillai, S. P., Morse, S. A., Khan, E., Hughes, M. A., Allard, M. W.,
and Sharma, S. K. (2018). Resistome of carbapenem- and colistin-resistant klebsiella pneumoniae clinical isolates. *PLOS ONE*, 13(6):e0198526.
15. Boonyasiri, A., Jauneikaite, E., Brinkac, L. M., Greco, C., Lerdlamyong, K., Tangkoskul, T., Nguyen, K., Thamlikitkul, V., and Fouts, D. E. (2021). Genomic and clinical
characterisation of multidrug-resistant carbapenemase-producing st231 and st16 klebsiella pneumoniae isolates colonising patients at siriraj hospital, bangkok, thailand from 2015 to 2017. *BMC Infectious Diseases*, 21(1):142.
16. Ngoi, S. T., Teh, C. S. J., Chong, C. W., Abdul Jabar, K., Tan, S. C., Yu, L. H., Leong, K. C., Tee, L. H., and AbuBakar, S. (2021). In vitro efficacy of flomoxef against extendedspectrum beta-lactamase-producing escherichia coli and klebsiella pneumoniae associated with urinary tract infections in malaysia. *Antibiotics*, 10(2):181.
17. Sofiana, E., Effendi, M., Plumeriastuti, H., and Pratama, J. (2021). CASES OF MULTIDRUG RESISTANCE (MDR) IN KLEBSIELLA PNEUMONIAE ISOLATED FROM HEALTHY PIGS. *Biochemical and Cellular Archives*, 21:1979–1985.
18. Yossapol, M., Suzuki, K., Odoi, J. O., Sugiyama, M., Usui, M., and Asai, T. (2020). Persistence of extended-spectrum beta-lactamase plasmids among enterobacteriaceae in
commercial broiler farms. *Microbiology and Immunology*, 64(10):712–718.
19. Peng, Y., Liang, S., Poonsuk, K., On, H., Li, S. W., Maurin, M. M. P., Chan, C. H., Chan, C. L., Sin, Z. Y., and Tun, H. M. (2021). Role of gut microbiota in travel-related acquisition of extended spectrum beta-lactamase-producing Enterobacteriaceae. *Journal of Travel Medicine*, 28(3).
20. Gorrie, C. L., Mirceta, M., Wick, R. R., Judd, L. M., Wyres, K. L., Thomson, N. R., Strugnell, R. A., Pratt, N. F., Garlick, J. S., Watson, K. M., Hunter, P. C., McGloughlin, S. A., Spelman, D. W., Jenney, A. W. J., and Holt, K. E. (2018). Antimicrobial-Resistant Klebsiella pneumoniae Carriage and Infection in Specialized Geriatric Care Wards Linked to Acquisition in the Referring Hospital. *Clinical Infectious Diseases*, 67(2):161–170.
21. Nakamura-Silva, R., Dias, L. L., Sousa, R. C., Fujimoto, R. Y., and Pitondo-Silva, A. (2022). Multidrug-resistant and potentially pathogenic enterobacteriaceae found
in a tertiary hospital sewage in southeastern brazil. *Environmental Monitoring and Assessment*, 194(10):782.

