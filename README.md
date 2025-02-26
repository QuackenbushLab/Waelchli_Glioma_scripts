# Single Cell Analysis Pipeline

1.  Download the PPI file and the motif file from GRAND: https://grand.networkmedicine.org/tissues/.
2.  Download the Waelchli dataset from GEO (GSE256493).
3.  Install the [SCORPION](https://github.com/kuijjerlab/SCORPION) and [NetZoo](https://github.com/netZoo/netZooR) R packages.
4.  Run **add_symbols_to_motif.R** to add gene symbols to the motif data, changing the following file paths:
    -   **motifFile:** The path to the file with the motif data.
    -   **motifOutFile:** The path to the file where you wish to save the motif data with the symbols.
5.  Modify the **run_scorpion.R** file to include the following file paths:
    -   **motifFile:** The path to the file with the motif data.
    -   **ppiFile:** The path to the file with the PPI data.
6.  Run **run_scorpion.sh** to run all SCORPION models. This will call the **run_scorpion.R** script for each data file. Change the following file paths as needed:
    -   **EXPRESSION_DIR:** The directory containing the downloaded GEO data.
    -   **OUTPUT_DIR:** The directory where you wish to save your output.
    -   **LOGGING_DIR:** The directory where you wish to save your logs
7.  Generate a null PANDA distribution using **netZooR::GenerateNullPANDADistribution**, setting **ppiFile** to your PPI file and **motifFile** to your motif file.
8. Run **run_BLOBFISH.R**, changing the following variables as needed:
    -   **sourceDirectory:** The directory containing the SCORPION output.
    -   **resultDirectory:** The directory where you wish to save your results.
    -   **allGenes:** The list of target genes of interest.
    -   **nullFile:** The file containing the null PANDA distribution generated in the previous step.
9. Run **figures_sox4_new_gene_set_HSPG2_new_immunosuppressive_set.R** to generate all plots, changing the following variables:
    -   **networkDir:** The directory where the BLOBFISH results are stored.
    -   **combMatDir:** The directory where the combinatorial matrices used in generating the UpSet plots should be stored.
