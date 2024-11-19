# SwarmTCRTools

SwarmTCRTools is a collection of R scripts intended to facilitate automated analysis of TCR sequence binding prediction with [SwarmTCR and SwarmTCRClassify](https://github.com/thecodingdoc/SwarmTCR). 

## Workflow Overview

The workflow for using SwarmTCRTools involves several key steps: preparation of experimental data, preparation of reference data, SwarmTCR pretraining, processing SwarmTCR results, SwarmTCRClassify, and processing SwarmTCRClassify results.

## Scripts

### [Experimental Data Preparation](experimental-preparation.Rmd)

The `experimental-preparation.Rmd` script prepares experimental data for use in SwarmTCRClassify. The data should be in the form of a data frame with columns for Clone ID, CDR sequences (CDR1a, CDR2a, CDR3a, CDR1b, CDR2b, CDR3b), and V and J gene segments (AV, AJ, BV, BJ). This script imports the necessary libraries and utility functions, reads the experimental data from a specified path, formats the data to ensure it meets the requirements for SwarmTCRClassify, and saves the formatted data as a space-delimited text file for further processing.

### [Reference Data Preparation](reference-preparation.Rmd)

The `reference-preparation.Rmd` script prepares reference data for use with SwarmTCR. The reference data should be in the form of a data frame with columns for Clone ID, CDR sequences, and Epitope. This script imports the necessary libraries and utility functions, reads the positive and negative reference data from specified paths, selects eligible epitopes based on predefined criteria, generates negative reference data by creating random sequences or using predefined negative sequences, splits the data into training and testing sets, and creates a reference file for SwarmTCRClassify by combining the positive and negative reference data.

### [SwarmTCR Pretraining](swarmTCR-pretraining.Rmd)

To perform SwarmTCR pretraining using the formatted reference data, run the following command in the terminal:

```sh
for f in ./*/; do cd $f; for fi in ./*/; 
do cd $fi; <path to>/swarmTCR -r Reference.txt -i Train.txt -t Test.txt -1 val_TCRdist_out.txt -2 val_SwarmTCR_out.txt -n 20 -s 25 > refWeights.txt;
cd ..; done; cd ..; done;
```

This command iterates through the directories containing the reference data, runs the SwarmTCR pretraining for each set of data, and saves the reference weights to a file. The parameters used in the command specify the reference file, training file, testing file, output files for TCRdist and SwarmTCR, the number of iterations, and the seed for random number generation.

### [Processing SwarmTCR Results](SwarmTCR-results-processing.Rmd)

The `SwarmTCR-results-processing.Rmd` script processes the results of SwarmTCR pretraining. It imports the necessary libraries and utility functions, reads the results files generated by SwarmTCR pretraining, processes the results to extract relevant information such as precision and recall values, combines the results into a master file for further analysis, visualizes the performance of the SwarmTCR model using precision-recall curves, and sets the per-epitope CDR sequence reference weights and score thresholds for SwarmTCRClassify. These weights and thresholds are set in an automated manner based on methods selected by the user; however, it is vital that these weights are verified for logical consistency and biological relevance. These weights may be changed manually if necessary.

### [SwarmTCRClassify](swarmTCRClassify.Rmd)

To classify new experimental data using the pre-trained SwarmTCR model, run the following command in the terminal:

```sh
for f in ./*_CDRs.txt; do <path to>/swarmTCRClassify -w <path to>/reference-data/SwarmTCR_refWeights.txt -i $f > $f.SwarmTCR.results; done;
```

This command iterates through the experimental data files, runs SwarmTCRClassify for each file, and saves the classification results to a file. The parameters used in the command specify the reference weights file and the input file.

### [Processing SwarmTCRClassify Results](SwarmTCRClassify-results-processing.Rmd)

The `SwarmTCRClassify-results-processing.Rmd` script processes the results of SwarmTCRClassify. It imports the necessary libraries and utility functions, reads the results files generated by SwarmTCRClassify, joins the results with the full receptor data to add additional information, applies score thresholds to filter the results based on predefined criteria, writes the processed results to files for further analysis, and visualizes the performance of the SwarmTCRClassify model using precision-recall curves and other metrics.

### Utility Functions

Several utility functions are used throughout the scripts to facilitate data processing and analysis. These functions are sourced from the TCRPaired repository or defined in the [SwarmTCRTools/util](./util/) directory and include functions for data import, processing, filtering, formatting, aligning, and plotting.

## Conclusion

SwarmTCRTools provides a comprehensive workflow for TCR sequence binding prediction using SwarmTCR and SwarmTCRClassify. By following the steps outlined in this README and using the provided scripts, users can prepare their data, perform SwarmTCR pretraining, classify new experimental data, and analyze the results.

For more detailed information on each step and the specific functions used, please refer to the individual R scripts and their documentation.
