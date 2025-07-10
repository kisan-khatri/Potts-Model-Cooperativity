# Using AlphaFold2 to Predict the Conformations of Side Chains in Folded Proteins
DOI: https://doi.org/10.1101/2025.02.10.637534

This repository contains the code and data used in our study on double mutant cycle (cooperativity) analysis in kinase proteins, as presented in the manuscript. It includes the code to reproduce the histogram shown in the manuscript and to identify the most cooperative double mutant pairs reported. The repository is shared in response to reviewer requests to promote transparency and reproducibility. Permission to share these materials was granted by my advisor, [Prof. Dr. Ronald M. Levy](https://scholar.google.com/citations?user=6CQ_uloAAAAJ&hl=en), who is the corresponding author of the manuscript.

We used a Potts model of kinase protein families to compute a histogram of double mutation energies, ΔΔE, and to identify the most cooperative mutation pairs in two kinases: ABL1 and PIM1.

# How to cite this code
If you think this has contributed to the work you are doing, consider citing it in the list of your references. Here is the recommended citation:

Khatri, K., Haldane, A., & Levy, R. M. (2025). kisan-khatri/Potts-Model-Cooperativity: v1.0.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.15856333

Zenodo. https://doi.org/10.5281/zenodo.15856333
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15856333.svg)](https://doi.org/10.5281/zenodo.15856333)
# Instruction

For a detailed flow path, please refer to the README.md file located in the Necessary_Files folder. However, I am providing a summary of the steps below.
# Step 1: Software Used
To perform this analysis, we used the [Mi3-GPU software](https://github.com/ahaldane/Mi3-GPU), which is publicly available on GitHub.  
Follow the installation and usage instructions provided in the 'README.md' file inside the folder 'Necessary_Files'.
# Step 2: Construct Kinase MSA and Fit the Potts Model
The multiple sequence alignment (MSA) used in this study is 'realigned259_aln', and the fitted Potts model coupling matrix is provided as 'J.npy'. These files, along with all other necessary data, scripts, and instructions, are provided in the folder 'Necessary_Files'.  
# Step 3: Calculate Cooperativity (ΔΔE)
To calculate the ΔΔE, we need a 'dde.npy' file, which can be computed using the command provided inside:
# calculate_dde.py. 
This will generate the file 'dde.npy', which is used in the main script cooperativity.py
# Step 4: Calculate Distance Matrix
Please take a look at 'Necessary Files/dist_abl1_2v7a_extra/README.md'  for details.
# Step 5: Run the final script: cooperativity.py
Put all the necessary files in the same directory and run this script, which will produce a histogram along with the 'max_dde_with_pdb.tsv' output file. Finally, the 10 most cooperative double mutant pairs are computed by running the code :

most_cooperative.py
# Note: All necessary data and files are provided inside the folder "Necessary Files". Some files have been provided directly (For file information, please see 'Necessary_Files/README.MD'.


#Code Contributions

The code builds upon original implementations by [Prof. Dr. Allan Haldane](https://scholar.google.com/citations?user=2MBqxWYAAAAJ&hl=en) , one of the co-authors of this manuscript, and incorporates prior modifications by [Dr. Joan Gizzio](https://scholar.google.com/citations?user=D5H_bWEAAAAJ&hl=en) , a former graduate student in our lab. Additional adaptations were made by [GKisan Khatri](https://scholar.google.com/citations?user=IXMrAFsAAAAJ&hl=en) and [Prof. Dr. Allan Haldane](https://scholar.google.com/citations?user=2MBqxWYAAAAJ&hl=en)) for the present analysis. All of these contributions were supervised by [Prof. Dr. Ronald M. Levy](https://scholar.google.com/citations?user=6CQ_uloAAAAJ&hl=en) who provided guidance throughout this project.

# Example Files:
The folder 'Example Files' contains all the necessary files to run the 'cooperativity.py' script for ABL1. It gives an ABL1 histogram. You can use the 'most_cooperative.py' script to identify the 10 most cooperative pairs for ABL1. Please see step 5.

Thank you for the attention.
