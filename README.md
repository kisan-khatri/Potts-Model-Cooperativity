# Using AlphaFold2 to Predict the Conformations of Side Chains in Folded Proteins
DOI: https://doi.org/10.1101/2025.02.10.637534

This repository contains the code and data used in our study on double mutant cycle (cooperativity) analysis in kinase proteins, as presented in the manuscript. It includes the code to reproduce the histogram shown in the manuscript and to identify the most cooperative double mutant pairs reported. The repository is shared in response to reviewer requests to promote transparency and reproducibility. Permission to share these materials was granted by my advisor, Prof. Dr. Ronald M. Levy [Google Scholar](https://scholar.google.com/citations?user=6CQ_uloAAAAJ&hl=en), who is also the corresponding author of the manuscript.

We used a Potts model of kinase protein families to compute a histogram of double mutation energies (specifically, |min(ΔΔE)|) and to identify the most cooperative mutation pairs in two kinases: ABL1 and PIM1.
# Instruction

For the detailed flow path, please refer to the README.md file located inside the folder Necessary_Files. However, I am providing a brief summary of the steps below.
# Step 1: Software Used
To perform this analysis, we used the [Mi3-GPU software](https://github.com/ahaldane/Mi3-GPU), which is publicly available on GitHub.  
Follow the installation and usage instructions provided in the `README.md` file inside the folder `Necessary_Files`.
# Step 2: Construct Kinase MSA and Fit the Potts Model
The multiple sequence alignment (MSA) used in this study is `realigned259_aln`, and the fitted Potts model coupling matrix is provided as `J.npy`. These files, along with all other necessary data, scripts, and documentation, are available in the folder `Necessary_Files`.  
# Step 3: Calculate Double Mutant Cycle Energies
To calculate the double mutant cycle energy differences (ΔΔE), run the script provided inside:
# calculate_dde.py. 
This will generate the file `dde.npy`, which contains the ΔΔE values used for downstream cooperativity analysis.
# Step 4: Calculate Distance Matrix
Please refer to 'Necessary Files/dist_abl1_2v7a_extra/README.md  for details.
# Step 5: Run the final script: cooperativity.py
Put all the necessary file in the same directory and run this script which will produce histogram along with 'max_dde_with_pdb.tsv'. Finally, 10 most cooperative double mutant pairs are computed by running code :

most_cooperative.py
# Note: All necessary data and files are provided inside folder "Necessary Files". Some files have been provided directly (For file information, please see Necessary_Files/README.MD.


#Code Contributions

The code builds upon original implementations by Prof. Dr. Allan Haldane (co-author of this manuscript) [Google Scholar](https://scholar.google.com/citations?user=2MBqxWYAAAAJ&hl=en) and incorporates prior modifications by Joan Gizzio [Google Scholar](https://scholar.google.com/citations?user=D5H_bWEAAAAJ&hl=en) , a former graduate student in our lab. Additional adaptations were made by Kisan Khatri and Prof. Dr. Allan Haldane for the present analysis.

# Example Files:
The folder, Example Files contains all the neccessar files to run cooperativity.py script for ABL1. It gives ABL1 histogram. You can use most_cooperative.py script to identify 10 most cooperative pairs for ABL1. Please see step 5.
