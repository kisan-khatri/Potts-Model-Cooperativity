This repository provides all relevant code and data to ensure full reproducibility of our results, as requested by the reviewer. We aim to support transparency and allow others to independently validate our findings.
We used a Potts model of Kinase Family proteins to compute a histogram of ΔΔE, specifically (|min(ΔΔE)|), and to identify the most cooperative mutation pairs for two Kinase Proteins, ABL1 and PIM1 as presented in our  AlphaFold Manuscript Entitled: "Using AlphaFold2 to Predict the Conformations of Side Chains in Folded Proteins". To reproduce our results, please follow the following steps.

Thi following information gives details about all the necessary steps.

### 1. Key Files and Tools

-MSA file : realigned259_aln {(used to construct the Potts model (J.npy- coupling file)}
- Couplings: 'J.npy' — obtained using Mi3GPU for kinase protein (We used the Mi3-GPU software package to train the Potts model in our study which is publicly available via its GitHub repository: https://github.com/ahaldane/Mi3-GPU. Please See details at point ###3).
- MSA and Sequence Files: ABL1.seq, PIM1.seq
- get_fitness.py
- Alignment Mapping file: alnMap.npz
- Cooperativity Data: dde.npy (make sure you use this file for the respective protein I have provided or you can compute for ABL1 and PIM1)
 The dde.npy file is required for the double mutant cycle analysis and can be computed using the following command:
 
 Python get_fitness.py  J.npy ABL1.seq  --mode 0 (make sure you have installed the Mi3-GPU software and necessary packages)

- Distance Matrix: dist.npy ((make sure you use this file for respective protein)
- Contact Frequency Data: freq1_IDs.npy
- Bivariate Marginals: bimarg0.4.npy
- PDB Structures: 2V7A_A.pdb, PIM1_A.pdb
-Final script to run: python cooperativity.py

### 2. How to Use
1. Run 'cooperativity.py' to generate the histogram of double mutant cycles.
    Output file: max_dde_with_pdb.tsv
2. Use max_dde_with_pdb.tsv to extract the top 10 most cooperative double mutations (using the script 'most_cooperative.py')
3. For details on how the distance matrix is computed, see the corresponding README.md in this folder: dist_abl1_2v7a_extra.

### 3. Reference for Fitting the Potts Model
Mi3GPU was used for Potts model inference:
We used the Mi3-GPU software package to train the Potts model in our study which is publicly available via its GitHub repository: https://github.com/ahaldane/Mi3-GPU. The details of the software can be found in the paper:

--> Haldane, Allan, and Ronald M. Levy. "Mi3-GPU: MCMC-based inverse Ising inference on GPUs for protein covariation analysis." *Computer Physics Communications* 260 (2021): 107312.

### 4. Notes
- All necessary files for reproducing our results are included.
- Scripts used to generate some files are not included, as the final outputs are provided directly.
- Make sure you are using the correct files for respective proteins (pdb file, distance matrix, bivariate marginal file,dde.npy file and so on (check files needed in 'cooperativity.py').

The example files are uploaded in Example Files* for a quick test.



For Citation: doi: https://doi.org/10.1101/2025.02.10.637534 
      Khatri, Kisan, Ronald M. Levy, and Allan Haldane. "Phylogenetic Corrections and Higher-Order Sequence Statistics in Protein Families: The Potts Model vs MSA Transformer." arXiv preprint arXiv:2503.00289 (2025).


We hope this helps with reproducing our results. Feel free to contact us for any inquiries.
Email: kisan.khatri@temple.edu
       kisankhatri11@gmail.com
