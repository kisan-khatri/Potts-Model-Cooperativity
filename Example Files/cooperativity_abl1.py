from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from numpy import load, array, arange, argwhere, argsort, mean, std, nanmax

# File paths
dde_file = 'dde_abl1.npy'
sequence_file = 'ABL1.seq'
alnMap_file = 'alnMap.npz'
bimarg_file = 'bimarg0.4.npy'
dist_file = 'dist_abl1.npy'
c_freq = 'freq1_IDs.npy'

# Load the DDE data
dde = -load(dde_file)  # change the sign convention so that positive DDE means higher cooperativity

# Load the sequence
with open(sequence_file, 'r') as f:
    seq = f.read().strip()
seq = array(list(seq))
L = len(seq)

# Generate pairs of indices
pairs = array([(i, j) for i in range(L-1) for j in range(i+1, L)])
assert(len(pairs) == dde.shape[0])  # Check to ensure that the length of the sequence is compatible with the shape of dde

# Filter out pairs containing gap characters
pairmsk = (seq[pairs[:, 0]] != '-') & (seq[pairs[:, 1]] != '-')
pairinds = arange(len(pairs))[pairmsk]

# Load the alignment map and PDB numbering
alnmap = load(alnMap_file)
pdb_numbering = alnmap['2V7A_A']

# Load distance and Bimarg data
dist = squareform(np.load(dist_file))
bim = np.load(bimarg_file).reshape((-1, 21, 21))
cf = np.load(c_freq)
ijinds = np.zeros((L, L), dtype=int)
ijinds[np.triu_indices(L, 1)] = np.arange(L*(L-1)//2)
ijinds = ijinds + ijinds.T

# Amino acid codes
alpha = list('-ACDEFGHIKLMNPQRSTVWY')

mutations = []
max_dde = []
distances = []
bimarg_values = []
con_freq = []

for n in pairinds:
    i = pairs[:, 0][n]
    j = pairs[:, 1][n]
    if np.abs(i-j) <= 6:
        continue
    i_pdb = pdb_numbering[i].decode('utf-8')
    j_pdb = pdb_numbering[j].decode('utf-8')
    dde_ind = argwhere(dde[n, 1:, 1:] == nanmax(dde[n, 1:, 1:]))  # dde[n,1:,1:] excludes double mutant cycles containing gap characters

    mutant_i = alpha[dde_ind[0][0] + 1]  # add 1 to the index to account for the offset after excluding gap characters from dde_ind
    mutant_j = alpha[dde_ind[0][1] + 1]
    if any(x in "GA" for x in [mutant_i, mutant_j, seq[i], seq[j]]):
        continue
    mutation_string_i = f"{seq[i]}{i_pdb}{mutant_i}"
    mutation_string_j = f"{seq[j]}{j_pdb}{mutant_j}"
    mutation_string_combined = f"{mutation_string_i}/{mutation_string_j}"
    mutations.append(mutation_string_combined)
    max_dde.append(dde[n, dde_ind[0][0] + 1, dde_ind[0][1] + 1])

    # Calculate distance and Bimarg value and contact frequencies
    d = dist[i, j]
    distances.append(d)
    c_f = cf[i, j]
    con_freq.append(c_f)

    a = alpha.index(mutant_i)
    b = alpha.index(mutant_j)
    if i < j:
        bi = bim[ijinds[i, j], a, b]
    else:
        bi = bim[ijinds[j, i], b, a]
    bimarg_values.append(bi)

max_dde = array(max_dde)
mutations = array(mutations)
distances = array(distances)
bimarg_values = array(bimarg_values)
con_freq = array(con_freq)
srt = argsort(max_dde)[::-1]  # Sort the cooperativities from largest to smallest

# Write to output file
with open('max_dde_with_pdb.tsv', 'w') as f:
    print('Mutation\tCooperativity\tDistance\tBimarg\tCF', file=f)
    for m, e, d, bi, cof in zip(mutations[srt], max_dde[srt], distances[srt], bimarg_values[srt],con_freq[srt]):
        print(f"{m}\t{round(e, 2)}\t{round(d, 2)}\t{100 * bi:5.2f}\t{round(cof, 2)}", file=f)



# Plot histogram
plt.hist(max_dde, bins=100, density=True)
plt.ylabel('Frequency')
plt.xlabel(r'$|min(\Delta \Delta E)|$')
plt.ylim(0, 400*0.007)
plt.xlim(0, 6)
plt.legend([f'Avg={mean(max_dde):.3f}\nStd={std(max_dde):.3f}'])
plt.title("$\Delta \Delta E$ for ABL1 (PDB 2V7A)")

ax2 = plt.gcf().add_axes([0.3, 0.3, 0.5, 0.4])
plt.hist(max_dde, bins=100, density=True)
plt.ylim(0, 3*0.007)
plt.xlim(1, 6)

#for i in  np.argsort(max_dde)[-10:]:
#    plt.text(max_dde[i], 0, mutations[i], rotation='vertical')

plt.show()

