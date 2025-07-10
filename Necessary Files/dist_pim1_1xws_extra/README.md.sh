./getAlignedContactMap.py 2v7a.pdb ".$(cat ABL1.seq)"  --oSCNHA dist.npy

# in python:
#         from scipy.spatial.distance import squareform
#         d = squareform(np.load('dist.npy'))
#         d[240,260]  # distance from 240 to 260

