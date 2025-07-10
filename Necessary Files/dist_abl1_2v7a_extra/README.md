Following are the information about how we can calculate dist.npy.

./getAlignedContactMap.py 2v7a.pdb ".$(cat ABL1.seq)"  --oSCNHA dist.npy

# in python:
#         from scipy.spatial.distance import squareform
#         d = squareform(np.load('dist.npy'))
#         d[240,260]  # distance from 240 to 260



getAlignedContactMap.py mypdb.pdb A kelvlaLYDYQEKSPREVTMKKGDILTLLNSTNKDWWKVEVNdRQGFVPaayvkkld --oSCNHA dist.npy



Once you have distances.npy, you can then plot it in python as:
"""
d = np.load('distances.npy')
from scipy.spatial.distance import squareform
plt.imshow(squareform(d < 6), cmap='gray_r', origin='lower')
