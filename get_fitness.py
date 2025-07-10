import numpy as np
from numpy.random import randint
import argparse, sys
# Adapted from original code by Allan Haldane (one of the coauthor of this manuscript) and Joan Gizzio (Former Graduate student of our Lab). Our advisor is Prof. Dr. Ronald M. Levy (Corresponding author of this manuscript).

def findall(cond):
    return np.argwhere(cond).flatten()

def getLq(J):
    L = int(((1+np.sqrt(1+8*J.shape[0]))/2) + 0.5)
    q = int(np.sqrt(J.shape[1]) + 0.5)
    return L, q

def energy(s, J):
    "Compute energy of sequence s"
    L, q = getLq(J)
    pairs = np.array([s[j] + q*s[i] for i in range(L-1) for j in range(i+1,L)])
    return np.sum(J[np.arange(len(pairs)), pairs])

def expandJ(J):
    "Convert J matrix from shape (L*(L-1)//2, q*q) to (L, L, q, q)"
    L, q = getLq(J)
    Jex = np.zeros((L, L, q, q))
    J = J.reshape((J.shape[0], q, q))
    for n,(i,j) in enumerate((i,j) for i in range(L-1) for j in range(i+1,L)):
        Jex[i,j,:,:] = J[n,:,:]
        Jex[j,i,:,:] = J[n,:,:].T
    return Jex

def dE(s, Jex):
    "Compute all single-mutant \Delta E for a sequence"
    L, q = Jex.shape[0], Jex.shape[2]
    Lrng = np.arange(L)
    dEs = np.empty((L,q), dtype=float)
    for i in range(L):
        # JJ is L x q matrix, where JJ[l,q] gives coupling J^{il}_{q,sl}
        JJ = Jex[i,Lrng,:,s]  
        dEs[i,:] = np.sum(JJ - JJ[:,s[i],None], axis=0)
    return dEs

def ddE(s, J):
    "Compute all \Delta \Delta E for a sequence"
    L, q = getLq(J)

    si, sj = zip(*((s[i],s[j]) for i in range(L-1) for j in range(i+1,L)))
    si, sj = np.array(si), np.array(sj)
    
    n = J.shape[0]
    J = J.reshape((n, q, q))
    N = np.arange(n)

    return ( + J
             - J[N,:,None,sj] 
             - J[N,None,si,:] 
             + J[N,None,si,None,sj] )

def get_all_dE(s, J):
    """
    Returns, for a sequence: 
        1. all \Delta\Delta E, 
        2. all double-mut \Delta E, 
        3. all single-mutant \Delta E
    """
    Jex = expandJ(J)
    L, q = getLq(J)

    ddEs = ddE(s, J)
    dEs = dE(s, Jex)
    dEab = np.array([np.add.outer(dEs[i,:], dEs[j,:])
                     for i in range(L-1) for j in range(i+1,L)])
    return ddEs, ddEs + dEab, dEs

def sanity_check(s, J, alldE, alpha):
    ddEs, dE2s, dEs = alldE

    L, q = getLq(J)
    N = J.shape[0]

    # choose two random positions
    rn = randint(N)
    pairs = ((i,j) for i in range(L-1) for j in range(i+1,L))
    i, j = ((i,j) for n,(i,j) in enumerate(pairs) if n == rn).next()
    # choose two random mutations
    qi, qj = randint(q), randint(q)
    
    # construct mutated sequences
    si, sj, sij = s.copy(), s.copy(), s.copy()
    si[i] = qi
    sj[j] = qj
    sij[i], sij[j] = qi, qj
    
    # compute sequence energies
    E = energy(s, J)
    Ei = energy(si, J)
    Ej = energy(sj, J)
    Eij = energy(sij, J)
    
    print("Testing mutation at ({}, {}) to  {},{}".format(
                                                    i, j, alpha[qi], alpha[qj]))
    print("             Expected      Computed")
    print("ddE     :  {:10.3f}    {:10.3f}".format(Eij-Ei-Ej+E, ddEs[rn,qi,qj]))
    print("dE ij   :  {:10.3f}    {:10.3f}".format(Eij - E, dE2s[rn,qi,qj]))
    print("dE i    :  {:10.3f}    {:10.3f}".format(Ei - E, dEs[i,qi]))
    print("dE j    :  {:10.3f}    {:10.3f}".format(Ej - E, dEs[j,qj]))
    print("")
    print("Energies:")
    print("E  : {:10.3f}".format(E))
    print("Ei : {:10.3f}".format(Ei))
    print("Ej : {:10.3f}".format(Ej))
    print("Eij: {:10.3f}".format(Eij))

def main():
    parser = argparse.ArgumentParser(description='Compute Sequence Energies from Couplings')
    parser.add_argument('J', help='Required: Potts model couplings. Make sure you are using the appropriate gauge for your calculations.')
    parser.add_argument('seq', help='Required: a file containing a single (aligned) sequence of length L, supplied as a string')
    # Note the addition of '--' before the argument name to make it optional
    parser.add_argument('--contacts', default=None, help="Optional: LxL numpy matrix of contact scores. If unspecified, all pairs of positions may contribute to the computed energies (classic 'fitness' type calculation). Otherwise, the user provides a .npy matrix of contact values. Boolean True/False or 1/0 values indicate the presence or absence of a contact, usually for a particular structure. Floating point values in the range [0, 1] are interpreted as contact frequencies, usually for an ensemble of structures. If contact frequency *differences* are used, usually between a pair of structures or two ensembles of structures, the values may be in the range [-1, 1].")
    parser.add_argument('--mode', type=int, choices=[0, 1, 2, 3], default=0, help="Specifies the type of energy calculation to perform (relative to the wildtype) for the positions listed in '-positions'. 0: compute all possible single and double mutant energy differences and the double-mutant cycle cooperativites. 1: single mutant energy differences only. 2: double mutant energy differences only. 3: double-mutant cycle cooperativities only.")
    parser.add_argument('--positions', default=None, help="Optional: a list of one or more comma-separated positions in the range [1, L], which are used to limit the output of the mutational scan. By default all positions are included in the mutational scan, and the only output will be in .npy format. Alternatively, if one or more positions are specified the output will be a text file. Note - the '-positions' argument only excludes positions from being mutated, which limits the output. It will not exclude any positions from contributing to the computed energies (for that you need to supply an argument for '-contacts').")
    parser.add_argument('--mutations', default=None, help="Optional: a list of one or more amino acid letters to limit the mutational scan. By default positions are mutated to all possible letters.")
    
    args = parser.parse_args()

#    args = parser.parse_args(sys.argv[1:])

    J = np.load(args.J)
    L, q = getLq(J)
    pairs = np.array([(i,j) for i in range(L-1) for j in range(i+1, L)]) #array of all possible position pairs, dim0 has size L*(L-1)/2, dim1 has size of 2
    #alpha = args.alpha.strip()
    assert(pairs.shape[0] == (L * (L-1)/2))
    assert(pairs.shape[1] == 2)
    pairinds = np.arange(pairs.shape[0])
    alpha = list('-ACDEFGHIKLMNPQRSTVWY') #all 20 amino acids plus a gap character (q=21)
    assert(q == len(alpha))
    alphainds = np.arange(len(alpha))
    alphaindpairs = np.array([(a,b) for a in alphainds for b in alphainds]) #get all 21*21 possible combinations of letter pairs 
    seq = open(args.seq, 'r').readlines()[0].strip()
    print(seq)
    seq = np.array([alpha.index(c) for c in seq], dtype='u1')
    contactmode = args.contacts
    if args.contacts:
        contacts = np.load(args.contacts)[np.triu_indices(L, k=1)]
        contacts = np.nan_to_num(contacts)
        es = (J.T * contacts).T
    else:
        es = J

    ddEs, dE2s, dEs = get_all_dE(seq, es)

    mutations = [] #create a list of mutations, if any are specified, to limit subsequent output
    for n,a in enumerate(alpha):
        if args.mutations:
            if a in args.mutations:
                mutations.append(n)
        else:
            mutations.append(n) #if the user hasn't supplied mutations with the -mutations argument, by default we will mutate to all possible amino acids
    mutations = np.array(mutations)
   
    if args.mode is None:
        mode = 0
    else:
        mode = int(args.mode)
    if mode == 0:
        if not args.positions:
            np.save('dde', ddEs)
            np.save('de', dEs)
            np.save('de2', dE2s)

    if mode == 1:
        if not args.positions:
            np.save('de', dEs)
        else:
            with open('output_de.tsv', 'w') as f:
                headers = ['i']
                for a in alphainds:
                    if a not in mutations:
                        continue
                    else:
                        headers.append(alpha[a])
                print('\t'.join(headers), file=f)
                for pos in args.positions.split(','): #loop through positions provided by the -positions argument
                    values = [] #keep track of the energy values for each amino acid at this position
                    i = int(pos)
                    for a in alphainds: #loop through all 20 amino acids (and a gap character) to accumulate mutational data
                        if a not in mutations: #check if the amino acid is missing from our list of mutations, to limit the output
                            continue
                        else:
                            values.append(str(round(dEs[i - 1,a],3))) #Important: Numpy indexes positions like [0,L-1] which is why I have supplied "i-1" as the index. For example, the first position in the MSA which, by conventioned, we refer to as "position 1", actually has an index of 0 in the numpy array. This offset in indexing needs to be accounted for carefully.
                    print('{}\t{}'.format(i, '\t'.join(values)), file=f) #print the accumulated mutational data for this position to the output file

    if mode == 2:
        if not args.positions:
            np.save('de2', dE2s)
        else:
            with open('output_de2.tsv', 'w') as f:
                headers = ['i', 'j']
                for (a,b) in alphaindpairs:
                    if (a not in mutations) or (b not in mutations):
                        continue
                    else:
                        headers.append('{}{}'.format(alpha[a], alpha[b]))
                print('\t'.join(headers), file=f)
                pos = np.array(args.positions.split(',')).astype(int)
                print(pos)
                ijs = np.array([(i, j) for n,i in enumerate(pos[:-1]) for j in pos[n+1:]]) #create a list of all possible pairs that can be made from the provided positions
                print(ijs)
                for pair in ijs: #loop through position pairs (i,j)
                    values = [] #keep track of the energy values for each pair of amino at this pair of positions
                    i = pair[0]
                    j = pair[1]
                    pairind = findall((pairs[:,0]==(i-1)) & (pairs[:,1]==(j-1)))[0] #get the index of this specific residue pair, so we can find it in the pairinds array defined earlier.
                    for (a,b) in alphaindpairs: #loop through all 21 * 21 possible amino acid pairs
                        if (a not in mutations) or (b not in mutations): #check if either mutation is missing from our list of mutations, to limit the output
                            continue
                        else:
                            values.append(str(round(dE2s[pairind, a, b],3)))
                    print('{}\t{}\t{}'.format(i, j, '\t'.join(values)), file=f) #print the mutational data for this position to the output file

    if mode == 3:
        if not args.positions:
            np.save('dde', ddEs)
        else:
            with open('output_dde.tsv', 'w') as f:
                headers = ['i', 'j']
                for (a,b) in alphaindpairs:
                    if (a not in mutations) or (b not in mutations):
                        continue
                    else:
                        headers.append('{}{}'.format(alpha[a], alpha[b]))
                print('\t'.join(headers), file=f)
                pos = np.array(args.positions.split(',')).astype(int)
                ijs = np.array([(i, j) for n,i in enumerate(pos[:-1]) for j in pos[n+1:]]) #create a list of all possible pairs that can be made from the provided positions
                for pair in ijs: #loop through position pairs (i,j)
                    values = [] #keep track of the energy values for each pair of amino at this pair of positions
                    i = pair[0]
                    j = pair[1]
                    pairind = findall((pairs[:,0]==(i-1)) & (pairs[:,1]==(j-1)))[0] #get the index of this specific residue pair, so we can find it in the pairinds array defined earlier.
                    for (a,b) in alphaindpairs: #loop through all 21 * 21 possible amino acid pairs
                        if (a not in mutations) or (b not in mutations): #check if either mutation is missing from our list of mutations, to limit the output
                            continue
                        else:
                            values.append(str(round(ddEs[pairind, a, b],3)))
                    print('{}\t{}\t{}'.format(i, j, '\t'.join(values)), file=f) #print the mutational data for this position to the output file

    for i in range(0):
        print("")
        print("Sanity check #{}:".format(i))
        sanity_check(seq, J, (ddEs, dE2s, dEs), alpha)
    
if __name__ == '__main__':
    main()

#python get_dde.py J.npy ..ABL1.seq --mode 0
