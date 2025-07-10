#!/usr/bin/env python3
import numpy as np
from Bio import Align
import pdbparse
from scipy.spatial.distance import cdist
import sys, os, argparse

code = pdbparse.residueCode
# some PDBs have nonstandard residues. PDBs are supposed to list them in MODRES
# records, but don't always. So we keep a "default" list of nonstandard
# residues which is overridden if MODRES is found (important since different
# PDBs sometimes assign the same name to different residues)
extracode = {'ABA': 'GLU', 'B3L': 'LEU', 'CAS': 'CYS', 'CIR': 'ARG', 'CME':
             'CYS', 'CMT': 'CYS', 'CSD': 'CYS', 'CSO': 'CYS', 'CSS': 'CYS',
             'CY0': 'CYS', 'FCL': 'PHE', 'KCX': 'LYS', 'L3O': 'LEU', 'LGY':
             'LYS', 'MEA': 'PHE', 'MHO': 'MET', 'MSE': 'MET', 'NMM': 'ARG',
             'OCS': 'CYS', 'OCY': 'CYS', 'PFF': 'TYR', 'PTR': 'TYR', 'SCS':
             'CYS', 'SEP': 'SER', 'TPO': 'THR', 'TYI': 'TYR'}
code.update(dict((k, code[v]) for k,v in extracode.items()))

pairs = lambda L: ((i,j) for i in range(L-1) for j in range(i+1,L))

def get_distances(c, CA):
    resids = CA.resid
    L = len(resids)

    # compute Nearest-heavy-atom (NHA) distances
    atoms = c[c.element != 'H']
    atomids = atoms.resid
    cc = [atoms.coords[np.argwhere(atomids == id).ravel()] for id in resids]
    NHAdist = np.array([min(cdist(cc[i], cc[j]).ravel()) for i,j in pairs(L)])

    # compute side-chain Nearest-heavy-atom (SCNHA) distances
    atoms = c[(c.element != 'H') & (c.atom != ' C  ') &
              (c.atom != ' O  ') & (c.atom != ' N  ') &
              ((c.atom != ' CA ') | (c.resname == 'GLY'))]
    atomids = atoms.resid
    cc = [atoms.coords[np.argwhere(atomids == id).ravel()] for id in resids]
    cc = [ci if ci.size != 0 else [[nan,nan,nan]] for ci in cc]
    SCNHAdist = np.array([min(cdist(cc[i], cc[j]).ravel()) for i,j in pairs(L)])

    # compute side-chain center (SCC) distances
    scc = np.array([np.mean(ci, axis=0) for ci in cc])
    SCCdist = cdist(scc, scc)[np.triu_indices(L,k=1)]

    # compute CA distance
    CAdist = cdist(CA.coords, CA.coords)[np.triu_indices(L,k=1)]

    # compute CB distance.  Count CA as CB for gly
    CB = c[(c.atom == ' CB ') |
           ((c.atom == ' CA ') & (c.resname == 'GLY'))]
    inds = [np.argwhere(CB.resid == id).ravel() for id in resids]
    cc = np.array([CB.coords[ind[0]] if len(ind) > 0 else [nan,nan,nan]
                  for ind in inds])
    CBdist = cdist(cc, cc)[np.triu_indices(L,k=1)]

    # sanity checks
    assert(len(CAdist) == len(CBdist))
    assert(len(CAdist) == len(NHAdist))
    # sanity check for clashes
    clashes = [(resids[i], resids[j]) for n,(i,j) in enumerate(pairs(L))
               if abs(i-j) > 4 and NHAdist[n] < 2.0]
    # 2.0 isdistance at which pymol draws bonds
    if len(clashes) > 0:
        print('Warning: clashes in residues {}'.format(clashes))

    return NHAdist, SCCdist, SCNHAdist, CAdist, CBdist

def print_alignment(seq, alnseq, alnpdb):

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    # First, we need to put inserts and gaps back in
    cseq = []
    seq_ind = 0
    for si in alnseq:
        if si == '-':
            cseq.append('.')
            continue
        oc = seq[seq_ind]
        if oc != '-':
            assert(oc.upper() == si)
        cseq.append(oc)
        seq_ind += 1
    cseq = "".join(cseq)

    print("")
    print('Alignment:   ("-" or uppercase in Seq included in contact map)')
    print("")
    for a,b in zip(chunks(cseq, 70), chunks(alnpdb, 70)):
        #print("Seq: ", aln[0])
        #print("PDB: ", aln[1])
        print("Seq: ", a)
        print("PDB: ", b)
        print("")

def alignedContacts(pdb, chain=None, seq=None):

    if chain is None:
        chain = pdb.chain[0]

    #get structure sequence from CA atoms
    c = pdb[(pdb.chain == chain)]
    c = c[((c.altloc == ' ') | (c.altloc == 'A'))]
    CA = c[c.atom == ' CA ']

    if len(CA) == 0:
        raise ValueError("No residues in chain {}".format(chain))

    # get letters and throw out unkown residue types (usually HET entries)
    CAseq = np.array([code.get(r, 'Z') for r in CA.resname])
    keep = CAseq != 'Z'
    CA = CA[keep]
    CAseq = CAseq[keep]

    # process PDB distances
    dists = get_distances(c, CA)
    #print('dists', dists)

    if seq is None:
        print("No alignment specified, using full pdb sequence:")
        print("PDB: ", "".join(CAseq))
        return dists

    inserts = [c.islower() for c in seq]
    # use "X" to distinguish parwise2 gaps from original gaps
    useq = seq.replace('-', 'X').upper() 
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1.5
    aligner.extend_gap_score = -0.1
    aln = aligner.align(useq, "".join(CAseq))[0]

    # print out alignment for instpection (sanity check)
    print_alignment(seq, aln[0], aln[1])
    
    # build up list of structural indices which map to MSA positions
    # could modify this to account for inserts in seq
    struct_ind = 0
    seq_ind = 0
    caind = []
    for si, cai in zip(aln[0], aln[1]):
        if si != '-':
            if not seq[seq_ind].islower():
                caind.append(struct_ind if (cai != '-' and si != 'X') else None)
            seq_ind += 1
        if cai != '-':
            struct_ind += 1
    print("Aligned length:", len(caind))

    # now get aligned subset of contact map
    L = len(caind)
    pdbL = len(CAseq)
    pdb_pairs = dict((ij,n) for n,ij in enumerate(pairs(len(CAseq))))
    aln_distmap = np.array([pdb_pairs[(i,j)] if i != None and j != None else -1
                         for i,j in ((caind[a], caind[b]) for a,b in pairs(L))],
                         dtype=int)
    nix = aln_distmap < 0

    out = []
    for d in dists:
        do = d[aln_distmap]
        do[nix] = np.nan
        out.append(do)
    return out

def getM(x, diag_fill=0):
    L = int(((1+np.sqrt(1+8*len(x)))/2) + 0.5)
    M = np.zeros((L,L), dtype=x.dtype)
    i, j = np.triu_indices(L,k=1)
    M[i,j] = x
    M[j,i] = x
    M[np.diag_indices(L)] = diag_fill
    return M

def main():
    parser = argparse.ArgumentParser(description="Compute a contact map",
    epilog=("Given a PDB and a sequence, this computes a contact map "
            "for the parts of the PDB sequence that align to the input "
            "sequence (or the whole sequence if no sequence is provided). "
            "The input sequence may contain gaps ('-'), and inserts "
            "in lowercase. The output contacts only correspond to pairs of "
            "positions which are gap or uppercase in the input sequence. "
            "Output is .npy file containing upper-triangle of contact map."
            "If no chain is supplied, the first chain will be used."
            ))
    parser.add_argument('pdbfile', help='PDB file')
    parser.add_argument('seq', nargs="?", help="sequence to align to")
    parser.add_argument('--chain', help='PDB chain id')
    parser.add_argument('--oNHA',
        help='output nearest-heavy-atom distance file')
    parser.add_argument('--oSCC',
        help='output side-chain-center distance file')
    parser.add_argument('--oSCNHA',
        help='output side-chain-nearest-heavy-atom distance file')
    parser.add_argument('--oCA',
        help='output C-alpha distance file')
    parser.add_argument('--oCB',
        help='output C-beta distance file')
    parser.add_argument('--plot', nargs='?', const='',
        help='make a plot of the last output distance file, '
             'using optional supplied distance cutoff')

    args = parser.parse_args(sys.argv[1:])

    pdbname = os.path.splitext(os.path.basename(args.pdbfile))[0]
    pdb = pdbparse.loadPDB(args.pdbfile)
    seq = args.seq.replace('.', '') if args.seq else None

    dists = alignedContacts(pdb, args.chain, seq)
    NHAdist, SCCdist, SCNHAdist, CAdist, CBdist = dists

    toplot = None
    if args.oNHA:
        np.save(args.oNHA, NHAdist)
        toplot = ('NHA', NHAdist, 6.0)
    if args.oSCC:
        np.save(args.oSCC, SCCdist)
        toplot = ('SCC', SCCdist, 8.0)
    if args.oSCNHA:
        np.save(args.oSCNHA, SCNHAdist)
        toplot = ('SCNHA', SCNHAdist, 8.0)
    if args.oCA:
        np.save(args.oCA, CAdist)
        toplot = ('CA', CAdist, 8.0)
    if args.oCB:
        np.save(args.oCB, CBdist)
        toplot = ('CB', CBdist, 8.0)

    if args.plot is not None:
        if toplot is None:
            raise ValueError("Please specify a score to compute, eg --SCNHA")

        name, dists, cutoff = toplot

        if args.plot != '': # use default cutoff
            cutoff = float(args.plot)
        import pylab as plt
        plt.imshow(getM(dists < cutoff), origin='lower',
                interpolation='nearest', cmap='gray_r')
        plt.title(f"{name} < {cutoff}")
        plt.show()

if __name__ == '__main__':
    main()
