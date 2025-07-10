#!/usr/bin/env python
import numpy as np
import sys, os, pickle
from numpy import rec

"""
ATOM      8  C  BASP A   2      76.283  57.361  60.309  0.50 84.80           C
ATOM  44035  H2  TIP3 9165     -42.898 -14.686  13.233  1.00  0.00      WT8
ATOM   3949  OXT ARG E 118       1.647   7.647  53.839  1.00 68.00           O
AAAAAAIIIII AAAA RRRRCNNNN    XXXXXXXXYYYYYYYYZZZZZZZZOOOOOOBBBBBB      SSSSEECC
                    ?     ?                                             ????
12345678901234567890123456789012345678901234567890123456789012345678901234567890
0        1         2         3         4         5         6         7         8

#I'm not quite consistent this table, see ? above.
#note: This is in base-1 indexing, but my arrays are base-0

COLUMNS      DATA TYPE        FIELD      DEFINITION
------------------------------------------------------
 1 -  6      Record name      "ATOM    "
 7 - 11      Integer          serial     Atom serial number.
13 - 16      Atom             name       Atom name.
17           Character        altLoc     Alternate location indicator.
18 - 20      Residue name     resName    Residue name.
22           Character        chainID    Chain identifier.
23 - 26      Integer          resSeq     Residue sequence number.
27           AChar            iCode      Code for insertion of residues.
31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angs.
39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angs.
47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angs.
55 - 60      Real(6.2)        occupancy  Occupancy.
61 - 66      Real(6.2)        tempFactor Temperature factor.
77 - 78      LString(2)       element    Element symbol, right-justified.
79 - 80      LString(2)       charge     Charge on the atom.

the atom name field (4 chars) is further subdivided: 1st two chars are element
(right justified, so " C  " is carbon, "CA  " is calcium. Next char is
'distance indicator', named by greek letters (ABGDEZH), so " CA " is an alpha
carbon. Last char is branch number, so eg " HA3" is the third hydrogen
attatched to an atom. Hydrogens are named specially: The 'distance indicator'
is labelled by which carbon the H is attached to, and if there are many H
attached to an atom, they are given names "1H", "2H" etc as the elemtn name, so
eg "1HA3" would be H #1 attached to alpha carbon, 1st branch (If that is even
possible).
"""

residueCode = { "GLY": "G", "PRO": "P", "ALA": "A", "VAL": "V", "LEU": "L",
                "ILE": "I", "MET": "M", "CYS": "C", "PHE": "F", "TYR": "Y",
                "TRP": "W", "HIS": "H", "LYS": "K", "ARG": "R", "GLN": "Q",
                "ASN": "N", "GLU": "E", "ASP": "D", "SER": "S", "THR": "T" }

weights = {' C': 12.01, ' H': 1.01,  ' O': 16.00, ' P': 30.97,
           ' N': 14.01, 'Na': 22.99, 'Cl': 35.45, ' S': 28.09}

# to speed up loading time, if a directory 'pickedPDBs' exists the data will
# be stored in pickled form in that directory, so it can be loaded faster the
# next time.
def loadPDB(filename):
    #in principle more careful error checking is needed but good enough for now
    if(os.path.exists('./pickledPDBs')):
        pickledName = os.path.join('./pickledPDBs', filename + '.pkl')

        if( os.path.exists(pickledName) ):
            with open(pickledName, 'rb') as f:
                pdb = pickle.load(f)
            return pdb
        else:
            print("creating pickle file")
            pdb = loadTextPDB(filename)
            with open(pickledName, 'wb') as f:
                pickle.dump(pdb, f, 2)
            return pdb
    else:
        return loadTextPDB(filename)

#loads a pdb from text. Returns a recarray, with entries as in pdbtype
def loadTextPDB(filename):
    with open(filename, 'r') as f:
        data = f.readlines()

    # note: Use of "U" a bit wasteful since uses utf-32, but is more convenient
    # than "S" since we don't have to encode/decode and add b' everywhere.
    pdbtype = [ ('id',          int),
                ('atom',      ('U4',
                   [('element', 'U2'), ('distance', 'U1'), ('branch', 'U1')])),
                ('altloc',     'U1'),
                ('resname',    'U4'), # longer than standard, stripped
                ('chain',      'U1'),
                ('resid',      'U5'), # longer than standard
                ('resnum',      int), # int version of resid. gro makes it hex?
                ('icode',      'U1'),
                ('coords',    float,3),
                ('occupancy', float),
                ('beta',      float),
                ('segment',    'U4'), # nonstandard, stripped
                ('element',    'U2'), # stripped
                ('charge',     'U2'), # stripped
                ('model',       int), # extra: model number
                ('hetatm',     bool)] # extra: whether a hetatm

    model = 0
    atomLines = []
    for line in data:
        if   line[0:6] == 'ATOM  ':
            atomLines.append((line, model, False))
        elif line[0:6] == 'HETATM':
            atomLines.append((line, model, True))
        elif line[0:6] == 'MODEL ':
            model += 1
        elif line[0:6] == 'ENDMDL':
            break
    
    # some prorgams (eg gro) can write weird resnum
    def tryint(x):
        try:
            return int(x)
        except:
            return -1

    atoms = np.empty(len(atomLines), np.dtype(pdbtype))
    for n,(line,model,hetatm) in enumerate(atomLines):
        atoms[n] = (int(line[6:11]),
                        line[12:16],
                        line[16],
                        line[17:21].strip(),
                        line[21],
                        line[22:27],
                 tryint(line[22:26]),
                        line[26],

                ( float(line[30:38]),
                  float(line[38:46]),
                  float(line[46:54]) ),

                  float(line[54:60]),
                  float(line[60:67]),
                        line[72:76].strip(),
                        line[76:78].strip(),
                        line[78:80].strip(),

                        model,
                        hetatm)

    pdb = rec.array(atoms)
    return pdb

def writeAtoms(atoms, file=None):

    last_model = atoms[0].model
    if atoms[-1].model != atoms[0].model:
        print("MODEL ", file=file)

    last_chain = atoms[0].chain

    for a in atoms:
        if a.model != last_model:
            print("TER   ", file=file)
            print("ENDMDL", file=file)
            print("MODEL ", file=file)
            last_ter = a.ter
        elif last_chain != a.chain:
            print("TER   ", file=file) #todo: should print resn/resi/etc too

        tp = 'ATOM  ' if not a.hetatm else 'HETATM'
        x,y,z = a.coords
        print(f"{tp}{a.id:5d} {a.name:4s}{a.altloc:1s}{a.resname:3s}"
              f"{a.chain:1s}{a.resid:5s}{a.icode:1s}  "
              f"{x:8.3f}{y:8.3f}{z:8.3f}"
              f"{a.occupancy:6.2f}{a.beta:6.2f}      "
              f"{a.segment:4s}{a.element:2s}{a.charge:2s}", file=file)

if __name__ == '__main__':
    loadPDB(sys.argv[1])



