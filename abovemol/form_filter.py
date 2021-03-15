#!/usr/bin/env python3
# coding: utf-8

'''
TODO: description
'''

from os import environ
import sys, getopt
import numpy as np
from MDAnalysis import *

infile = ''
bt_file = ''
tprfile = ''
lipid_residue = ''

# parse command line options
try:
    opts, args = getopt.getopt(sys.argv[1:], 'hi:s:r:o:')
except getopt.error as msg:
    print(msg)
    print('for help use -h')
    sys.exit(2)
# process options
for o, a in opts:
    if o == '-h':
      print('''
./form_filter.py -s sd0.pqr -i input.xtc -r POPC -o polymer [ -x index.ndx:LYS1:LYS2:LYS3  --  not implemented yet]
            ''')
      sys.exit(0)
    elif o == '-i':
        infile = a
    elif o == '-s':
        tprfile = a
    elif o == '-o':
        bt_file = a+'_bt.txt'
        up_file = a+'_up.txt'
    elif o == '-r':
        lipid_residue = a
    # todo:
    # implement -x

if infile == '' or bt_file == '' or tprfile == '' or lipid_residue == '':
    print('for help use -h')
    sys.exit(2)

# now we use explicitly determined molecule borders
# todo: replace with -x option
index_pl_b = np.array([34789, 35232, 70463, 70906, 106137, 106580, 141811, 142254]) - 1
index_pl_e = index_pl_b + 442

# THRESHOLD (Rxy - com_xy)^2 < Rthr^2 == THRESHOLD
# default 400 i.e. Rthr = 2 nm
# TODO: make parameter
THRESHOLD = 400 # 2^2

maxframes = 1e20

un = Universe(tprfile, infile)

natm = un.atoms.n_atoms
natm2 = int(np.ceil(natm/8)*8)

bilayer_atoms = un.select_atoms('resname %s' % lipid_residue)
if bilayer_atoms.n_atoms == 0:
    print('Incorrect bilayer residue!')
    sys.exit(2)
zc = bilayer_atoms.center_of_geometry()[2]

# determine which molecule is in lower ML and which is in upper ML
maskmols = [ un.select_atoms('index %d:%d' % (index_pl_b[i], index_pl_e[i]) ) for i in range(index_pl_b.shape[0]) ]

maskmols_up = []
maskmols_bt = []
for i in range(len(maskmols)):
    mol = maskmols[i]
    if (mol.center_of_mass()[2] > zc):
        maskmols_up.append(mol)
    else:
        maskmols_bt.append(mol)

print("""
Found molecules:
    ABOVE - %d
    UNDER - %d
""" % (len(maskmols_up), len(maskmols_bt)) )

mask_up = np.zeros(un.atoms.n_atoms, dtype=np.bool)
mask_bt = np.zeros(un.atoms.n_atoms, dtype=np.bool)

finished = False
# auxilary variable for encountering PBC
pbcs = np.array([[0,0],[0,1],[1,0],[1,1],[-1,0],[0,-1],[-1,-1],[-1,1],[1,-1]])

with open(up_file, 'wb') as upfd:
 with open(bt_file, 'wb') as btfd:
    co = 0
    while not finished:
        _ts = un.trajectory.ts
        xxx = _ts._pos

        # bilayer center (OZ)
        zc = bilayer_atoms.center_of_geometry()[2]

        up_coms = [mol.center_of_mass() for mol in maskmols_up]
        bt_coms = [mol.center_of_mass() for mol in maskmols_bt]

        mask_up = (xxx[:,2] > zc)
        totnup = np.sum(mask_up)
        mask_bt = np.invert(mask_up)
        totnbt = np.sum(mask_bt)

        for com in up_coms:
            com_ = com
            cmar = np.zeros((natm,2))

            for i in range(pbcs.shape[0]):
                cmar *= 0.0
                com_ = com[0:2] + _ts.dimensions[0:2]*pbcs[i,:]
                cmar[:,0] = com_[0]; cmar[:,1] = com_[1]
                mask_up *= (np.sum( (xxx[:,0:2]-cmar)**2, axis=1) > THRESHOLD )

        for com in bt_coms:
            com_ = com
            cmar = np.zeros((natm,2))

            for i in range(pbcs.shape[0]):
                cmar *= 0.0
                com_ = com[0:2] + _ts.dimensions[0:2]*pbcs[i,:]
                cmar[:,0] = com_[0]; cmar[:,1] = com_[1]
                mask_bt *= (np.sum( (xxx[:,0:2]-cmar)**2, axis=1) > THRESHOLD )

        btfd.write(np.packbits(mask_bt).tobytes())
        upfd.write(np.packbits(mask_up).tobytes())
        btfd.flush(); upfd.flush();

        print('.',end='')
        co += 1
        if co % 60 == 0:
            print('')
        try:
            if (not un.trajectory.next() or _ts.frame > maxframes):
                finished = True
        except StopIteration:
            finished = True
        sys.stdout.flush()
    # end while
# end double with
print("""
Pack size: %d
Frames converted: %d
Bye!
""" % (int(natm2/8), co) )

