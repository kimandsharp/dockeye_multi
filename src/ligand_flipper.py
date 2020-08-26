"""
subset st of make_conformers_flip.py that just
generates the three poses created by
180 degree flips around axes of MOI
esp useful for ligands with exact/almost 2-fold symmetry
"""
import sys
import math
import matrix3 as mt
import numpy as np
from numpy import linalg as la
import pdb_methods

def crd_flip(crds,natom):
   """ for any ligand conformer, generate 3 flipped versions around x,y,z 
       principal moment of inertia axes """
   #
   # find center, and set to 0,0,0
   gcen = np.zeros(3)
   for i in range(natom):
     for k in range(3):
       gcen[k] += crds[i][k]
   for k in range(3):
     gcen[k] /= natom
   #print('center: ',gcen)
   for i in range(natom):
     for k in range(3):
       crds[i][k] -= gcen[k]
   #
   # find moi
   moi = np.zeros((3,3))
   for i in range(natom):
     for k in range(3):
       for j in range(3):
         moi[k][j] += crds[i][k]*crds[i][j]
   #print(moi)
   #
   # find principal axes (eigenvectors of moi)
   # which will form rotation matrix to rotate into xyz frame
   # then flip around x, y, and then z
   e,ev = la.eigh(moi)
   #print(e)
   #print(ev)
   crds_rot = np.zeros((natom,3,3))
   xyz_rot = np.zeros(3)
   xyz_flip = np.zeros((3,3))
   for i in range(natom):
     for k in range(3):
       #
       # rotate original with major axes aligned along xyz
       xyz_rot[k] = 0.
       for j in range(3):
         xyz_rot[k] += crds[i][j]*ev[j][k]
       #
       # flip around x, y, and z
     # x rotation
     xyz_flip[0][0] = +xyz_rot[0]
     xyz_flip[1][0] = -xyz_rot[1]
     xyz_flip[2][0] = -xyz_rot[2]
     # y rotation
     xyz_flip[0][1] = -xyz_rot[0]
     xyz_flip[1][1] = +xyz_rot[1]
     xyz_flip[2][1] = -xyz_rot[2]
     # z rotation
     xyz_flip[0][2] = -xyz_rot[0]
     xyz_flip[1][2] = -xyz_rot[1]
     xyz_flip[2][2] = +xyz_rot[2]
     for l in range(3):
       #
       # rotate and translate back
       for k in range(3):
         crds_rot[i][k][l] = gcen[k]
         for j in range(3):
           crds_rot[i][k][l] += xyz_flip[j][l]*ev[k][j]
   return crds_rot
# =======================
#
print('\ngenerates the three poses created by')
print('180 degree flips around axes of MOI of a ligand')
print('esp useful for ligands with exact/almost 2-fold symmetry\n')
#
# input file
if(len(sys.argv) < 2):
  pdbfile = input('pdb file containing you want to flip>> ')
else:
  pdbfile = sys.argv[1]
pdb1 = pdb_methods.pdb_struct()
pdb1.readfile(pdbfile)
#
# output file
fileout = pdbfile[:-4] + '_flip.pdb'
pdb_out = open(fileout,'w')
#
# write out & make working copy of original coords
xyz = [0.,0.,0.]
crds = []
nflip = 1
pdb_out.write('MODEL %3d \n' % (nflip))
for i in range(pdb1.natom):
  for k in range(3):
    xyz[k] = pdb1.coords[i][k]
  crds.append([xyz[0],xyz[1],xyz[2]])
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
  pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
  pdb1.radius[i],pdb1.bfact[i])
  pdb_out.write(string)
pdb_out.write('ENDMDL\n')
#
# generate and write the 3 flipped conformers
crds1 = crd_flip(crds,pdb1.natom)
for l in range(3):
  nflip = l + 2
  pdb_out.write('MODEL %3d \n' % (nflip))
  for i in range(pdb1.natom):
    #
    for k in range(3):
      xyz[k] = crds1[i][k][l]
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
    pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
pdb_out.close()
sys.exit()
