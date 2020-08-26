#############################################
#  Author: Kim Sharp
#  Date:   7/21/2017
#  
#############################################
import sys
import math
#from pdb_methods import *
from dockeye_methods import *
print('=================')
print('dockeye_prep: position two molecules nicely for pymol/dockeye ')
print('=================')
if(len(sys.argv) <3):
  print ('usage: dockeye_prep.py protein_pdbfile ligand_pdbfile')
  sys.exit()
pdbfile1 = sys.argv[1]
pdbfile2 = sys.argv[2]
pdb1 = pdb_struct()
pdb1.readfile(pdbfile1)
pdb2 = pdb_struct()
pdb2.readfile(pdbfile2)
#
# find residues near ligand
#
rcut = 5.
rcut2 = rcut**2
resnam = '****' # dummy name to start off
resnum = '****'
reschn = '*'
for i in range(pdb1.natom):
  if((pdb1.res[i] != resnam)or(pdb1.resnum[i] != resnum) or(pdb1.chain[i] != reschn)):
    for j in range(pdb2.natom):
      d2 = 0.
      for k in range(3):
        d2 += (pdb1.coords[i][k] - pdb2.coords[j][k])**2
      if(d2 < rcut2):
        resnam = pdb1.res[i]
        resnum = pdb1.resnum[i]
        reschn = pdb1.chain[i]
        print(resnam,resnum,reschn,' < ',rcut, ' from ligand')
        break
#
# position ligand standard distance in -X direction wrt proteins
gcen1 = [0.,0.,0.]
gmax1 = [pdb1.coords[0][0],pdb1.coords[0][1],pdb1.coords[0][2]]
gmin1 = [pdb1.coords[0][0],pdb1.coords[0][1],pdb1.coords[0][2]]
for i in range(pdb1.natom):
  for k in range(3):
    gcen1[k] += pdb1.coords[i][k]
    gmax1[k] = max(gmax1[k],pdb1.coords[i][k])
    gmin1[k] = min(gmin1[k],pdb1.coords[i][k])
for k in range(3):
  gcen1[k] /= pdb1.natom
#
gcen2 = [0.,0.,0.]
gmax2 = [pdb2.coords[0][0],pdb2.coords[0][1],pdb2.coords[0][2]]
gmin2 = [pdb2.coords[0][0],pdb2.coords[0][1],pdb2.coords[0][2]]
for i in range(pdb2.natom):
  for k in range(3):
    gcen2[k] += pdb2.coords[i][k]
    gmax2[k] = max(gmax2[k],pdb2.coords[i][k])
    gmin2[k] = min(gmin2[k],pdb2.coords[i][k])
for k in range(3):
  gcen2[k] /= pdb2.natom
print( '# of atoms 1: %6d   2: %6d' % (pdb1.natom,pdb2.natom))
print( 'geometric centers: ')
print( '1:  %8.3f %8.3f %8.3f ' % (gcen1[0],gcen1[1],gcen1[2]))
print( '2:  %8.3f %8.3f %8.3f ' % (gcen2[0],gcen2[1],gcen2[2]))
print( 'min max coords: ')
print( '1 min:  %8.3f %8.3f %8.3f ' % (gmin1[0],gmin1[1],gmin1[2]))
print( '1 max:  %8.3f %8.3f %8.3f ' % (gmax1[0],gmax1[1],gmax1[2]))
print( '2 min:  %8.3f %8.3f %8.3f ' % (gmin2[0],gmin2[1],gmin2[2]))
print( '2 max:  %8.3f %8.3f %8.3f ' % (gmax2[0],gmax2[1],gmax2[2]))
xsep = [0.,0.,0.]
xsep[0] = ((gmax2[0] - gmin2[0]) + (gmax1[0] - gmin1[0]))
print('offsetting molecule2 by %8.3f %8.3f %8.3f ' % (xsep[0],xsep[1],xsep[2]))
print('and rotating by Z-120o')
cs = -0.5
sn = math.sqrt(0.75)
rmt = [[cs,sn,0.],[-sn,cs,0.],[0.,0.,1.]]
xyz = [0.,0.,0.]
for i in range(pdb2.natom):
  #for k in range(3):
  #  pdb2.coords[i][k] = pdb2.coords[i][k] - xsep[k] - gcen2[k] + gcen1[k] # translate
  for k in range(3):
    xyz[k] = pdb2.coords[i][k] - gcen2[k]
  for k in range(3):
    pdb2.coords[i][k] = - xsep[k] + gcen1[k]
    for j in range(3):
      pdb2.coords[i][k] += xyz[j]*rmt[j][k]
qtot1 = 0.
qtot2 = 0.
for i in range(pdb1.natom):
  qtot1 += pdb1.bfact[i]
for i in range(pdb2.natom):
  qtot2 += pdb2.bfact[i]
print( 'net charge 1:  %8.3f 2:  %8.3f ' % (qtot1,qtot2))
#
# open and write atm files
#
pdb1_out = open('de_prot.pdb','w')
for i in range(pdb1.natom):
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
  pdb1.chain[i], pdb1.resnum[i], pdb1.coords[i][0], \
  pdb1.coords[i][1],pdb1.coords[i][2],pdb1.occ[i],pdb1.bfact[i])
  pdb1_out.write(string)
pdb1_out.close()
#
pdb2_out = open('de_ligand.pdb','w')
for i in range(pdb2.natom):
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb2.name[i],pdb2.res[i], \
  pdb2.chain[i], pdb2.resnum[i], pdb2.coords[i][0], \
  pdb2.coords[i][1],pdb2.coords[i][2],pdb2.occ[i],pdb2.bfact[i])
  pdb2_out.write(string)
pdb2_out.close()

