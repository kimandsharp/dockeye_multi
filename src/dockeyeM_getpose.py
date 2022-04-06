#!/usr/bin/env python
#############################################
#  Author: Kim Sharp
#  Date:   7/21/2018
# dockeyeM_getpose.py 
# extract pose from dockeye logfile
# add generation of movie file of ligand
# branch off dockeye_getpose.py to read models from input pdb file
# and best model out of logfile
# now applies inverse of protein rot/trans to ligand, so can use original protein
# file, avoid writing a new protein file sept. 2019
# and writes movie file first
#############################################
import sys
import math
from dockeye_methods import *
#=======================================
# MAIN
#=======================================
print('\n extract pose from dockeye logfile ')
print(' also generate movie file of ligand poses \n')
if(len(sys.argv) < 2):
  print('Usage: python dockeyeM_getpose.py dockeyeM_logfile')
  sys.exit()
file_name = sys.argv[1]
dockeye_log = open(file_name,'r')
contents = dockeye_log.readlines()
nline = len(contents)
dockeye_log.close()
print('lines read: ',nline)
#
# pdb files
#
iline = 0
fields = contents[iline].split()
pdbfile1 = fields[2]
iline += 1
fields = contents[iline].split()
pdbfile2 = fields[2]
print('pdb files used by dockeye: ')
print(pdbfile1,pdbfile2)
#
# get # of atoms
#
iline += 1
fields = contents[iline].split()
nat1 = int(fields[4])
nat2 = int(fields[6])
if(nat1 < nat2):
  print('warning: protein needs to be first pdb in file')
  sys.exit()
print('# atoms in protein, ligand: ',nat1,nat2)
#
# geometric centers
#
iline += 2
gcen1 = [0.,0.,0.]
fields = contents[iline].split()
for k in range(3):
  gcen1[k] = float(fields[k+1])
iline += 1
gcen2 = [0.,0.,0.]
fields = contents[iline].split()
for k in range(3):
  gcen2[k] = float(fields[k+1])
print('geometric centers: ',gcen1,gcen2)
header = True
iline += 1
while(header):
  head = contents[iline][0:7]
  if(head == 'new min'):
    header = False
    istart = iline
    #print('istart: ',istart)
  else:
    print(contents[iline][0:-1])
    iline += 1
#
#poses
#
npose = 0
et = []
ee = []
ev = []
nbest = []
while(iline < nline):
  head = contents[iline][0:7]
  #print(head)
  if(head == 'new min'):
    #print(head)
    fields = contents[iline].split()
    ee.append(float(fields[2]))
    ev.append(float(fields[3]))
    et.append(float(fields[4]))
    nbest.append(int(fields[6]))
    npose +=1
  iline += 1
print('# of poses: ',npose)
#
# generate movie file of ligand poses
#
print('generating movie of poses...')
pdbl = pdb_struct()
gcenp = [0.,0.,0.]
gcenl = [0.,0.,0.]
xyz = [0.,0.,0.]
pdbl.readligand(pdbfile2)
atm_file = pdbfile2[0:-4]+'_allposes.pdb'
for k in range(3):
  gcenp[k] = gcen1[k]
  gcenl[k] = gcen2[k]
pdb_out = open(atm_file,'w')
print('# of atoms: ',pdbl.natom)
iline = istart
nframe = 0
rmtp = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
trnp = [0.,0.,0.]
rmtl = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
trnl = [0.,0.,0.]
while(iline < nline):
  head = contents[iline][0:7]
  #print(head)
  if(head == 'new min'):
    fields = contents[iline].split()
    imod = int(fields[6])
    #print(head)
    nframe +=1
    #print ('f: ',nframe)
    pdb_out.write('MODEL %3d\n' % (nframe))
    #iline += 4
    #
    # extract rot, trans for this frame
    #
    for i in range(3):
      iline +=1
      fields = contents[iline].split()
      trnp[i] = float(fields[3])
      for j in range(3):
        rmtp[i][j] = float(fields[j])
    iline +=1
    for i in range(3):
      iline +=1
      fields = contents[iline].split()
      trnl[i] = float(fields[3])
      for j in range(3):
        rmtl[i][j] = float(fields[j])
    #print(rmtp)
    #print(trnp)
    #print(rmtl)
    #print(trnl)
    iline += 1
    #
    # find rotated origin in molc. local coords
    # then displacement of origin subtracted from apparant trans to get actual trans
    #
    gcenl_rot = rot_vec(rmtl,gcenl)
    gcenp_rot = rot_vec(rmtp,gcenp)
    for k in range(3):
      trnl[k] = trnl[k] - (gcenl[k] - gcenl_rot[k])
      trnp[k] = trnp[k] - (gcenp[k] - gcenp_rot[k])
    for i in range(pdbl.natom):
      #
      # apply rotations and translations
      # and inverse of protein rot/trans to ligand
      # in case user moved protein too- now ligand should be in
      # coord frame of original protein pdb
      i1 = i + imod*pdbl.natom
      for k in range(3):
        xyz[k] = pdbl.coords[i1][k] - gcenl[k]
      xyz1 = rot_vec(rmtl,xyz)
      for k in range(3):
        xyz[k] = xyz1[k] + gcenl[k] + trnl[k] - gcenp[k] - trnp[k]
      xyz2 = rot_vec(rmtp,xyz,inv=1)
      for k in range(3):
        xyz[k] = xyz2[k] + gcen1[k]
      string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdbl.name[i],pdbl.res[i], \
      pdbl.chain[i], pdbl.resnum[i], xyz[0],xyz[1],xyz[2], \
      pdbl.radius[i],pdbl.bfact[i])
      pdb_out.write(string)
    # end of pdb write
  # next frame
  pdb_out.write('ENDMDL\n')
  iline += 1
print('# of pose movie frames written: ',nframe)
print(' to file: ',atm_file)
print(' view poses in PyMol to help in selection')
pdb_out.close()
#
# choose pose
#
print(' ')
print('pose #      Ee           Ev           Et     model')
for i in range(npose):
  print('(%4d) %12.5f %12.5f %12.5f %6d' % (i+1,ee[i],ev[i],et[i],nbest[i]))
ipose = -1
while((ipose <0) or(ipose > npose)):
  print('select pose # (1 - ',npose,')',)
  print(' enter 0 to exit')
  ipose = int(input('>> '))
  #ipose = npose # debug
if(ipose > 0):
  indx = 9*(ipose - 1) + istart  + 1
  imod = nbest[ipose - 1]
  print('pose, model: ',ipose,imod)
  #
  # open pose pdb file
  if(ipose <10):
   numb = ('%1d'%(ipose))
  elif(ipose <100):
   numb = ('%2d'%(ipose))
  else:
   numb = ('%3d'%(ipose))
  atm_file = pdbfile2[0:-4]+'_pose_'+numb+'.pdb'
  print('writing pose to ',atm_file)
  pdb_out = open(atm_file,'w')
  #
  for i in range(3):
    fields = contents[indx].split()
    trnp[i] = float(fields[3])
    for j in range(3):
      rmtp[i][j] = float(fields[j])
    indx +=1
  indx +=1
  for i in range(3):
    fields = contents[indx].split()
    trnl[i] = float(fields[3])
    for j in range(3):
      rmtl[i][j] = float(fields[j])
    indx +=1
  print('\nprotein rotation matrix ')
  print(rmtp)
  print('\nprotein translation ')
  print(trnp)
  print('\nligand rotation matrix ')
  print(rmtl)
  print('\nligand translation ')
  print(trnl)
  # find rotated origin in molc. local coords
  # displacement of origin subtracted from apparant trans to get actual trans
  gcenl_rot = rot_vec(rmtl,gcenl)
  gcenp_rot = rot_vec(rmtp,gcenp)
  for k in range(3):
    trnl[k] = trnl[k] - (gcenl[k] - gcenl_rot[k])
    trnp[k] = trnp[k] - (gcenp[k] - gcenp_rot[k])
  pdb_out.write('MODEL %3d\n' % (ipose))
  for i in range(pdbl.natom):
    #
    # apply rotations and translations
    # and inverse of protein rot/trans to ligand
    # in case user moved protein too- now ligand should be in
    # coord frame of original protein pdb
    i1 = i + imod*pdbl.natom
    for k in range(3):
      xyz[k] = pdbl.coords[i1][k] - gcenl[k]
    xyz1 = rot_vec(rmtl,xyz)
    for k in range(3):
      xyz[k] = xyz1[k] + gcenl[k] + trnl[k] - gcenp[k] - trnp[k]
    xyz2 = rot_vec(rmtp,xyz,inv=1)
    for k in range(3):
      xyz[k] = xyz2[k] + gcen1[k]
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdbl.name[i],pdbl.res[i], \
    pdbl.chain[i], pdbl.resnum[i], xyz[0],xyz[1],xyz[2], \
    pdbl.radius[i],pdbl.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
  pdb_out.close()
