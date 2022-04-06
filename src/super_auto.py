#!/usr/bin/env python
"""
python implementation of super_auto.f
c aligns atom in pdb file so as to minimize
c rms atom distance compared to reference pdb file
c branch off super.f that automatically seraches for matching atom pairs
21aug2020: added calculation of mean absolute dedivation (MAD)
july 2021 - for pdb inputs with duplicate atom names - e.g. from gold
prevent atoms being used > once for a match
"""
import sys
import math
import pdb_methods as pm
import numpy as np
#=====================================================
def supq(co1,co2,nat):
  """ find rotation matrix to minimix rmsd - we assume centers
  have already  been moved to the origin using method of Kearsley
  Acta Cryst (1989) A45:208 usng quaternions
  """
  #
  # calculate omega matrix of Kearsley et al.
  omega = np.zeros((4,4),'float')
  qmin = np.zeros(4,'float')
  for i in range(nat):
    # diagonal terms
    omega[0][0] += (co1[i][0] - co2[i][0])**2 + (co1[i][1] - co2[i][1])**2 + (co1[i][2] - co2[i][2])**2
    omega[1][1] += (co1[i][0] - co2[i][0])**2 + (co1[i][1] + co2[i][1])**2 + (co1[i][2] + co2[i][2])**2
    omega[2][2] += (co1[i][0] + co2[i][0])**2 + (co1[i][1] - co2[i][1])**2 + (co1[i][2] + co2[i][2])**2
    omega[3][3] += (co1[i][0] + co2[i][0])**2 + (co1[i][1] + co2[i][1])**2 + (co1[i][2] - co2[i][2])**2
    # off-diagonal terms
    omega[1][0] += (co1[i][1]+co2[i][1])*(co1[i][2]-co2[i][2]) - (co1[i][1]-co2[i][1])*(co1[i][2]+co2[i][2])
    omega[2][0] += (co1[i][2]+co2[i][2])*(co1[i][0]-co2[i][0]) - (co1[i][2]-co2[i][2])*(co1[i][0]+co2[i][0])
    omega[3][0] += (co1[i][0]+co2[i][0])*(co1[i][1]-co2[i][1]) - (co1[i][0]-co2[i][0])*(co1[i][1]+co2[i][1])
    omega[2][1] += (co1[i][0]-co2[i][0])*(co1[i][1]-co2[i][1]) - (co1[i][0]+co2[i][0])*(co1[i][1]+co2[i][1])
    omega[3][1] += (co1[i][0]-co2[i][0])*(co1[i][2]-co2[i][2]) - (co1[i][0]+co2[i][0])*(co1[i][2]+co2[i][2])
    omega[3][2] += (co1[i][1]-co2[i][1])*(co1[i][2]-co2[i][2]) - (co1[i][1]+co2[i][1])*(co1[i][2]+co2[i][2])
  #
  for i in range(4):
    for j in range(i+1,4):
      omega[i][j] = omega[j][i]
  #print('\nOmega matrix: \n',omega,'\n')
  eign,eigv = np.linalg.eig(omega)
  #print('\nEigen values: \n',eign)
  #print('\nEigen vectors: \n',eigv)
  imn = 0
  eign_mn = eign[0]
  imx = 0
  eign_mx = eign[0]
  for i in range(4):
    if(eign[i] < eign_mn):
      imn = i
      eign_mn = eign[i]
    if(eign[i] > eign_mx):
      imx = i
      eign_mx = eign[i]
  if(eign_mn < 1.e-6):
    print('Warning: eigen value close to or below 0: ',eign_mn)
  #print('min,max eign: ',eign[imn],eign[imx])
  rms_mn = math.sqrt(eign_mn/nat)
  rms_mx = math.sqrt(eign_mx/nat)
  #print('min, max rmsd possible by rotating: ',rms_mn,rms_mx)
  #print('Final RMS Deviation over selected set (A):: ',rms_mn)
  qmin[0] = - eigv[0][imn]
  for k in range(1,4):
    qmin[k] = eigv[k][imn]
  #print('min rmsd quaternion: ',qmin)
  return rms_mn,qmin

def qttomt(qt):
  """ convert quaternion to roation matrix"""
  rmt = np.zeros((3,3),'float')
  rmt[0][0] = qt[0]**2+qt[1]**2-qt[2]**2-qt[3]**2
  rmt[1][1] = qt[0]**2-qt[1]**2+qt[2]**2-qt[3]**2
  rmt[2][2] = qt[0]**2-qt[1]**2-qt[2]**2+qt[3]**2
  #
  rmt[1][0] = 2.*(qt[1]*qt[2]+qt[0]*qt[3])
  rmt[0][1] = 2.*(qt[1]*qt[2]-qt[0]*qt[3])
  #
  rmt[2][0] = 2.*(qt[1]*qt[3]-qt[0]*qt[2])
  rmt[0][2] = 2.*(qt[1]*qt[3]+qt[0]*qt[2])
  #
  rmt[2][1] = 2.*(qt[2]*qt[3]+qt[0]*qt[1])
  rmt[1][2] = 2.*(qt[2]*qt[3]-qt[0]*qt[1])
  trace = rmt[0][0] + rmt[1][1] + rmt[2][2]
  cosang = (trace - 1.)/2.
  angle = 180.*math.acos(cosang)/math.pi
  #print('angle: ',angle)
  return rmt,angle
#
#=====================================================
#main
#=====================================================
if(len(sys.argv) < 3):
  print('Usage: python super_auto.py ref_pdbfile mov_pdbfile')
  sys.exit()
#print(sys.argv)
#
# read pdb files
#
pdb_ref = pm.pdb_struct()
pdb_ref.readfile(sys.argv[1])
nref = pdb_ref.natom
print(sys.argv[1],' has ',nref,' atoms')
#
pdb_mov = pm.pdb_struct()
pdb_mov.readfile(sys.argv[2])
nmov = pdb_mov.natom
print(sys.argv[2],' has ',nmov,' atoms')
#
# write header for aligned file
#
pdb_out = open('super_auto_py.pdb','w')
head = 'REMARK  pdb file      ' + sys.argv[2] + '\n'
pdb_out.write(head)
head = 'REMARK  aligned with  ' + sys.argv[1] + '\n'
pdb_out.write(head)
#
# find matching atom pairs
#
nmatch = 0
indx1 = []
indx2 = []
used = []
for j in range(nmov):
  used.append(False)
for i in range(nref):
  for j in range(nmov):
    if(used[j]): continue # prevent using atom twice when atom names not unique
    if(i not in indx2):
      ifind = 1
      if(pdb_ref.name[i]   != pdb_mov.name[j]):   ifind = 0
      if(pdb_ref.res[i]    != pdb_mov.res[j]):    ifind = 0
      if(pdb_ref.resnum[i] != pdb_mov.resnum[j]): ifind = 0
      if(ifind == 1):
        nmatch += 1
        used[j] = True # prevent using atom twice when atom names not unique
        #print('matched ref ',i,' moving ',j)
        indx1.append(j)
        indx2.append(i)
        break
print('# of matched atom pairs: ',nmatch)
pdb_out.write('REMARK # of matched atom pairs: %6d\n'%(nmatch))
print('----')
#
# make working copy of coords
#
crd_ref = np.zeros((nmatch,3),'float')
crd_mov = np.zeros((nmatch,3),'float')
xyz = [0.,0.,0.]
#rmsd = 0.
for i in range(nmatch):
  iref = indx2[i]
  imov = indx1[i]
  print(imov,pdb_mov.name[imov],iref,pdb_ref.name[iref])
  for k in range(3):
    crd_ref[i][k]= pdb_ref.coords[iref][k]
    crd_mov[i][k]= pdb_mov.coords[imov][k]
    #del2 = (crd_mov[i][k] - crd_ref[i][k])**2
    #rmsd += del2
  #print(crd_ref[i],crd_mov[i],del2)
print('----')
#rmsd = math.sqrt(rmsd/nmatch)
#print('rmsd: ',rmsd) # check
#print(crd_ref)
#print(crd_mov)
#
# initial rmsdev, and centroids
#
rms = 0.
mad = 0.
cen_ref = [0.,0.,0.]
cen_mov = [0.,0.,0.]
trn_vec = [0.,0.,0.]
for i in range(nmatch):
  dist2 = 0.
  for k in range(3):
    cen_ref[k] += crd_ref[i][k]
    cen_mov[k] += crd_mov[i][k]
    dist2 += (crd_ref[i][k] - crd_mov[i][k])**2
    #rms += (crd_ref[i][k] - crd_mov[i][k])**2
  rms += dist2
  mad += math.sqrt(dist2)
rms = math.sqrt(rms/nmatch)
mad /= nmatch
print('Initial Rms Deviation, MAD over selected set (A):: %8.3f %8.3f'%(rms,mad))
pdb_out.write('REMARK Initial Rms Deviation, MAD over selected set (A):: %8.3f %8.3f\n'%(rms,mad))
for k in range(3):
  cen_ref[k] = cen_ref[k]/nmatch
  cen_mov[k] = cen_mov[k]/nmatch
  trn_vec[k] = cen_ref[k] - cen_mov[k]
#print('ref molc centroid: ',cen_ref)
#print('mov molc centroid: ',cen_mov)
#print('translation vector: ',trn_vec)
pdb_out.write('REMARK translation vector: %8.3f %8.3f %8.3f\n'%(trn_vec[0],trn_vec[1],trn_vec[2]))
trn_mag = math.sqrt(trn_vec[0]**2 + trn_vec[1]**2 + trn_vec[2]**2)
#
#  move molecules to the origin before finding rotation
#
rms = 0.
for i in range(nmatch):
  for k in range(3):
    crd_ref[i][k] -= cen_ref[k]
    crd_mov[i][k] -= cen_mov[k]
    rms += (crd_ref[i][k] - crd_mov[i][k])**2
rms = math.sqrt(rms/nmatch)
print('Rms Deviation after translation (A):: %8.3f'%(rms))
#
# find best rotation matrix
rms_mn,qmin = supq(crd_ref,crd_mov,nmatch)
#print('min rmsd quaternion: ',qmin)
#
# convert quaternion to rotation matrix and rotation magnitude
rotmat,angle = qttomt(qmin)
#print(rotmat)
print('Magnitude of Translation (A) & rotation (o): %8.2g %8.2f\n' % (trn_mag,angle))
pdb_out.write('REMARK Magnitude of Translation (A) & rotation (o): %8.2g %8.2f\n' % (trn_mag,angle))
#print('Final RMS Deviation over selected set (A):: ',rms_mn)
#pdb_out.write('REMARK Final Rms Deviation over selected set (A):: %8.3f\n'%(rms_mn))
#
# apply rotation to moving set so we can compute final MAD 
# (and final rms again, but directly, not by minimum eigen value)
#
xyz1 = [0.,0.,0.]
rms_mn = 0.
mad_mn = 0.
for n in range(nmatch):
  for j in range(3):
    xyz1[j] = 0.
    for k in range(3):
      xyz1[j] += rotmat[j][k]*crd_mov[n][k]
  dist2 = 0.
  for k in range(3):
    dist2 += (crd_ref[n][k] - xyz1[k])**2
  rms_mn += dist2
  mad_mn += math.sqrt(dist2)
rms_mn = math.sqrt(rms_mn/nmatch)
mad_mn /= nmatch
print('Final Rms Deviation, MAD over selected set (A):: %8.3f %8.3f'%(rms_mn,mad_mn))
pdb_out.write('REMARK Final Rms Deviation, MAD over selected set (A):: %8.3f %8.3f\n'%(rms_mn,mad_mn))
#
# apply rotation and translation to entire moving set and write out new pdb file
#
for n in range(nmov):
  for k in range(3):
    xyz[k] = pdb_mov.coords[n][k] - cen_mov[k]
  for j in range(3):
    xyz1[j] = 0.
    for k in range(3):
      xyz1[j] += rotmat[j][k]*xyz[k]
  for k in range(3):
    pdb_mov.coords[n][k] =xyz1[k] +  cen_ref[k]
pm.pdb_write(pdb_out,pdb_mov)
pdb_out.close()
