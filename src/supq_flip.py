import numpy as np
from numpy import linalg as la
import math
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
  if(eign_mn < 1.e-3):
    print('Warning: eigen value close to or below 0: ',eign_mn)
    eign_mn = 0.
  #print('min,max eign: ',eign[imn],eign[imx])
  rms_mn = math.sqrt(eign_mn/nat)
  rms_mx = math.sqrt(eign_mx/nat)
  print('min, max rmsd possible by rotating: ',rms_mn,rms_mx)
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
  print('angle: ',angle)
  return rmt,angle

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
