"""
read an auto dock format pdb file (*.pdbqt)
and generate conformers, write out in MODEL/ENDMDL format
for pymol and dockeye multi
14 oct 2019, add 180 degree flips around axes of MOI
"""
import sys
import math
import matrix3 as mt
import numpy as np
from numpy import linalg as la

#
# include pdb_methods so code selfcontained
# and we need special code for parsing the *.pdbqt ROOT/BRANCH
# torsion tree records
"""implement specific class and methods for handling pdb with torsion
records file"""
#import string
class pdb_struct:
    def __init__(self):
        """create some placeholders for coords, radius, atom name, atom type etc"""
        self.coords = []
        self.radius = []
        self.name = []
        self.res = []
        self.element = []
        self.resnum = []
        self.chain = []
        self.bfact = []
        self.atom_rads = {' C':1.8,' S':1.9,' O':1.6,' N':1.4,' P':1.8,' H':0.0,'ZN':1.4,
        ' Z':1.4,' B':2.46}
        self.root = [1,0]
    def readfile(self,pdb_name):
        try:
            pdb_file = open(pdb_name,"r") # python3
            #pdb_file = file(pdb_name,"r")  # python2
            #print 'opening file'
        except:
            # print "can't find file", pdb_name
            sys.exit("can't find file")
        contents = pdb_file.readlines()
        pdb_file.close()
        self.natom = 0
        xyz = [0.,0.,0.]
        branch_stack = []
        nbranch = 0
        self.branchend = []
        self.branch_atoms = []
        self.ntor = 0
        for entry in contents:
            if(entry[0:7] == 'TORSDOF'):
              fields = entry.split()
              ntor_expt = int(fields[1])
              print ('expected # of torsions: ',ntor_expt)
              if((ntor_expt != nbranch) or (ntor_expt != self.ntor)):
                print('ERROR: found ',nbranch)
                sys.exit()
            if(entry[0:6] == 'REMARK'):
              if('between atoms' in entry):
                self.ntor += 1
            if(entry[0:7] == 'ENDROOT'):
              self.root[1] = self.natom
              print('found root: ',self.root)
            if(entry[0:6] == 'BRANCH'):
              nbranch += 1
              fields = entry.split()
              m1 = int(fields[1])
              m2 = int(fields[2])
              self.branch_atoms.append([m1,m2])
              branch_stack.insert(0,nbranch)
              self.branchend.append(0)
            if(entry[0:9] == 'ENDBRANCH'):
              i = branch_stack.pop(0) - 1
              self.branchend[i] = self.natom
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
                self.natom +=1
                xyz[0] = float(entry[30:38])
                xyz[1] = float(entry[38:46])
                xyz[2] = float(entry[46:54])
                self.coords.append([xyz[0],xyz[1],xyz[2]])
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17]
                #atom_name_number[atname] = self.natom
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                    atype = entry[12:14]
                    #atype = ' ' + entry[12]
                else:
                    atype = ' ' + entry[12]
                #if ( self.atom_rads.has_key(atype)):
                if ( atype in self.atom_rads):
                    self.radius.append(self.atom_rads[atype])
                else:
                    print ("atom radius not in dictionary", atype)
                    self.atom_rads[atype] = 0.0
                    self.radius.append(0.0)
        print('# of branches: ',nbranch)
        if(len(branch_stack) !=0):
          print('ERROR: branches: ',branch_stack)
        print('branchends: ',self.branchend)
        print('branch atoms: ',self.branch_atoms)
        """
        for i in range(self.ntor):
          a1 = self.branch_atoms[i][0]
          a2 = self.branch_atoms[i][1]
          print(i,'twisting about ',a1,a2)
          for j in range(a2,self.branchend[i]):
            print('%4d ' % (j+1),end='')
          print(' \n')
        """
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
# input file
if(len(sys.argv) < 2):
  pdbfile = input('pdb file containing adt torsion records>> ')
else:
  pdbfile = sys.argv[1]
pdb1 = pdb_struct()
pdb1.readfile(pdbfile)
#
# output file
pdbout = pdbfile[0:-4] + '_torf.pdb'
print('writing conformers to file ',pdbout)
pdb_out = open(pdbout,'w')
#
# setup
pi = 3.14159
axis = [0.,0.,0.]
rmt1 = mt.Rotmat()
#
# space for working copy of coords
xyz = [0.,0.,0.]
crds = []
for i in range(pdb1.natom):
  for k in range(3):
    crds.append([xyz[0],xyz[1],xyz[2]])
#
# handle no torsion case
if(pdb1.ntor == 0):
  print('no torsions, just 1 conformer')
  #
  # make working copy of coords
  for i in range(pdb1.natom):
    for k in range(3):
      crds[i][k] = pdb1.coords[i][k]
  nconf = 1
  nflip = 1
  #
  # write conformer
  pdb_out.write('MODEL %3d   %3d\n' % (nconf,nflip))
  for i in range(pdb1.natom):
    #
    for k in range(3):
      xyz[k] = crds[i][k]
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
    pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
  #
  # generate and write the 3 flipped conformers
  crds1 = crd_flip(crds,pdb1.natom)
  for l in range(3):
    nflip = l + 1
    nconf += 1
    pdb_out.write('MODEL %3d   %3d\n' % (nconf,nflip))
    for i in range(pdb1.natom):
      #
      for k in range(3):
        xyz[k] = crds1[i][k][l]
      string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
      pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
      pdb1.radius[i],pdb1.bfact[i])
      pdb_out.write(string)
    pdb_out.write('ENDMDL\n')
  sys.exit()
#
# generate rotamers using odometer method
#
rotamers = []
#
#pdb1.ntor = 2 # debug
nconf = 0
for i in range(pdb1.ntor):
  rotamers.append(0)
ic = 0
while(ic < pdb1.ntor):
  nconf += 1
  print(rotamers)
  #===============
  #
  # make working copy of coords
  #
  for i in range(pdb1.natom):
    for k in range(3):
      crds[i][k] = pdb1.coords[i][k]
  #
  # trace torsion tree applying rotations
  #
  for i in range(pdb1.ntor):
    #
    # get vector along rotating bond
    a1 = pdb1.branch_atoms[i][0] - 1
    a2 = pdb1.branch_atoms[i][1] - 1
    #print(i,'twisting about ',a1,a2)
    for k in range(3):
      #print(crds[a2][k] , crds[a1][k])
      axis[k] = crds[a2][k] - crds[a1][k]
    mt.vnorm(axis)
    #print(axis)
    # debug
    #axis[0] = 0.
    #axis[1] = 0.
    #axis[2] = 1.
    #
    # generate rotation matrix for current rotation about this bond
    #chi = 180.
    chi = rotamers[i]*120.
    psi = 180.*math.asin(axis[2])/pi
    xy = math.sqrt(axis[0]**2 + axis[1]**2)
    if(xy < 1.e-6):
      phi = 90.
    else:
      phi = 180.*math.acos(axis[0]/xy)/pi
    if(axis[1] < 0.):phi = -1.*phi
    print('phi,psi,chi: ',phi,psi,chi)
    rmt1.polar_rot(phi,psi,chi)
    rmt1.printm()
    #
    # apply rotation to all atoms in this branch
    for j in range(a2+1,pdb1.branchend[i]):
      #print('%4d ' % (j),end='')
      for k in range(3):
        #
        # shift atom so end of bond is origin
        xyz[k] = crds[j][k] - crds[a2][k]
      # apply rotation and shift back
      xyz1 = rmt1.rot_vec(xyz)
      #print('xyz:  ',xyz)
      #print('xyz1: ',xyz1)
      for k in range(3):
        crds[j][k] = xyz1[k] + crds[a2][k]
    #print(' \n')
  #
  # write current conformer
  pdb_out.write('MODEL %3d\n' % (nconf))
  for i in range(pdb1.natom):
    #
    for k in range(3):
      xyz[k] = crds[i][k]
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
    pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
  #
  # generate and write the 3 flipped conformers
  crds1 = crd_flip(crds,pdb1.natom)
  for l in range(3):
    nflip = l + 1
    nconf += 1
    pdb_out.write('MODEL %3d   %3d\n' % (nconf,nflip))
    for i in range(pdb1.natom):
      #
      for k in range(3):
        xyz[k] = crds1[i][k][l]
      string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
      pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
      pdb1.radius[i],pdb1.bfact[i])
      pdb_out.write(string)
    pdb_out.write('ENDMDL\n')
  #===============
  icarry = 1
  while(icarry !=0):
    rotamers[ic] += 1
    if(rotamers[ic] < 3):
      icarry = 0
    else:
      rotamers[ic] = 0
      ic += 1
      if(ic == pdb1.ntor): break
  if(ic == pdb1.ntor): break
  ic = 0
print('# of conformers: ',nconf)
pdb_out.close()
