"""
read pdb file with auto dock  torsion information from pdbqt
and sample conformers, write out in MODEL/ENDMDL format
for pymol and dockeye multi
"""
import sys
import math
import matrix3 as mt
import numpy as np

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
        self.atom_rads = {' C':1.8,' S':1.9,' O':1.6,' N':1.4,' P':1.8,' H':1.0,'ZN':1.4,
        ' Z':1.4,' B':2.46, '1':1.0, '2':1.0, '3':1.0}
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
                print('ERROR: found ',nbranch,self.ntor)
                sys.exit()
            if(entry[0:6] == 'REMARK'):
              if('A    between atoms' in entry):
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
def make_bonds(pdb1):
  """ create bonded list based on overalp ov vdw radii """
  margin = 1.3
  bond_list = []
  nbond = 0
  for i in range(pdb1.natom):
    for j in range(i+1,pdb1.natom):
      dist2 = 0.
      for k in range(3):
        dist2 += (pdb1.coords[i][k] - pdb1.coords[j][k])**2
      dist = math.sqrt(dist2)
      overlap = dist + margin - pdb1.radius[i] - pdb1.radius[j]
      if (overlap < 0.):
        nbond +=1
        bond_list.append([i,j])
  #print('# of bonds: ',nbond)
  return bond_list
  
def make_nonbonds(bond_list,natom):
  exclusion_list = np.zeros((natom,natom),'I')
  bond_order = np.zeros((natom),'I')
  atom_bonds = []
  # first make a list of bonded atoms for each atom
  for n in range(natom):
    atom_bonds.append([])
  for n in range(len(bond_list)):
    i = bond_list[n][0]
    j = bond_list[n][1]
    bond_order[i] += 1
    bond_order[j] += 1
    atom_bonds[i].append(j)
    atom_bonds[j].append(i)
  #print(bond_order)
  #print(atom_bonds)
  for n in range(len(bond_list)):
    i = bond_list[n][0]
    j = bond_list[n][1]
    #print('bond ',i+1,j+1)
    exclusion_list[i][j] = 1 # add 1-2 pairs 
    exclusion_list[j][i] = 1 # add 1-2 pairs 
    for i1 in range(bond_order[i]):
      j1 = atom_bonds[i][i1]
      #print('atom 3i ',j+1,j1+1)
      exclusion_list[j][j1] = 1 # add 1-3 pairs 
      exclusion_list[j1][j] = 1 # add 1-3 pairs 
    for j1 in range(bond_order[j]):
      i1 = atom_bonds[j][j1]
      #print('atom 3j ',i+1,i1+1)
      exclusion_list[i1][i] = 1 # add 1-3 pairs 
      exclusion_list[i][i1] = 1 # add 1-3 pairs 
    for i1 in range(bond_order[i]):
      j2 = atom_bonds[i][i1]
      for j1 in range(bond_order[j]):
        i2 = atom_bonds[j][j1]
        #print('atoms 1-4',i2+1,j2+1)
        exclusion_list[i2][j2] = 1 # add 1-4 pairs 
        exclusion_list[j2][i2] = 1 # add 1-4 pairs 
  # finally create list of nonbond interactions- whats  not 1-2,1-3,14
  non_bonds = []
  for i in range(natom):
    for j in range(i+1,natom):
      if(exclusion_list[i][j] == 0):
        non_bonds.append([i,j])
  return(non_bonds)
   

# =======================
#         MAIN
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
pdbout = pdbfile[0:-4] + '_cnf.pdb'
print('writing conformers to file ',pdbout)
pdb_out = open(pdbout,'w')
#
# handle no torsion case
#
if(pdb1.ntor == 0):
  print('no torsions, just 1 conformer')
  nconf = 1
  pdb_out.write('MODEL %3d\n' % (nconf))
  for i in range(pdb1.natom):
    #
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], pdb1.coords[i][0],pdb1.coords[i][1],pdb1.coords[i][2], \
    pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
  sys.exit()
#
# create list of bonds
#
bond_list = make_bonds(pdb1)
print('# of bonds: ',len(bond_list))
#print(bond_list)
#
# create list of nonbonds
#
non_bonds = make_nonbonds(bond_list,pdb1.natom)
print('# of nonbond interactions: ',len(non_bonds))
#print(non_bonds)
sys.exit()

pi = 3.14159
xyz = [0.,0.,0.]
axis = [0.,0.,0.]
crds = []
for i in range(pdb1.natom):
  for k in range(3):
    crds.append([xyz[0],xyz[1],xyz[2]])
rmt1 = mt.Rotmat()
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
