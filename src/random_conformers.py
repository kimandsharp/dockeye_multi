"""
read pdb file with auto dock  torsion information from pdbqt
and sample conformers randomly, keep nconf with lowest energy
write out in MODEL/ENDMDL format
for pymol and dockeye multi
"""
import sys
import math
import matrix3 as mt
import numpy as np
import random as rn
import supq_flip as supq
rn.seed(7777777)

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
        ' Z':1.4,' B':2.46, '1':1.0, '2':1.0, '3':1.0, 'X':0.0}
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
  """ create bonded list based on overlap of vdw radii """
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
   
def conf_energy(pdb1,crds,non_bonds):
  """ energy of current config- use conformer crds not original
  pdb crds. assumes atm fomat input file (radii, charges in occ, bfact fields)
  """
  efact = 332./80. # dielectric of water
  vfact = 4.*0.1  # vdw depth of 0.1 kcal/mole
  dfact = 140 # approx value in kcal/proton charge^2/A^4 from Gilson/Honig (1991) JCAMD 5:5-20
  dxyz = [0.,0.,0.]
  evdw = 0.
  eelect = 0.
  edsolv = 0.
  for n in range(len(non_bonds)):
    i = non_bonds[n][0]
    j = non_bonds[n][1]
    sigma = pdb1.radius[i] + pdb1.radius[j]
    qi = pdb1.bfact[i]
    qj = pdb1.bfact[j]
    d2 = 0.
    for k in range(3):
      dxyz[k] = crds[i][k] - crds[j][k]
      d2 += dxyz[k]**2
    dd = math.sqrt(d2)
    rr = sigma/dd
    rr2 = rr*rr
    rr4 = rr2*rr2
    rr6 = rr2*rr2*rr2
    rr12 = rr6*rr6
    evdw += vfact*(rr12 - rr6) 
    eelect += efact*qi*qj/dd
    #edsolv += dfact*rr4*(qi**2 + qj**2)
    #print(i+1,j+1,sigma,dd,evdw,eelect)
  etot = eelect + evdw + edsolv
  #print('energies: ',evdw,eelect,edsolv,etot)
  return etot

def rmsd_conf(iconf_ref,iconf,conf_coords,natom):
  sd2 = 0.
  for i in range(natom):
    for k in range(3):
      sd2 += (conf_coords[iconf_ref][i][k] - conf_coords[iconf][i][k])**2
  rmsd = math.sqrt(sd2/natom)
  return rmsd
#
def sort_1_by_2(x,y,rev=False):
  """
  sort one list by elements in another list
  """
  #print('reverse',rev)
  if(len(x) == len(y)):
    y_x = zip(y,x)
    y_x_sorted = sorted(y_x,reverse=rev)
    y = [z[0] for z in y_x_sorted]
    x = [z[1] for z in y_x_sorted]
    return x,y
  else:
    print('lists of different length- not sorting')
# =======================
#         MAIN
# =======================
# input file
print('\nread pdb file with auto dock  torsion information from pdbqt')
print('and sample conformers randomly, option to generate 3 180o-flipmers, and ')
print('write out in MODEL/ENDMDL format for pymol and dockeye multi\n')
if(len(sys.argv) < 3):
  print('USAGE: python random_conformers.py pdbfile nconf {flipit (T/F)}')
  sys.exit()
pdbfile = sys.argv[1]
nconfmax = int(sys.argv[2])
print('# of conformers: ',nconfmax)
flipit = False
if(len(sys.argv) == 4):
  flipit = True
#
if(flipit):
  print('Generating 3 flipmers for each conformer')
else:
  print('Not generating flipmers')
print('reading input pdb file...')
pdb1 = pdb_struct()
pdb1.readfile(pdbfile)
qnet = 0.
for i in range(pdb1.natom):
  qnet += pdb1.bfact[i]
print('Net charge from bfactor column: %8.3f' %(qnet))
#
# output file
pdbout = pdbfile[0:-4] + '_ran.pdb'
print('writing conformers to file ',pdbout)
pdb_out = open(pdbout,'w')
#
# handle no torsion case
#
if(pdb1.ntor == 0):
  print('no torsions, just 1 conformer')
  nconfmax = 1
if(nconfmax == 1):
  nsample = 0
else:
  nsample = 200
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
#
conf_coords = np.zeros((nconfmax,pdb1.natom,3),'f')
conf_econf = np.zeros(nconfmax,'f')
conf_index = np.zeros(nconfmax,'int')
#
# make working copies of coords
#
xyz = [0.,0.,0.]
crds_old = []
crds_new = []
cen_old = [0.,0.,0.]
for i in range(pdb1.natom):
  for k in range(3):
    xyz[k] = pdb1.coords[i][k]
    cen_old[k] += xyz[k]
  crds_old.append([xyz[0],xyz[1],xyz[2]])
  crds_new.append([xyz[0],xyz[1],xyz[2]])
for k in range(3):
  cen_old[k] /= pdb1.natom
print('center of original conformer: ',cen_old)
#
# store original coords and its energy
#
nconf = 0
econf_old = conf_energy(pdb1,crds_old,non_bonds)
print('Initial energy: ',econf_old)
#conf_econf[nconf] = econf_old
#conf_index[nconf] = nconf
#for i in range(pdb1.natom):
#  for k in range(3):
#    conf_coords[nconf][i][k] = crds_old[i][k]
#nconf += 1
#
# init variables for sampling
pi = 3.14159
dangles = [-120.,120.,180.]
axis = [0.,0.,0.]
rmt1 = mt.Rotmat()
nrep = 0
for isample in range(nsample):
  # select random torsion, change
  # when list full, sort, replace worst with current if better
  itor = rn.randint(1,pdb1.ntor) - 1
  #print('selecting torsion: ',itor)
  #
  # get vector along rotating bond
  a1 = pdb1.branch_atoms[itor][0] - 1
  a2 = pdb1.branch_atoms[itor][1] - 1
  #print(itor,'twisting about ',a1,a2)
  for k in range(3):
    #print(crds[a2][k] , crds[a1][k])
    axis[k] = crds_new[a2][k] - crds_new[a1][k]
  mt.vnorm(axis)
  #print(axis)
  # debug
  #axis[0] = 0.
  #axis[1] = 0.
  #axis[2] = 1.
  #
  # generate rotation matrix for current rotation about this bond
  #chi = 180.
  iran = rn.randint(0,len(dangles)-1)
  chi = dangles[iran]
  #print('chi: ',chi)
  psi = 180.*math.asin(axis[2])/pi
  xy = math.sqrt(axis[0]**2 + axis[1]**2)
  if(xy < 1.e-6):
    phi = 90.
  else:
    phi = 180.*math.acos(axis[0]/xy)/pi
  if(axis[1] < 0.):phi = -1.*phi
  #print('phi,psi,chi: ',phi,psi,chi)
  rmt1.polar_rot(phi,psi,chi)
  #rmt1.printm()
  #
  # apply rotation to all atoms in this branch
  for j in range(a2+1,pdb1.branchend[itor]):
    #print('%4d ' % (j),end='')
    for k in range(3):
      #
      # shift atom so end of bond is origin
      xyz[k] = crds_new[j][k] - crds_new[a2][k]
    # apply rotation and shift back
    xyz1 = rmt1.rot_vec(xyz)
    #print('xyz:  ',xyz)
    #print('xyz1: ',xyz1)
    for k in range(3):
      crds_new[j][k] = xyz1[k] + crds_new[a2][k]
  #print(' \n')
  #
  # score current conformer
  econf_new = conf_energy(pdb1,crds_new,non_bonds)
  #print('new E: ',itor,econf_new,nconf)
  #
  # add to the list
  if(nconf < nconfmax):
    conf_econf[nconf] = econf_new
    conf_index[nconf] = nconf
    for i in range(pdb1.natom):
      for k in range(3):
        conf_coords[nconf][i][k] = crds_new[i][k]
    nconf += 1
  else:
    #print('reached max confs')
    # sort energies,
    conf_index_sort, conf_econf_sort = sort_1_by_2(conf_index,conf_econf,rev=False)
    iworst = conf_index_sort[nconf-1]
    eworst = conf_econf[iworst]
    #print('worst conf: ',iworst,eworst)
    if(econf_new < eworst): # replace worst conf with current one
      print('replacing ',iworst,econf_new,isample)
      nrep += 1
      conf_econf[iworst] = econf_new
      for i in range(pdb1.natom):
        for k in range(3):
          conf_coords[iworst][i][k] = crds_new[i][k]
    #sys.exit()
# end sample loop
print('number of times a better conf was found: ',nrep)
print(conf_econf,conf_index)
#print(conf_econf_sort,conf_index_sort)
# align all conformers before output
co1 = np.zeros((pdb1.natom,3),'f')
co2 = np.zeros((pdb1.natom,3),'f')
for i in range(pdb1.natom):
  for k in range(3):
    co1[i][k] = crds_old[i][k] - cen_old[k]
cen2 = [0.,0.,0.]
xyz = [0.,0.,0.]
xyz1 = [0.,0.,0.]
for nconf in range(nconfmax):
  #
  # copy and center coords of each conformaer
  for k in range(3):
    cen2[k] = 0.
  for i in range(pdb1.natom):
    for k in range(3):
      co2[i][k] = conf_coords[nconf][i][k]
      cen2[k] += co2[i][k]
  for k in range(3):
    cen2[k] /= pdb1.natom
  for i in range(pdb1.natom):
    for k in range(3):
      co2[i][k] -= cen2[k]
  #
  # rigid body align
  rms_min,qmin = supq.supq(co1,co2,pdb1.natom)
  rotmat,angle = supq.qttomt(qmin)
  #print('min rmsd, quat: ',rms_min,qmin,rotmat,angle)
  #
  # rotate and translate whole molecule
  for i in range(pdb1.natom):
    for k in range(3):
      xyz[k] = conf_coords[nconf][i][k] - cen2[k]
    for k in range(3):
      xyz1[k] = 0.
      for j in range(3):
        xyz1[k] += xyz[j]*rotmat[k][j]
      conf_coords[nconf][i][k] = xyz1[k] + cen_old[k]
# end alignment
# write conformers
pdb_out.write('MODEL %3d %12.5f\n' % (0,econf_old))
for i in range(pdb1.natom):
  for k in range(3):
    xyz[k] = crds_old[i][k]
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
  pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
  pdb1.radius[i],pdb1.bfact[i])
  pdb_out.write(string)
pdb_out.write('ENDMDL\n')
if(flipit):
  crd_rot = supq.crd_flip(crds_old,pdb1.natom)
  for nf in range(3):
    pdb_out.write('MODEL %3d %12.5f\n' % (0,econf_old))
    for i in range(pdb1.natom):
      for k in range(3):
        xyz[k] = crd_rot[i][k][nf]
      string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
      pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
      pdb1.radius[i],pdb1.bfact[i])
      pdb_out.write(string)
    pdb_out.write('ENDMDL\n')
#
for nconf in range(nconfmax):
  pdb_out.write('MODEL %3d %12.5f\n' % (nconf+1,conf_econf[nconf]))
  for i in range(pdb1.natom):
    for k in range(3):
      xyz[k] = conf_coords[nconf][i][k]
    string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
    pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
    pdb1.radius[i],pdb1.bfact[i])
    pdb_out.write(string)
  pdb_out.write('ENDMDL\n')
  if(flipit):
    crd_rot = supq.crd_flip(conf_coords[nconf],pdb1.natom)
    for nf in range(3):
      pdb_out.write('MODEL %3d %12.5f\n' % (nconf,conf_econf[nconf]))
      for i in range(pdb1.natom):
        for k in range(3):
          xyz[k] = crd_rot[i][k][nf]
        string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb1.name[i],pdb1.res[i], \
        pdb1.chain[i], pdb1.resnum[i], xyz[0],xyz[1],xyz[2], \
        pdb1.radius[i],pdb1.bfact[i])
        pdb_out.write(string)
      pdb_out.write('ENDMDL\n')
print('# of conformers: ',nconfmax)
pdb_out.close()
