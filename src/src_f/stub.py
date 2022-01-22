#############################################
# code stub to check python-fortran inteface
# just read pdb files, send atom data and rotation matrices
#############################################
import sys
import math
from dockeye_methods import *
import dockeyeM_energy
import numpy  as np

#=======================================
#def de(pdbfile1="IL1B.atm",pdbfile2="MIM_tor.atm",charges=True,logscale=True,dielectric=80.,eps=0.1):
#def de(pdbfile1="ab.atm",pdbfile2="ag_tor.atm",charges=True,logscale=True,dielectric=80.,eps=0.1):
def de(pdbfile1="test1a.pdb",pdbfile2="test2a.pdb",charges=True,logscale=True,dielectric=80.,eps=0.0):
  # extract names, and create the 'object' name in the pymol menu window
  pdbobj1 = pdbfile1[:-4]
  #
  # Dockeye class reads pdb file upon init
  pdbobj2 = 'dockeye_lig'
  obj = Dockeye(pdbfile1,pdbfile2,pdbobj1,pdbobj2,charges,logscale,dielectric,eps)
  #
  return obj

#=======================================
class Dockeye():
  def __init__(self,pdbfile1,pdbfile2,pdbobj1,pdbobj2,charges,logscale,dielectric,eps):
    # calling arguments
    self.pdbfile1 = pdbfile1
    self.pdbfile2 = pdbfile2
    self.pdbobj1 = pdbobj1
    self.pdbobj2 = pdbobj2
    self.logscale = logscale
    self.dielectric = dielectric
    self.eps = eps
    print('Initializing Dockeye...')
    print('pdb file 1: ',pdbfile1)
    print('pdb file 2: ',pdbfile2)
    print('charges, logscale energy bars: ',charges,logscale)
    print('energy parameters dielectric: %8.3f VDW depth: %8.3f' % (dielectric,eps))
    # read original pdb coords
    self.pdb1 = pdb_struct()
    self.pdb1.readfile(self.pdbfile1)
    self.pdb2 = pdb_struct()
    self.pdb2.readligand(self.pdbfile2)
    self.iconf = self.pdb2.nmodel
    #self.iconf = 2
    #self.objmat1 = [1.,0.,0.,1., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    #self.objmat2 = [1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    pdbmat1 = [1.,0.,0.,1., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    pdbmat2 = [1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.]
    self.my_view = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    #print name
    #print self.pdb1.coords
    #print self.pdb2.coords
    qtot1 = 0.
    qtot2 = 0.
    self.gcen1 = [0.,0.,0.]
    self.energy = [0.,0.,0.]
    self.nbest = [0.0]
    for i in range(self.pdb1.natom):
      for k in range(3):
        self.gcen1[k] += self.pdb1.coords[i][k]
    for k in range(3):
      self.gcen1[k] /= self.pdb1.natom
    self.gcen2 = [0.,0.,0.]
    for i in range(self.pdb2.natom):
      for k in range(3):
        self.gcen2[k] += self.pdb2.coords[i][k]
    for k in range(3):
      self.gcen2[k] /= self.pdb2.natom
    if(not charges):
      for i in range(self.pdb1.natom):
        self.pdb1.bfact[i] = 0.
      for i in range(self.pdb2.natom):
        self.pdb2.bfact[i] = 0.
    if(charges):
      for i in range(self.pdb1.natom):
        qtot1 += self.pdb1.bfact[i]
      for i in range(self.pdb2.natom):
        qtot2 += self.pdb2.bfact[i]
    print('# of atoms 1: %6d   2: %6d' % (self.pdb1.natom,self.pdb2.natom))
    print('geometric centers: ')
    print('1:  %8.3f %8.3f %8.3f ' % (self.gcen1[0],self.gcen1[1],self.gcen1[2]))
    print('2:  %8.3f %8.3f %8.3f ' % (self.gcen2[0],self.gcen2[1],self.gcen2[2]))
    print('net charge 1:  %8.3f 2:  %8.3f ' % (qtot1,qtot2))
    #
    do_mm = True
    self.energy_min = 0.
    cgo_obj = pdb_interaction(pdbmat1,pdbmat2,self.pdb1,self.pdb2,self.gcen1,self.gcen2,
        self.energy,do_mm,self.logscale,self.dielectric,self.eps,self.nbest,self.energy_min)
#=======================================
#=======================================
def pdb_interaction(pdbmat1,pdbmat2,pdb1,pdb2,gcen1,gcen2,energy,do_mm,logscale,dielectric,eps,
                    nbest,emin):
  """
  # this is where we calculated interaction energy between
  # two molecules
  #
  # extract rotation matrices and translation vectors
  # note that if center of molecule not at origin
  # then rotation of a molecule moves the origin relative to its
  # local center, and that this displacment appears in translation so
  # this apparent trans needs to be removed to get actual translation
  """
  #=========================================
  # global parameters, these are the nominal values used by C
  # subroutine, energies are then post-scaled for actual parameters given as arguments to de()
  DIEL = 80.
  efact = 332./DIEL # dielectric constant factor gives kcal/mole for p+ unit charge, Angstroms
  EPS = 0.1 # depth of vdw potl. kcal/mole
  maxobjdata = 20000
  energy_obj = np.zeros(maxobjdata,float)
  #=========================================
  # extract rot mat
  rmt1 = [[pdbmat1[0],pdbmat1[1],pdbmat1[2]],
         [pdbmat1[4],pdbmat1[5],pdbmat1[6]],
         [pdbmat1[8],pdbmat1[9],pdbmat1[10]]]
  rmt2 = [[pdbmat2[0],pdbmat2[1],pdbmat2[2]],
         [pdbmat2[4],pdbmat2[5],pdbmat2[6]],
         [pdbmat2[8],pdbmat2[9],pdbmat2[10]]]
  #
  #
  xyz = [0.,0.,0.]
  #atom_data = [0.,0.] # stuff # atoms (as floats), coords, radii and charges in one long array to pass to energy_c
  atom_data = [0.,0.,0.] # stuff # atoms (as floats), models, coords, radii and 
  # charges in one long array to pass to energy_c
  nat1 = pdb1.natom
  nat2 = pdb2.natom
  nmod = pdb2.nmodel
  #print('data length 1: ',len(pdb1.coords))
  #print('data length 2: ',len(pdb2.coords))
  #
  # molecule 1
  #
  for i in range(nat1):
    atom_data.append(pdb1.radius[i])
  for i in range(nat1):
    atom_data.append(pdb1.bfact[i])
  for i in range(nat1):
    # apply rotations and translations
    for k in range(3):
      atom_data.append(pdb1.coords[i][k])
  #
  # molecule 2
  #
  for i in range(nat2):
    atom_data.append(pdb2.radius[i])
  for i in range(nat2):
    atom_data.append(pdb2.bfact[i])
  i = 0
  for nm in range(nmod):
    for j in range(nat2):
      # apply rotations and translations
      for k in range(3):
        atom_data.append(pdb2.coords[i][k])
      i += 1
  print('i: ',i)
  print('beg, end: ',pdb2.coords[0][0],pdb2.coords[i-1][2])
  #sys.exit()
  #
  #==========================================
  # C subroutine version
  #==========================================
  atom_data[0] = float(nat1) # now we know # of atoms, put in front of data array
  atom_data[1] = float(nat2)
  atom_data[2] = float(nmod)
  print('length of data: ',len(atom_data))
  #energy_obj = dockeyeM_energy.energy_f(energy_obj,atom_data)
  dockeyeM_energy.energy_f(energy_obj,atom_data)
  ndata = int(energy_obj[0])
  print('in calling program: ')
  il = int(energy_obj[0])
  for i in range(il):
    print(i,energy_obj[i])
  #print(energy_obj[0],energy_obj[1],energy_obj[2],energy_obj[3])
  #print(energy_obj[ndata-1],energy_obj[ndata-2],
  #          energy_obj[ndata-3],energy_obj[ndata-4],energy_obj[ndata-5])

  """
  ndata = int(energy_obj[0])
  #print('ndata: ',type(ndata))
  #print('ndata: ',ndata)
  # slice energy terms off end of data
  #energy[0] = energy_obj[ndata-3]
  nbest[0] = int(energy_obj[ndata-1])
  #print('from energy_c best model is: ',nbest[0])
  energy[1] = energy_obj[ndata-3]*DIEL/dielectric
  energy[2] = energy_obj[ndata-2]*eps/EPS
  energy[0] = energy[1] + energy[2]
  """
#
# now actually run 
de()
