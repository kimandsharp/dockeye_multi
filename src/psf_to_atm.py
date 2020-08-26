"""
getting charge from charmm psf file and assign radii
"""
import sys
import dockeye_methods as pdbm
#
print('get charges from charmm psf file')
print('assign radii and paste into occupancy, bfactor columns of pdb file')
if(len(sys.argv)<3):
  print('USAGE: psf_to_atm.py psf_file pdb_file')
  sys.exit()
#
# read pdb file
#
pdb_file = sys.argv[2]
pdb = pdbm.pdb_struct()
pdb.readfile(pdb_file)
print('pdb input # of atoms: ',pdb.natom)
#
# read psf file, extract charges
#
pdb_charge = []
psf_file = sys.argv[1]
pfile = open(psf_file,'r')
line = pfile.readline()
while (line[10:15] != 'NATOM'):
  line = pfile.readline()
natom = int(line[0:8])
print('psf input # of atoms: ',natom)
if(natom != pdb.natom):
  print('number of atoms in psf, pdb do not match:')
netcrg = 0.
for i in range(natom):
  line = pfile.readline()
  atmcrg = float(line[34:44])
  pdb_charge.append(atmcrg)
  netcrg += atmcrg
print('# of atoms: ',natom,' netcrg: ',netcrg)
#
# open and write atm file
#
atm_file = pdb_file[0:-4]+'_psf.pdb'
print('writing output to ',atm_file)
pdb_out = open(atm_file,'w')
for i in range(natom):
  #print('q: %8.3f  r: %8.3f ' % (pdb_charge[i],pdb.radius[i]))
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb.name[i],pdb.res[i], \
  pdb.chain[i], pdb.resnum[i], pdb.coords[i][0], \
  pdb.coords[i][1],pdb.coords[i][2],pdb.radius[i],pdb_charge[i])
  pdb_out.write(string)
pdb_out.close()
