#!/usr/bin/env python
"""
create *.atm from pdb file by:
getting charge from charmm psf file
assigning radii
"""
import sys
#
"""implement class and methods for handling pdb file"""
import sys
import string
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
        self.occ = []
        self.bfact = []
        self.atom_rads = {' C':1.8,' S':1.9,' O':1.6,' N':1.4,' P':1.8,' H':0.0,'ZN':1.4,' B':2.46,' F':1.4}
    def readfile(self,pdb_name):
        try:
            pdb_file = open(pdb_name,"r")
            #print 'opening file'
        except:
            # print "can't find file", pdb_name
            sys.exit("can't find file")
        contents = pdb_file.readlines()
        pdb_file.close()
        #print(type(contents))
        #print('lines read: ', len(contents))
        self.natom = 0
        xyz = [float(0.),float(0.),float(0.)]
        for entry in contents:
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
                self.natom +=1
                xyz[0] = float(entry[30:38])
                xyz[1] = float(entry[38:46])
                xyz[2] = float(entry[46:54])
                self.coords.append([xyz[0],xyz[1],xyz[2]])
                atname = entry[11:17]
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                self.occ.append(float(entry[54:60]))
                self.bfact.append(float(entry[60:67]))
                if(entry[12] != 'H'):
                    atype = entry[12:14]
                    #atype = ' ' + entry[12]
                else:
                    atype = ' ' + entry[12]
                #if(atype==' '):
                #    atype = entry[13]
                self.element.append(atype)
                if ( atype in self.atom_rads):
                    self.radius.append(self.atom_rads[atype])
                else:
                    print ("atom radius not in dictionary", atype)
                    self.radius.append(0.0)
#        print self.resnum[self.natom-1]

def pdb_write(pdb_file,pdb):
    for i in range(pdb.natom):
        string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f \n' % (i,pdb.name[i],pdb.res[i], \
                 pdb.chain[i], pdb.resnum[i], pdb.coords[i][0], \
                 pdb.coords[i][1],pdb.coords[i][2])
        pdb_file.write(string)

def res_write(pdb_file,pdb,ibeg,ifin):
    # note it inludes 'last' atom index ifin
    for i in range(ibeg,ifin+1):
        string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f \n' % (i,pdb.name[i],pdb.res[i], \
                 pdb.chain[i], pdb.resnum[i], pdb.coords[i][0], \
                 pdb.coords[i][1],pdb.coords[i][2])
        pdb_file.write(string)

#
print('create *.atm from pdb file by:')
print('getting charge from charmm psf file')
print('assigning radii')
if(len(sys.argv)<3):
  print('USAGE: psfpdb_to_atm.py psf_file pdb_file')
  sys.exit()
#
# read pdb file
#
pdb_file = sys.argv[2]
pdb = pdb_struct()
pdb.readfile(pdb_file)
print('pdb input # of atoms: ',pdb.natom)
#
# read psf file, extract charges
#
pdb_charge = []
psf_file = sys.argv[1]
pfile = open(psf_file,'r')
line = pfile.readline()
#while (line[10:15] != 'NATOM'):
while ('NATOM' not in line):
  line = pfile.readline()
natom = int(line[0:8])
print('psf input # of atoms: ',natom)
if(natom != pdb.natom):
  print('number of atoms in psf, pdb do not match:')
netcrg = 0.
for i in range(natom):
  line = pfile.readline()
  fields = line.split()
  atmcrg = float(fields[6])
  #atmcrg = float(line[34:44])
  pdb_charge.append(atmcrg)
  netcrg += atmcrg
print('# of atoms: ',natom,' netcrg: ',netcrg)
#
# open and write atm file
#
atm_file = pdb_file[0:-4]+'.atm'
print(atm_file)
pdb_out = open(atm_file,'w')
for i in range(natom):
  #print('q: %8.3f  r: %8.3f ' % (pdb_charge[i],pdb.radius[i]))
  string = 'ATOM %6d%6s%4s%1s%4s    %8.3f%8.3f%8.3f%6.2f%7.3f \n' % (i,pdb.name[i],pdb.res[i], \
  pdb.chain[i], pdb.resnum[i], pdb.coords[i][0], \
  pdb.coords[i][1],pdb.coords[i][2],pdb.radius[i],pdb_charge[i])
  pdb_out.write(string)
pdb_out.close()
