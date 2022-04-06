#=======================================
# include pdb_methods so code selfcontained
"""implement class and methods for handling pdb file"""
import sys
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
        self.occ = []
        self.bfact = []
        self.atom_rads = {' C':1.8,' S':1.9,' O':1.6,' N':1.4,' P':1.8,' H':0.0,'ZN':1.4,
        ' Z':1.4,' B':2.46,'MG':2.00,' F':1.56, '1H':0.0, '2H':0.0, '3H':0.0,}
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
                self.occ.append(float(entry[54:60]))
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17]
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                    atype = entry[12:14]
                    #atype = ' ' + entry[12]
                else:
                    atype = ' ' + entry[12]
                # if ( self.atom_rads.has_key(atype)): #python2
                if ( atype in self.atom_rads): # python3
                    self.radius.append(self.atom_rads[atype])
                else:
                    print ("atom radius not in dictionary", atype)
                    self.atom_rads[atype] = 0.0
                    self.radius.append(0.0)
    def readligand(self,pdb_name,handle=False):
        #
        # read multiple entries, delimited by MODEL/ENDMDL records
        #
        try:
            pdb_file = open(pdb_name,"r") # python3
            #pdb_file = file(pdb_name,"r")  # python2
            #print 'opening file'
        except:
            # print "can't find file", pdb_name
            sys.exit("can't find file")
        contents = pdb_file.readlines()
        pdb_file.close()
        #print(type(contents))
        #print('lines read: ', len(contents))
        self.natom = 0
        self.nmodel = 0
        xyz = [0.,0.,0.]
        nline = 0
        while(nline < len(contents)):
          entry = contents[nline]
          while(entry[0:5] != 'MODEL'):
            nline += 1
            if(nline == len(contents)):break
            entry = contents[nline]
          if(nline == len(contents)):break
          #print('found model: ',nline)
          self.nmodel += 1
          #
          # if first model, create temporary pdb file for pymol
          if(self.nmodel == 1):
            tmp_pdb = open('dockeye_lig.pdb','w')
            gcen = [0.,0.,0.]
            dcen = 2.
          nline += 1
          entry = contents[nline]
          #print(entry)
          while(entry[0:6] != 'ENDMDL'):
            if (entry[0:4] == 'ATOM') or (entry[0:6] == 'HETATM'):
              #
              xyz[0] = float(entry[30:38])
              xyz[1] = float(entry[38:46])
              xyz[2] = float(entry[46:54])
              self.coords.append([xyz[0],xyz[1],xyz[2]])
              # only do names, charge, radius for first model
              if(self.nmodel == 1):
                # this code will write ligand conformer #1 as graphical manipulation 'handle'
                if(not handle): tmp_pdb.write(entry)  
                for k in range(3):
                  gcen[k] += xyz[k]
                self.natom +=1
                self.occ.append(float(entry[54:60]))
                self.bfact.append(float(entry[61:67]))
                atname = entry[11:17]
                self.name.append(atname)
                self.res.append(entry[17:21])
                self.chain.append(entry[21:22])
                self.resnum.append(entry[22:26])
                if(entry[12] != 'H'):
                  atype = entry[12:14]
                  #atype = ' ' + entry[12]
                else:
                  atype = ' ' + entry[12]
                # if ( self.atom_rads.has_key(atype)): #python2
                if ( atype in self.atom_rads): # python3
                  self.radius.append(self.atom_rads[atype])
                else:
                  print ("atom radius not in dictionary", atype)
                  self.atom_rads[atype] = 0.0
                  self.radius.append(0.0)
              # end if model 1
            # end if atom entry
            nline += 1
            entry = contents[nline]
          # end of model
          if((self.nmodel == 1) and (handle)):
            # this code will generate a 3-d cross as a graphical manipulation handle
            print('making handle...')
            for k in range(3): # create molecular handle for rot/trans
              gcen[k] /= self.natom
            strhead = 'HETATM    1  C1  HND A   1    '
            strtail = '  1.00  0.00           C'
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[0] += dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[0] -= 2.*dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[0] += dcen
            gcen[1] += dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[1] -= 2.*dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[1] += dcen
            gcen[2] += dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            gcen[2] -= 2.*dcen
            tmp_pdb.write("%s%8.3f%8.3f%8.3f%s\n" % (strhead,gcen[0],gcen[1],gcen[2],strtail))
            # end of code for a 3-d cross as a graphical manipulation handle
        tmp_pdb.close()
        # end of file
        if(self.nmodel == 0):
          print('ERROR: ligand file must have at least one conformer')
          print('bracketed by MODEL, ENDMDL records')
          sys.exit()
        else:
          print('# of models: ',self.nmodel)
          nrec = len(self.coords)
          nrec_need = self.nmodel*self.natom
          if(nrec == nrec_need):
            print('found all coordinates: ',nrec_need)
          else:
            print('mismatch between expected & found # of coords ',nrec_need,nrec)
            sys.exit()
#=======================================
def rot_vec(rmt,vec,inv=0):
  vec_rot = [0.,0.,0.]
  if(inv == 0):
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[i][j]*vec[j]
  else:
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[j][i]*vec[j]
  return vec_rot

def vdot(v1,v2):
  dot = 0.
  for k in range(3):
    dot += v1[k]*v2[k]
  return dot
#
def vcross(v1,v2):
  cross = [0.,0.,0.]
  cross[0] = v1[1]*v2[2] - v1[2]*v2[1]
  cross[1] = v1[2]*v2[0] - v1[0]*v2[2]
  cross[2] = v1[0]*v2[1] - v1[1]*v2[0]
  return cross
#
def vnorm(v1):
  dot = vdot(v1,v1)
  dot = math.sqrt(dot)
  if(dot > 0.):
    for k in range(3):
      v1[k] /= dot
  return dot
#
def vperp(v1):
  v2 = [v1[1],v1[2],v1[0]]
  v3 = vcross(v1,v2)
  return v3
#=======================================
