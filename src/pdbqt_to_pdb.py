"""
read pdbqt format coordinate file output from autodock tools,
move autdock charges into bfactor field, and assign radii
"""
import sys
#===============================
atom_rads = {' C':' 1.80',' S':' 1.90',' O':' 1.60',' N':' 1.40',
        ' P':' 1.80',' H':' 0.00','ZN':' 1.40',
        ' Z':' 1.40',' B':' 2.46',' F':' 1.56'}
#===============================
print('\n read pdbqt format coordinate file output from autodock tools,')
print(' move autdock charges into bfactor field, and assign radii\n')
if(len(sys.argv) <2):
  print('Usage: python pdbqt_to_pdb.py autodock_pdbqt_file')
  sys.exit()
pdb_in = sys.argv[1]
pdbqt_file = open(pdb_in,'r')
contents = pdbqt_file.readlines()
pdb_out = pdb_in[:-6] + '_qr.pdb'
#print(pdb_out)
pdb_file = open(pdb_out,'w')
pdb_file.write('MODEL\n')
nrec = 0
netcrg = 0.
for entry in contents:
  if((entry[0:6] == 'HETATM') or (entry[0:6] == 'ATOM  ')):
    nrec += 1
    if(entry[12] != 'H'):
      atype = entry[12:14]
    else:
      atype = ' ' + entry[12]
    if ( atype in atom_rads): # python3
      atrad = atom_rads[atype]
    else:
      print ("atom radius not in dictionary", atype)
      atom_rads[atype] = ' 0.00'
      atrad = ' 0.00'
    charge = entry[69:76]
    netcrg += float(charge)
    entry1 = entry[0:55] + atrad + charge + '\n'
    pdb_file.write(entry1)
  else:
    pdb_file.write(entry)
pdb_file.write('ENDMDL\n')
pdbqt_file.close()
print('atom records read: ',nrec,' net charge: ',netcrg)
