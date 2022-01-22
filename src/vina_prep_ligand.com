#!/bin/csh
# produces pdb file with torsion records, autodock atom types for vina
#=================================
if ($#argv == 0) then
  echo "usage vina_prep.com pdbfile"
  exit
else
  set pdbfile=$1
  echo $pdbfile
endif
preplig:
set pyshell='~/MGLTools-1.5.6/bin/python2.5'
set preplig='~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
# don't merge apolar H atoms: need to keep same # of atoms as am1bcc/Bqequil files!
$pyshell $preplig -v -U '' -l $pdbfile
