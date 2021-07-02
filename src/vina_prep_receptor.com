#!/bin/csh
# produces pdb file with torsion records, autodock atom types for vina
#=================================
if ($#argv == 0) then
  echo "usage vina_prep_receptor.com pdbfile"
  exit
else
  set pdbfile=$1
  echo $pdbfile
endif
preplig:
set pyshell='~/MGLTools-1.5.6/bin/python2.5'
set prepprot='~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'
$pyshell $prepprot -v  -A hydrogens -r $pdbfile
