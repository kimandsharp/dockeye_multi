#------------------------------------------------
#run $PYTHONPATH/dockeyeM_c.py
de("ab.pdb","dockeye_lig_mark_1.pdb")
#optional view settings
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------
