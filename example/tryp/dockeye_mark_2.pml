#------------------------------------------------
#run $PYTHONPATH/dockeyeM_c.py
de("de_prot_qr.pdb","dockeye_lig_mark_2.pdb")
#optional view settings
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------
