#------------------------------------------------
run $HOME/source/dockeye_multi/src/dockeyeM_c.py
de("1cnx_prot_qr.pdb","dockeye_lig_mark_3.pdb")
#optional view settings
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------
