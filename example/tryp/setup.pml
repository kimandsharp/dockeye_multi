#------------------------------------------------
run $HOME/source/dockeye_multi/src/dockeyeM_c.py
de("de_prot_qr.pdb","de_ligand_qr_torf.pdb")
#optional view settings
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------

