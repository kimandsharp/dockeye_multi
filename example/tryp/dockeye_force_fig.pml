#------------------------------------------------
run $HOME/source/dockeye_multi/src/dockeyeM_c.py
de("de_prot_qr.pdb","dockeye_lig_force_fig.pdb",forcerep=True)
#optional view settings
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------
