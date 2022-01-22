#------------------------------------------------
#set PYTHONPATH=( $HOME/source/dockeye_multi/src )
#run $PYTHONPATH/dockeyeM_c.py
run $HOME/source/dockeye_multi/src/dockeyeM_c.py
de("1cnx_prot_qr.pdb","1cnx_eg2_qr_cnf.pdb")
#de("1cnx_prot_qr.pdb","1cnx_eg2_qr_cnf_pose_67.pdb")
#optional view settings
hide lines
spectrum b, red_white_blue
select site, resi  92+94+96+106+119+121+131+135+143+198+200+204+209
show sticks, site
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
#------------------------------------------------
