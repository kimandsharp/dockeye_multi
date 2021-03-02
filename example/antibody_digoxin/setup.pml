#run $PYTHONPATH/dockeyeM_c.py
run $HOME/source/dockeye_multi/src/dockeyeM_c.py
de('ab.pdb','ag_tor.pdb',forcerep=False)
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
show surface, dockeye_prt
set transparency, 0.4
