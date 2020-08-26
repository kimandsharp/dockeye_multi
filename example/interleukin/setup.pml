run $PYTHONPATH/dockeyeM_c.py
de('IL1B.pdb','MIM_tor.pdb')
hide lines
spectrum b, red_white_blue
show sticks, dockeye_lig
select site, resi 18+39+40+149
show sticks, site
show surface, dockeye_prt
set transparency, 0.4
