# eventually figure out correct lib-dynload directory on each machine to insert link to *.so
# and each PyMol inst.
#DYNLOAD = /home/kim/Downloads/pymol/ext/lib/python2.4/lib-dynload
# cd /Applications/PyMOL.app/Contents/lib/python2.7/lib-dynload
# better: set PYTHONPATH=/home/kim/source/dockeye_multi/src
# feb 23 this src now in subdirectory src/src_c
#
MODULE = ../dockeyeM_energy.so
PYFILE = ../dockeyeM_c.py
SRC_C = src_c


# rules
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f

#
# home computer
home:
	gcc -O2 -fPIC -m32 -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/home/kim/anaconda2/include/python2.7/
	gcc -m32 -shared -s dockeyeM_energy.o -lm  -L/home/kim/anaconda2/lib/python2.7/config/ -lpython2.7 -u MAIN__ -o dockeyeM_energy_linux14_python27_32bit.so
	ln  -s -f $(SRC_C)/dockeyeM_energy_linux14_python27_32bit.so $(MODULE)
	ln -s -f $(SRC_C)/dockeyeM_c_v1.5.py $(PYFILE)
#
# mac book air
macbook_fink:
	gcc -O2 -shared -I/sw/include/python2.7 -lpython2.7 -L/sw/lib/python2.7/config/ -o dockeyeM_energy_MacOS_python27_32bit.so dockeyeM_energy.c
	ln  -s -f $(SRC_C)/dockeyeM_energy_MacOS_python27_32bit.so $(MODULE)
	ln -s -f $(SRC_C)/dockeyeM_c_v1.5.py $(PYFILE)
#
# mac book air with schrodinger PyMol
macbook:
	gcc -O2 -shared -I/sw/include/python2.7 -lpython2.7 -L/Applications/PyMOL_py27.app/Contents/lib/ -o dockeyeM_energy_MacOS_python27_64bit.so dockeyeM_energy.c
	ln  -s -f $(SRC_C)/dockeyeM_energy_MacOS_python27_64bit.so $(MODULE)
	ln -s -f $(SRC_C)/dockeyeM_c_v1.5.py $(PYFILE)
#
# office
office:
	gcc -O2 -fPIC -m32 -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/home/sharp/anaconda2/include/python2.7/
	gcc -m32 -shared -s dockeyeM_energy.o -lm  -L/home/sharp/anaconda2/lib/python2.7/config/ -lpython2.7 -u MAIN__ -o dockeyeM_energy_linux16_python27_32bit.so
	ln  -s -f $(SRC_C)/dockeyeM_energy_linux16_python27_32bit.so $(MODULE)
	ln -s -f $(SRC_C)/dockeyeM_c_v1.5.py $(PYFILE)
#
# redtrooper
redtrooper:
	gcc -O2 -fPIC -c dockeyeM_energy.c -o dockeyeM_energy.o -I/usr/include -I/usr/include/python2.7/
	gcc -shared -s dockeyeM_energy.o -lm  -Lpython2.7 -lpython2.7 -u MAIN__ -o dockeyeM_energy_redtrooper_python27_64bit.so
	ln  -s -f $(SRC_C)/dockeyeM_energy_redtrooper_pythonn27_64bit.so $(MODULE)
	ln -s -f $(SRC_C)/dockeyeM_c_v1.5.py $(PYFILE)
