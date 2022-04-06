FFLAGS = -g -C -Wall -ffixed-line-length-132
MODULE = ../dockeyeM_energy.so
PYFILE = ../dockeyeM_c.py
#F2PY = /home/sharp/miniconda3/bin/f2py
F2PY = /home/sharp/Downloads/PyMOL-2.4.1_198-Linux-x86_64-py37/pymol/bin/f2py_local
SRC_F = src_f
#CWD = `pwd`
#$(CWD)
#
# check fortran code for compile errors
code:
	gfortran $(FFLAGS) -c dockeyeM_energy.f
# use simple f2py make
main:
	$(F2PY)  -c --f77flags="-O" dockeyeM_energy.f -m dockeyeM_energy
	#python -m numpy.f2py -c dockeyeM_energy.f -m dockeyeM_energy
	ln -s -f $(SRC_F)/dockeyeM_f_v1.4.py $(PYFILE)
	# macos
	#ln -s -f $(SRC_F)/dockeyeM_energy.cpython-37m-darwin.so $(MODULE)
      # linux
	ln -s -f $(SRC_F)/dockeyeM_energy.cpython-37m-x86_64-linux-gnu.so $(MODULE)

# use  smarter method with in/out intent statements in *.pyf file- not working right now
# but is it really needed?
smart:
	#python -m numpy.f2py dockeyeM_energy.f -m dockeyeM_energy -h dockeyeM_energy.pyf
	python -m numpy.f2py -c dockeyeM_energy.pyf dockeyeM_energy.f
