# rules
.f.o:
	$(FCOMP) -c $(FFLAGS) $*.f

build:
	swig -python dockeyeMS_energy.i
	gcc -fpic -c dockeyeMS_energy.c dockeyeMS_energy_wrap.c
	gcc -shared dockeyeMS_energy.o dockeyeMS_energy_wrap.o -o _dockeyeMS_energy.so
	#gcc -m32 -shared -s dockeyeM_energy.o -o dockeyeMS.so
	#ln  -s -f dockeyeM_energy_linux14_python27_32bit.so dockeyeM_energy.so
