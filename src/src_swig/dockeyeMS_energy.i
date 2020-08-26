/* swig interface file */
%module dockeyeMS_energy

%inline %{
extern double energy_c(int n1, int n2, int nmod, double *atom_data);
%}

