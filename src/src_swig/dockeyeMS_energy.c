/*
  reimplementation of dockeyeM_energy.c to use swig to generate
  pyton-C interface with dockeyeMS_c.py
  kas 13 feb 2020
*/

#include <stdio.h>
#include <math.h>
//===============================================
// MAIN
//===============================================

/*
    this is what will be 'called' from Python program as energy_c(atom_data)
*/
double energy_c(int n1, int n2, int nmod, double *atom_data){
    double x;
    printf("# atoms 1 %d \n ",n1);
    printf("# atoms 2 %d \n ",n2);
    printf("# models  %d \n ",nmod);
    printf("atom      %f \n ",*atom_data);
    x = 4.*atom_data[2];
    return x;
    /*
    int i,j,k,nm;
    float *ioff;
    float array[6];
    float crd1[MAXATOMP][3];
    float crd2[MAXMOD][MAXATOML][3];
    float rad1[MAXATOMP];
    float rad2[MAXATOML];
    float q1[MAXATOMP];
    float q2[MAXATOML];
    float circle_mid[MAXCIRCLE][3],circle_perp[MAXCIRCLE][3];
    float circle_color[MAXCIRCLE],circle_rad[MAXCIRCLE];
    int ninter,ninter_at,ndata;
    float et,ee,ev;
    float et_best;
    int nbest;
    float dxyz[3],xyzmid[3],dxyz_at[3];
    float rcut2,d2,efact,dd,sigma;
    float rr,rr2,rr6,rr12,evdw,eelect,rr9;
    float e_at_max,e_at_tot,e_at,e_abs,color;
    float radius,pi,tpi,angle,ca,sa;
    float dsd[MAXDATA];
    float vbeg[3],vend[3],vp[3],v1[3],v2[3],rgb[3];
    */
}
