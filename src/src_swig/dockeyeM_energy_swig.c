/*
  Kim Sharp, 1/19/2014, taken
  from C wrapper called by dockeyeM_c.py 
  now pruned back to to basic c for use by swig- kas 14 feb 2020
*/

#include <stdio.h>
#include <math.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
int MAXATOMP = 6000;
int MAXATOML = 6000;
int MAXCIRCLE = 500;
int NUM_POINTS = 10;
int MAXDATA = 70000; /* should be at least 12*MAXCIRCLE*NUM_POINTS+5 */
int MAXMOD = 100;
/* pymol cgo keywords */
float BEGIN = 2.0;
float END = 3.0;
float VERTEX = 4.0;
float POINTS =  0.0;
float LINE = 1.0;
float LINES = 1.0;
float TRIANGLES =  4.0;
float COLOR =  6.0;
float LINEWIDTH =  10.0;
float RCUT=5., EPS=0.1, DIEL=80.;
//===============================================
void color_map(float color, float *rgb){
    *(rgb+0) = 1.;
    *(rgb+1) = 1.;
    *(rgb+2) = 1.;
    if((color >= 0.) && (color < 0.5)){
      *(rgb+0) = 0.;
      *(rgb+1) = 2.*color;
      *(rgb+2) = (1. - 2.*color);
    }else{
      if ((color >= 0.5) && (color <= 1.0)){
        *(rgb+0) = 2.*color - 1.;
        *(rgb+1) = 2. - 2.*color;
        *(rgb+2) = 0.;
      }
    }
}
float vdot(float *v1, float *v2){
    float rdot;
    int k;
    rdot = 0.;
    for (k=0;k<3;k++){
      rdot = rdot + *(v1+k)* *(v2+k);
    }
    return rdot;
}
void vcross(float *v1, float *v2, float *v3){
    *(v3+0) = *(v1+1)* *(v2+2) - *(v1+2)* *(v2+1);
    *(v3+1) = *(v1+2)* *(v2+0) - *(v1+0)* *(v2+2);
    *(v3+2) = *(v1+0)* *(v2+1) - *(v1+1)* *(v2+0);
    /*
    printf("v1 in vcross %f %f %f \n",*(v1+0),*(v1+1),*(v1+2));
    printf("v2 in vcross %f %f %f \n",*(v2+0),*(v2+1),*(v2+2));
    printf("v3 in vcross %f %f %f \n",*(v3+0),*(v3+1),*(v3+2));
    */
}
void vperp(float *v1, float *v3){
    float v2[3];
    v2[0] = *(v1+1);
    v2[1] = *(v1+2);
    v2[2] = *(v1+0);
    // printf("in vperp %f %f %f \n",v2[0],v2[1],v2[2]);
    vcross(v1,v2,v3);
    // printf("v3 in vperp %f %f %f \n",*(v3+0),*(v3+1),*(v3+2));
}
void vnorm(float *v1){
    float dot;
    int k;
    dot = vdot(v1,v1);
    dot = sqrt(dot);
    if(dot > 0.){
      for (k=0;k<3;k++){
        *(v1+k) = *(v1+k)/dot;
      }
    }

}
//===============================================
// MAIN
//===============================================


/*
    this is what will be 'called' from Python program as energy_c(atom_data)
*/
static float * energy_c(float *atom_data)
{
    int i,j,k,n1,n2,nmod,nm;
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
    /*
    get # of atoms stored as floats in 1st two elements
    and # of models
    */
    ioff = atom_data;
    n1 = *ioff;
    *ioff++;
    printf("# atoms 1 %d \n ",n1);
    /*
    n2 = *ioff;
    *ioff++;
    printf("# atoms 1 %d \n ",n2);
    nmod = *ioff
    *ioff++;
    printf("# of models %d \n ",nmod);
    */
}
