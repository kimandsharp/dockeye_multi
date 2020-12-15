/*
  Kim Sharp, 1/19/2014, C wrapper called by dockeyeM_c.py 
  much of this taken verbatim from:
  http://docs.python.org/2/extending/extending.html#a-simple-example
  nov 18th, 2020: return translational and torque forces
*/

#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
int MAXATOMP = 12000;
int MAXATOML = 300;
int MAXCIRCLE = 500;
int NUM_POINTS = 10;
int MAXDATA = 70000; /* should be at least 12*MAXCIRCLE*NUM_POINTS+5 */
int MAXMOD = 729;
/* pymol cgo keywords */
float BEGIN = 2.0;
float END = 3.0;
float VERTEX = 4.0;
float POINTS =  0.0;
float LINE = 1.0;
float LINES = 1.0;
float TRIANGLES =  4.0;
float COLOR =  6.0;
float LINEWIDTH =  5.0;
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
static PyObject * energy_c(PyObject *self, PyObject *args)
{
    PyObject* atom_data;
    PyObject* disp_obj = PyList_New(MAXDATA+5); /* list object to pass back results */
    int i,j,k,n1,n2,nmod,nm;
    int ioff;
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
    float rcut2,d2,efact,dd,sigma,vfact;
    float rr,rr2,rr6,rr12,evdw,eelect,rr9;
    float e_at_max,e_at_tot,e_at,e_abs,color;
    float radius,pi,tpi,angle,ca,sa;
    float dsd[MAXDATA];
    float vbeg[3],vend[3],vp[3],v1[3],v2[3],rgb[3];
    float ftrans[3],frot[3],gcen[3],fat[3],dr[3],dg[3],df[3];
    /*
      parse the calling arguments- lists etc are read as python objects by format 'O'
    */
    if (!PyArg_ParseTuple(args, "O", &atom_data)){
        return NULL;
    }
    /*
    get # of atoms stored as floats in 1st two elements
    and # of models
    */
    ioff = 0;
    n1 = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
    ioff++;
    n2 = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
    ioff++;
    nmod = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
    ioff++;
    // printf("# of models %d \n ",nmod);
    if (nmod > MAXMOD){
         printf ("too many models - will only do %d \n",MAXMOD);
         nmod = MAXMOD;
    }

    if (n1 > MAXATOMP){
         printf ("# atoms %d > MAXATOMP %d: will truncate \n",n1, MAXATOMP);
         n1 = MAXATOMP;
    }
    if (n2 > MAXATOML){
         printf ("# atoms %d > MAXATOML %d: will quit \n",n2, MAXATOML);
         exit(1);
    }
    /*
    get rest of data
    */
    for (i=0 ; i<n1; i++){
      rad1[i] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
      ioff++;
    }
    for (i=0 ; i<n1; i++){
      q1[i] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
      ioff++;
    }
    for (i=0 ; i<n1; i++){
      for (k=0; k<3; k++){
        crd1[i][k] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
        ioff++;
      }
    }
    ioff = 5*n1+3; // because we might have truncated # of atoms in molc 1
    for (i=0 ; i<n2; i++){
      rad2[i] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
      ioff++;
    }
    for (i=0 ; i<n2; i++){
      q2[i] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
      ioff++;
    }
    for (nm=0 ; nm<nmod; nm++){
      for (i=0 ; i<n2; i++){
        for (k=0; k<3; k++){
          crd2[nm][i][k] = PyFloat_AsDouble(PyList_GET_ITEM(atom_data, (Py_ssize_t)ioff));
          ioff++;
        }
      }
    }
    // printf ("# atoms %d %d ioff %d \n",n1,n2,ioff);
    // printf(" x1  %f xN %f \n",crd1[0][0],crd1[n1-1][2]);
    // printf(" x1  %f xN %f \n",crd2[0][0],crd2[n2-1][2]);
    // printf(" r1  %f rN %f \n",rad1[0],rad2[n2-1]);
    // printf(" q1  %f qN %f \n",q1[0],q2[n2-1]);
    /*
    calculate energies
    */
    rcut2 = RCUT*RCUT;
    efact = 332./DIEL;
    vfact = 27./4.;  // for 9-6 potl.
    // start model selection
    et_best = 1.e6;
    nbest = 0;
    for (nm=0 ; nm<nmod; nm++){
      et = 0.;
      for (i=0;i<n1;i++){
        for (j=0;j<n2;j++){
          d2 = 0.;
          for (k=0;k<3;k++){
            dxyz[k] = (crd1[i][k] - crd2[nm][j][k]);
            d2 = d2 + dxyz[k]*dxyz[k];
          }
          if (d2 < rcut2){
            sigma = rad1[i] + rad2[j];
            dd = sqrt(d2) + 1.e-3;
            rr = sigma/dd;
            rr2 = rr*rr;
            rr6 = rr2*rr2*rr2;
            //
            // 12-6 potl
            //
            // rr12 = rr6*rr6;
            // evdw = 4.*EPS*(rr12 - rr6);
            //
            // 9-6 potl
            //
            rr9 = rr6*rr2*rr;
            evdw = vfact*EPS*(rr9 - rr6);
            eelect = efact*q1[i]*q2[j]/dd;
            et = et + evdw + eelect;
          }
        }
      }
      if(et < et_best){
        et_best = et;
        nbest = nm;
        // printf("best model %d %f \n",nbest,et_best);
      }
    }
    nm = nbest;
    // printf("using model %d \n",nm);
    // end model selection
    ninter = 0;
    et = 0.;
    ee = 0.;
    ev = 0.;
    for (k=0;k<3;k++){
      dxyz_at[k] = 0;
      ftrans[k] = 0.;
      frot[k] = 0.;
      gcen[k] = 0.;
    }
    // for net force, torque
    for (j=0;j<n2;j++){
      for (k=0;k<3;k++){
        gcen[k] = gcen[k] + crd2[nm][j][k];
      }
    }
    for (k=0;k<3;k++){
      gcen[k] = gcen[k]/n2;
    }
    for (i=0;i<n1;i++){
      e_at_max = 0.;
      e_at_tot = 0.;
      ninter_at = 0;
      for (j=0;j<n2;j++){
        d2 = 0.;
        for (k=0;k<3;k++){
          dxyz[k] = (crd1[i][k] - crd2[nm][j][k]);
          d2 = d2 + dxyz[k]*dxyz[k];
        }
        if (d2 < rcut2){
          ninter_at++;
          sigma = rad1[i] + rad2[j];
          dd = sqrt(d2) + 1.e-3;
          rr = sigma/dd;
          rr2 = rr*rr;
          rr6 = rr2*rr2*rr2;
          // for net force, torque
          for (k=0;k<3;k++){
            dg[k] = (crd2[nm][j][k] - gcen[k]);
            dr[k] = dxyz[k];
          }
          vnorm(dr);
          //
          // 12-6 potl
          //
          // rr12 = rr6*rr6;
          // evdw = 4.*EPS*(rr12 - rr6);
          //
          // 9-6 potl
          //
          rr9 = rr6*rr2*rr;
          evdw = vfact*EPS*(rr9 - rr6);
          ev = ev + evdw;
          eelect = efact*q1[i]*q2[j]/dd;
          ee = ee + eelect;
          e_at = eelect + evdw;
          e_at_tot = e_at_tot + e_at;
          // for net force, torque
          for (k=0;k<3;k++){
            fat[k] =  - dr[k]*(vfact*EPS*(9.*rr9/dd - 6.*rr6/dd) +efact*q1[i]*q2[j]/(dd*dd));
            ftrans[k] = ftrans[k] + fat[k];
          }
          vcross(fat,dg,df);
          for (k=0;k<3;k++){
            frot[k] = frot[k] + df[k];
          }
          /* store max interaction for atom i */
          if (fabs(e_at) >= fabs(e_at_max)){
            e_at_max = e_at;
            for (k=0;k<3;k++){
              xyzmid[k] = crd2[nm][j][k] + dxyz[k]/2.;
              dxyz_at[k] = dxyz[k];
            }
          }
        }
      }
      if (ninter_at > 0){
        d2 = 0;
        for (k=0;k<3;k++){
          d2 = d2 + dxyz_at[k]*dxyz_at[k];
        }
        dd = sqrt(d2);
        for (k=0;k<3;k++){
          dxyz_at[k] = dxyz_at[k]/dd;
        }
        // printf(" dxyz_at  %f %f %f \n",dxyz_at[0],dxyz_at[1],dxyz_at[2]);
        e_abs = fabs(e_at_tot)/2.;
        if (e_at_tot >= 0.){
          color = 0.7 + MIN(0.3,e_abs);
        }else{ 
          color = 0.3 - MIN(0.3,e_abs);
        }
        // printf("e_abs, e_at_tot, color: %f %f %f \n",e_abs,e_at_tot,color);
        radius = 1. - dd/RCUT;
        for (k=0;k<3;k++){
          circle_mid[ninter][k] = xyzmid[k];
          circle_perp[ninter][k] = dxyz_at[k];
        }
        circle_color[ninter] = color;
        circle_rad[ninter] = radius;
        ninter++;
      }
    }
    et = ee + ev;
    /* printf("%d interactions \n",ninter); */
    /* printf("ftrans %f %f %f \n",*(ftrans+0),*(ftrans+1),*(ftrans+2)); */
    /* printf("frot   %f %f %f \n",*(frot+0),*(frot+1),*(frot+2)); */
    ndata = 12*NUM_POINTS*ninter + 5;
    if(ndata > MAXDATA){
      printf("warning: too many circles, increase maxdata \n");
      ninter = MAXDATA/(12*NUM_POINTS)-1;
      printf("# of interactions reduced to %d \n",ninter);
    }
    /* printf(" Ee  %g Ev %g Et %g \n",ee,ev,et); */
    /*
    generate display object
    */
    // header
    ndata = 0;
    dsd[ndata] = LINEWIDTH;
    ndata++;
    dsd[ndata] = 4.0; // linewidth
    ndata++;
    dsd[ndata] = BEGIN;
    ndata++;
    // printf("LINE %8.3f \n",LINE);
    // printf("LINES %8.3f \n",LINES);
    dsd[ndata] = LINE;
    ndata++;
    tpi = 2.*3.14159;
    for (i=0;i<ninter;i++){
      for (k=0;k<3;k++){
        vp[k] = circle_perp[i][k];
      }
      // printf(" vp  %f %f %f \n",vp[0],vp[1],vp[2]);
      vperp(vp,v1);
      vnorm(v1);
      // printf(" v1  %f %f %f \n",v1[0],v1[1],v1[2]);
      vcross(v1,vp,v2);
      // printf(" v2  %f %f %f \n",v2[0],v2[1],v2[2]);
      // printf(" vmid, rad  %f %f %f %f \n",circle_mid[i][0],circle_mid[i][1],circle_mid[i][2],circle_rad[i]);
      // printf(" color  %f \n",circle_color[i]);
      for (k=0;k<3;k++){
        vbeg[k] = circle_rad[i]*v1[k] + circle_mid[i][k];
      }
      for (j=1;j<NUM_POINTS+1;j++){
        angle = j*tpi/NUM_POINTS;
        ca = cos(angle);
        sa = sin(angle);
        for (k=0;k<3;k++){
          vend[k] = circle_rad[i]*(ca*v1[k] + sa*v2[k]) + circle_mid[i][k];
        }
        color_map(circle_color[i],rgb);
        dsd[ndata] = COLOR;
        ndata++;
        for (k=0;k<3;k++){
          dsd[ndata] = rgb[k];
          ndata++;
        }
        dsd[ndata] = VERTEX;
        ndata++;
        for (k=0;k<3;k++){
          dsd[ndata] = vbeg[k];
          ndata++;
        }
        dsd[ndata] = VERTEX;
        ndata++;
        for (k=0;k<3;k++){
          dsd[ndata] = vend[k];
          ndata++;
          vbeg[k] = vend[k];
        }
      }
    }
    dsd[ndata] = END;
    ndata++;
    // put energy terms at end
    dsd[ndata] = et;
    ndata++;
    dsd[ndata] = ee;
    ndata++;
    dsd[ndata] = ev;
    ndata++;
    dsd[ndata] = nbest;
    ndata++;
    // put force/torque terms at the end
    dsd[ndata] = ftrans[0];
    ndata++;
    dsd[ndata] = ftrans[1];
    ndata++;
    dsd[ndata] = ftrans[2];
    ndata++;
    dsd[ndata] = frot[0];
    ndata++;
    dsd[ndata] = frot[1];
    ndata++;
    dsd[ndata] = frot[2];
    ndata++;
    dsd[0] = ndata; // we temp overwrite LINEWIDTH with ndata, 
    // printf ("display object length %d \n",ndata);
    j = 0;
    /*
    for (i=0;i<ndata;i++){
      printf(" %8.3f ",dsd[i]);
      j++;
      if(j = 10){
        printf("\n");
        j = 0;
      }
    }
    */

    /* construct Python list object to pass back results */
    for (i=0;i<ndata; i++){
      PyObject *num = PyFloat_FromDouble(dsd[i]);
      if (!num){
        Py_DECREF(disp_obj);
        return NULL;
      }
      PyList_SET_ITEM(disp_obj, i, num);
    }
    return disp_obj;
    /* example from web:
    for (i=1;i<6; i++){
      array[i] = 1.*i;
      PyObject *num = PyFloat_FromDouble(array[i]);
      if (!num){
        Py_DECREF(lst);
        return NULL;
      }
      PyList_SET_ITEM(lst, i, num);
    }
    return lst;
    */
    /* 
    finished
    */
    /* return Py_BuildValue("l", MAXATOM); */
}

/* create Python's methods table. 2nd argument must be same as
   that of the PyObject defined above */
static PyMethodDef energy_methods[] = {
    {"energy_c", energy_c, METH_VARARGS,"energy_c"},
    {NULL, NULL,0,NULL}
};
/* call to initialize- the function must have the exact module name 
   (i.e. this file's name) preceded by string 'init' 
   and the 1st argument to Py_InitModule must be the exact module name too.
   2nd argument must be same as the PyMethodDef function above */
PyMODINIT_FUNC initdockeyeM_energy(void)
{
    Py_InitModule("dockeyeM_energy", energy_methods);
}

