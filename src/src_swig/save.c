    if (nmod > MAXMOD){
         printf ("too many models - will only do %d \n",MAXMOD);
         nmod = MAXMOD;
    }

    if (n1 > MAXATOMP){
         printf ("# atoms %d > MAXATOMP %d: will truncate \n",n1, MAXATOMP);
         n1 = MAXATOMP;
    }
    if (n2 > MAXATOML){
         printf ("# atoms %d > MAXATOML %d: will truncate \n",n2, MAXATOML);
         n2 = MAXATOML;
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
            evdw = 27.*EPS*(rr9 - rr6)/4.;
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
          //
          // 12-6 potl
          //
          // rr12 = rr6*rr6;
          // evdw = 4.*EPS*(rr12 - rr6);
          //
          // 9-6 potl
          //
          rr9 = rr6*rr2*rr;
          evdw = 27.*EPS*(rr9 - rr6)/4.;
          ev = ev + evdw;
          eelect = efact*q1[i]*q2[j]/dd;
          ee = ee + eelect;
          e_at = eelect + evdw;
          e_at_tot = e_at_tot + e_at;
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
}
