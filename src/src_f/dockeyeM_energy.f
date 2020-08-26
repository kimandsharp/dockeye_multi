c
c implementation of dockeyeM_energy.c in fortran
c interfaced using numpy.f2py
c
      subroutine energy_f(energy_obj,atom_data,n)
      implicit none
      integer  maxobjdata,maxcircle
      integer num_points
      parameter (maxobjdata = 20000)
      parameter (maxcircle = 500)
      parameter (num_points = 10)
      real*8 atom_data(n)
      integer n,i,j,k,m
c
      integer n1, n2,nmod
      integer indqi,indqj,indri,indrj,indci,indcj
      integer m_best,ninter,ninter_at,ndata
c
      real*8 rcut2, efact, vfact, et, ev, ee, et_best
      real*8 e_at,e_at_tot,e_at_max,e_abs
      real*8 color,radius
      real*8 q1,q2,r1,r2,d2,rr,rr2,rr6,rr9
      real*8 x1(3),x2(3)
      real*8 sigma,dd,evdw,eelect
      real*8 dxyz(3),xyzmid(3),dxyz_at(3)
      real*8 circle_mid(3,maxcircle),circle_perp(3,maxcircle)
      real*8 circle_color(maxcircle),circle_rad(maxcircle)
      real*8 energy_obj(maxobjdata),vp(3),v1(3),v2(3),angle,ca,sa
      real*8 vbeg(3),vend(3),rgb(3)
      real*8 vdot
c
      real*8 RCUT, EPS, DIEL,tpi
      real*8 BEGIN,END,VERTEX,LINE,KOLOR,LINEWIDTH
      data RCUT, EPS, DIEL / 5., 0.1, 80. /
      data BEGIN,END,VERTEX,LINE,KOLOR,LINEWIDTH 
     &     / 2.0, 3.0, 4.0, 1.0, 6.0, 10.0 /
      data tpi / 6.2832 /
c========================================

      n1 = int(atom_data(1))
      n2 = int(atom_data(2))
      nmod = int(atom_data(3))
      print *,'n1: ',n1,n2,nmod
c
c calculate energies find best model
c
      rcut2 = RCUT*RCUT
      efact = 332./DIEL
      vfact = 27.*EPS/4.
      et_best = 1.e8
      m_best = 0
      do m = 0,nmod-1
        et = 0.
        do j = 1,n2
          indrj = 3 + 5*n1 + j
          indqj = 3 + 5*n1 + n2 + j
          indcj = 3 + 5*n1 + 2*n2 + 3*m*n2 + 3*(j-1)
          r2 = atom_data(indrj)
          q2 = atom_data(indqj)
          do k = 1,3
            indcj = indcj + 1
            x2(k) = atom_data(indcj)
          end do
c          print *,m,j,indcj,r2,q2,x2
          do i =1,n1
            indri = 3 + i
            indqi = 3 + n1 + i
            indci = 3 + 2*n1 + 3*(i-1)
            r1 = atom_data(indri)
            q1 = atom_data(indqi)
            do k = 1,3
              indci = indci + 1
              x1(k) = atom_data(indci)
            end do
            d2 = 0.
            do k = 1,3
              d2 = d2 + (x1(k) - x2(k))**2
            end do
            if (d2 .lt. rcut2)then
              sigma = r1 + r2
              dd = sqrt(d2) + 1.e-3
              rr = sigma/dd
              rr2 = rr*rr
              rr6 = rr2*rr2*rr2
              rr9 = rr6*rr2*rr
              evdw = vfact*(rr9 - rr6)
              eelect = efact*q1*q2/dd
              et = et + evdw + eelect
            end if
          end do
        end do
        if(et < et_best)then
          et_best = et
          m_best = m
        end if
      end do
      print *,'best model, energy: ',m_best, et_best
c
c create dockeye object from best model
c
      ninter = 0
      et = 0.
      ee = 0.
      ev = 0.
      do k = 1,3
        xyzmid(k) = 0.
        dxyz_at(k) = 0.
      end do
      do i = 1,n1 ! for each protein target atom
        e_at_max = 0.
        e_at_tot = 0.
        ninter_at = 0
        indri = 3 + i
        indqi = 3 + n1 + i
        indci = 3 + 2*n1 + 3*(i-1)
        r1 = atom_data(indri)
        q1 = atom_data(indqi)
        do k = 1,3
          indci = indci + 1
          x1(k) = atom_data(indci)
        end do
        do j = 1,n2 ! for each ligand atom
          indrj = 3 + 5*n1 + j
          indqj = 3 + 5*n1 + n2 + j
          indcj = 3 + 5*n1 + 2*n2 + 3*m_best*n2 + 3*(j-1)
          r2 = atom_data(indrj)
          q2 = atom_data(indqj)
          d2 = 0.
          do k = 1,3 ! distance and interatomic vector
            indcj = indcj + 1
            x2(k) = atom_data(indcj)
            dxyz(k) = x1(k) - x2(k)
            d2 = d2 + dxyz(k)*dxyz(k)
          end do
          if(d2 .lt. rcut2)then ! within distance cutoff
            ninter_at = ninter_at + 1
            sigma = r1 + r2
            dd = sqrt(d2) + 1.e-3
            rr = sigma/dd
            rr2 = rr*rr
            rr6 = rr2*rr2*rr2
            rr9 = rr6*rr2*rr
            evdw = vfact*(rr9 - rr6)
            ev = ev + evdw
            eelect = efact*q1*q2/dd
            ee = ee + eelect
            e_at = eelect + evdw
            e_at_tot = e_at_tot + e_at
            if(abs(e_at) .ge. abs(e_at_max))then
              e_at_max = e_at
              do k = 1,3
                xyzmid(k) = x2(k) + 0.5*dxyz(k)
                dxyz_at(k) = dxyz(k)
              end do
            end if
          end if
        end do
        if(ninter_at > 0)then ! create midline circle for strongest pair interaction
          d2 = 0.
          do k = 1,3
            d2 = d2 + dxyz_at(k)*dxyz_at(k)
          end do
          dd = sqrt(d2)
          do k = 1,3
            dxyz_at(k) = dxyz_at(k)/dd
          end do
          e_abs = 0.5*abs(e_at_tot)
          if(e_at_tot .ge. 0.)then
            color = 0.7 + min(0.3,e_abs)
          else
            color = 0.3 - min(0.3,e_abs)
          end if
c          print *,'color: ',color,e_abs,e_at_tot
          radius = 1. - dd/RCUT
c
c store circle data
c
          ninter = ninter + 1
          if(ninter .gt. MAXCIRCLE)then
            print *,'ERROR: > max pairwise interactions: ',maxcircle
            stop
          end if
          do k = 1,3
            circle_mid(k,ninter) = xyzmid(k)
            circle_perp(k,ninter) = dxyz_at(k)
          end do
          circle_color(ninter) = color
          circle_rad(ninter) = radius
        end if
      end do
      et = ee + ev
      print *,'number of interactions: ',ninter
      ndata = 12*num_points*ninter + 5
      if(ndata .gt. maxobjdata)then
        print *,'warning: too many circles, 
     & increase maxobjdata: ',maxobjdata
        ninter = int(maxobjdata/12/num_points - 5)
        print *,'reducing # of interactions to ',ninter
      end if
c
c generate pymol display object of dockeye object
c
      ndata = 1
      energy_obj(ndata) = LINEWIDTH
      ndata = ndata + 1
      energy_obj(ndata) = 4. ! linewidth
      ndata = ndata + 1
      energy_obj(ndata) = BEGIN
      ndata = ndata + 1
      energy_obj(ndata) = LINE
      do i = 1,ninter
        do k = 1,3
          vp(k) = circle_perp(k,i)
        end do
        call vperp(vp,v1)
        call vnorm(v1)
        call vcross(v1,vp,v2)
        do k = 1,3
          vbeg(k) = circle_rad(i) + circle_mid(k,i)
        end do
        do j = 1,num_points
          angle = j*tpi/num_points
          ca = cos(angle)
          sa = sin(angle)
          do k = 1,3
            vend(k) = circle_rad(i)*(ca*v1(k) + sa*v2(k))
     &      + circle_mid(k,i)
          end do
          call color_map(circle_color(i),rgb)
          ndata = ndata + 1
          energy_obj(ndata) = KOLOR
          do k = 1,3
            ndata = ndata + 1
            energy_obj(ndata) = rgb(k)
          end do
c
          ndata = ndata + 1
          energy_obj(ndata) = VERTEX
          do k = 1,3
            ndata = ndata + 1
            energy_obj(ndata) = vbeg(k)
          end do
c
          ndata = ndata + 1
          energy_obj(ndata) = VERTEX
          do k = 1,3
            ndata = ndata + 1
            energy_obj(ndata) = vend(k)
          end do
        end do
      end do
      ndata = ndata + 1
      energy_obj(ndata) = END
c
      ndata = ndata + 1
      energy_obj(ndata) = et
      ndata = ndata + 1
      energy_obj(ndata) = ee
      ndata = ndata + 1
      energy_obj(ndata) = ev
      ndata = ndata + 1
      energy_obj(ndata) = m_best
      energy_obj(1) = ndata
      print *,'# of return data: ',ndata
      print *,energy_obj(1),energy_obj(2),energy_obj(3),energy_obj(4)
      print *,energy_obj(ndata),energy_obj(ndata-1),energy_obj(ndata-2),
     &     energy_obj(ndata-3),energy_obj(ndata-4)
      return
      end
c
c==========================
      function vdot(v1,v2)
      implicit none
      real*8 v1(3),v2(3),vdot
      vdot = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      return
      end
c
      subroutine vcross(v1,v2,v3)
      implicit none
      real*8 v1(3),v2(3),v3(3)
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      return
      end
c
      subroutine vnorm(v1)
      implicit none
      real*8 v1(3),rdot,vdot
      rdot = vdot(v1,v1)
      rdot = sqrt(rdot) + 1.e-6
      v1(1) = v1(1)/rdot
      v1(2) = v1(2)/rdot
      v1(3) = v1(3)/rdot
      return
      end
c
      subroutine vperp(v1,v3)
      implicit none
      real*8 v1(3),v2(3),v3(3)
      v2(1) = v1(2)
      v2(2) = v1(3)
      v2(3) = v1(1)
      call vcross(v1,v2,v3)
      return
      end
c
      subroutine color_map(color,rgb)
      implicit none
      real*8 color,rgb(3)
      rgb(1) = 1.
      rgb(2) = 1.
      rgb(3) = 1.
      if((color .gt. 0.) .and. (color .lt. 0.5))then
        rgb(1) = 0.
        rgb(2) = 2.*color
        rgb(3) = 1. - 2.*color
      else
        if((color .ge. 0.5) .and. (color .le. 1.0))then
          rgb(1) = 2.*color - 1.
          rgb(2) = 2. - 2.*color
          rgb(3) = 0.
        end if
      end if
      return
      end
