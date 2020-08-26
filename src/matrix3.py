"""implement matrix rotations as
a class"""
import math
PI = 3.1415927
DTOR = PI/180.
def arctan(x,y):
    if(x == 0.):
      if(y > 0.):
        angle = 90.
      else:
        angle = -90.
    else:
      angle = math.atan(abs(y/x))/DTOR
      if(x < 0.): angle = 180. - angle
      if(y < 0.): angle = -1.*angle
    return angle

class Rotmat:
    def __init__(self):
        self.terms = [[1., 0., 0.],[0.,1.,0.],[0.,0.,1.]]
    def printm(self,s1='\n'):
        print(s1)
        print( "| %9.4f %9.4f %9.4f |" % (self.terms[0][0], self.terms[0][1], self.terms[0][2]))
        print( "| %9.4f %9.4f %9.4f |" % (self.terms[1][0], self.terms[1][1], self.terms[1][2]))
        print( "| %9.4f %9.4f %9.4f |" % (self.terms[2][0], self.terms[2][1], self.terms[2][2]))
    def term(self,i,j,value):
        self.terms[i][j] = value
    def transpose(self):
        val = self.terms[0][1]
        self.terms[0][1] = self.terms[1][0]
        self.terms[1][0] = val
        #
        val = self.terms[1][2]
        self.terms[1][2] = self.terms[2][1]
        self.terms[2][1] = val
        #
        val = self.terms[0][2]
        self.terms[0][2] = self.terms[2][0]
        self.terms[2][0] = val
    def __add__(self,m2):
        terms_out = []
        for i  in range(3):
            row = []
            for j in range(3):
                row.append(self.terms[i][j] + m2.terms[i][j])
            terms_out.append(row)
        return terms_out
    def __mul__(self,m2):
        terms_out = []
        for i  in range(3):
            row = []
            for j in range(3):
                tmp = 0.
                for k in range(3):
                    tmp+= self.terms[i][k]*m2.terms[k][j]
                row.append(tmp)
            terms_out.append(row)
        return terms_out
    def rotx(self,angle):
        for i in range(3):
            for j in range(3):
                self.terms[i][j] = 0.
            self.terms[0][0] = 1.
            self.terms[1][1] = math.cos(angle*math.pi/180.)
            self.terms[2][2] = math.cos(angle*math.pi/180.)
            self.terms[1][2] = -math.sin(angle*math.pi/180.)
            self.terms[2][1] = math.sin(angle*math.pi/180.)
    def roty(self,angle):
        for i in range(3):
            for j in range(3):
                self.terms[i][j] = 0.
            self.terms[1][1] = 1.
            self.terms[2][2] = math.cos(angle*math.pi/180.)
            self.terms[0][0] = math.cos(angle*math.pi/180.)
            self.terms[2][0] = -math.sin(angle*math.pi/180.)
            self.terms[0][2] = math.sin(angle*math.pi/180.)
    def rotz(self,angle):
        for i in range(3):
            for j in range(3):
                self.terms[i][j] = 0.
            self.terms[2][2] = 1.
            self.terms[0][0] = math.cos(angle*math.pi/180.)
            self.terms[1][1] = math.cos(angle*math.pi/180.)
            self.terms[0][1] = -math.sin(angle*math.pi/180.)
            self.terms[1][0] = math.sin(angle*math.pi/180.)
    def rot_vec(self,vector):
        vector_rot = []
        for i in range(3):
            tmp = 0.
            for j in range(3):
                tmp +=self.terms[i][j]*vector[j]
            vector_rot.append(tmp)
        return vector_rot
    def polar_rot(self,phi,psi,chi):
        # TO rotation matrix from rotation by angle chi 
        # about axis in polar coords phi (angle to x-axis)
        # & psi (angle to x-y plane)
        cx = math.cos(DTOR*phi)
        cy = math.cos(DTOR*psi)
        cz = math.cos(DTOR*chi)
        sx = math.sin(DTOR*phi)
        sy = math.sin(DTOR*psi)
        sz = math.sin(DTOR*chi)
        self.terms[0][0] = cx**2*cy**2*(1.-cz) + cz
        self.terms[1][0] = cx*sx*cy**2*(1.-cz) - sz*sy
        self.terms[2][0] = cx*cy*sy*(1.-cz) + sz*sx*cy
        self.terms[0][1] = cx*sx*cy**2*(1.-cz) + sz*sy
        self.terms[1][1] = sx**2*cy**2*(1.-cz) + cz
        self.terms[2][1] = sx*cy*sy*(1.-cz) - sz*cx*cy
        self.terms[0][2] = cx*cy*sy*(1.-cz) - sz*sx*cy
        self.terms[1][2] = sx*cy*sy*(1.-cz) + sz*cx*cy
        self.terms[2][2] = sy**2*(1.-cz) + cz
    def rot_polar_old(self):
        # FROM rotation matrix to rotation by angle chi 
        # about axis in polar coords phi (angle to x-axis)
        # & psi (angle to x-y plane) 
        EPS = 1.e-6
        tr = self.terms[0][0]+self.terms[1][1]+self.terms[2][2]
        csk = (tr - 1.)/2.
        csk = min(csk,1.)
        csk = max(csk,-1.)
        chi = math.acos(csk)
        snk = math.sin(chi)
        chi = chi/DTOR
        if((1.-csk) < EPS):
          # print('kappa = 0 = no rotation')
          psi = 0.
          phi = 0.
          angles = [phi,psi,chi]
          return angles
        # endif
        if(1. - self.terms[2][2] < EPS):
          if(self.terms[1][0] > 0.):
            psi = -90.
          else:
            psi = 90.
          phi = 0.
          angles = [phi,psi,chi]
          return angles
        #end if
        if((1.+csk) < EPS):
          # write(6,*)'c kappa = 180'
          sns = math.sqrt((self.terms[2][2]+1.)/2.)
          sns = min(sns,1.)
          sns = max(sns,-1.)
          psi = math.asin(sns)
          css = math.cos(psi)
          psi = psi/DTOR
          if(abs(css) < EPS):
            # write(6,*)'c psi = 90 or 270'
            # doesnt matter, also phi angle is immaterial
            phi = 0.
          elif(abs(sns) < EPS):
            # write(6,*)'c psi = 0, 180 (doesnot matter which)'
            csp = math.sqrt((self.terms[0][0]+1.)/2.)
            if(abs(csp) < EPS):
              # write(6,*)'c phi = 90, 270 (doesnot matter which)'
              phi = 90.
            else:
              # write(6,*)'phi not 90, 270'
              snp = self.terms[1][0]/(2.*csp)
              phi = arctan(csp,snp)
            #end if 
          else:
            # write(6,*)'c psi not 0,90,180,270'
            csp = self.terms[2][0]
            snp = self.terms[2][1]
            phi = arctan(csp,snp)
          #end if
        else:
          #   write(6,*)'c 0 < kappa < 180' 
          d12 = (self.terms[0][1] - self.terms[1][0])/2.
          sns = d12/snk
          sns = min(sns,1.)
          sns = max(sns,-1.)
          psi = math.asin(sns)
          css = math.cos(psi)
          psi = psi/DTOR
          if(abs(css) < 100.*EPS):
            # write(6,*)'c psi = 90 or 270'
            phi = 0.
          else:
            #   write(6,*)'c psi not 90 or 270'
            d13 = (self.terms[2][0] - self.terms[0][2])/2.
            d23 = (self.terms[1][2] - self.terms[2][1])/2.
            phi = arctan(d23,d13)
          #end if
        #end if
        angles = [phi,psi,chi]
        return angles
    def rot_cart(self):
        # FROM rotation matrix to cartesian rotations alpha, beta and gamma
        EPS = 1.e-6
        sy = self.terms[0][2]
        sy = min(sy,1.)
        sy = max(sy,-1.)
        beta = math.asin(sy)
        cy = math.cos(beta)
        beta = beta/DTOR
        if(abs(cy) < EPS):
          # y = 90 or 270
          alpha = 0.
          cz = self.terms[1][1]
          sz = -1*self.terms[1][0]
          gamma = arctan(cz,sz)
        else:
          # y not 90 or 270
          cysx = self.terms[1][2]
          cycx = self.terms[2][2]
          alpha = arctan(cycx,cysx)
          czcy = self.terms[0][0]
          szcy = self.terms[1][1]
          gamma = arctan(czcy,szcy)
        #end if
        angles = [alpha,beta,gamma]
        return angles
    def rot_polar(self):
        # FROM rotation matrix to rotation by angle chi 
        # about axis in polar coords phi (angle to x-axis)
        # & psi (angle to x-y plane) 
        tr = self.terms[0][0]+self.terms[1][1]+self.terms[2][2]
        csk = (tr - 1.)/2.
        csk = min(csk,1.)
        csk = max(csk,-1.)
        chi = math.acos(csk)
        snk = math.sin(chi)
        chi = chi/DTOR
        rotx = (self.terms[1][2] - self.terms[2][1])
        roty = (self.terms[2][0] - self.terms[0][2])
        rotz = (self.terms[0][1] - self.terms[1][0])
        rotmag = math.sqrt(rotx*rotx + roty*roty + rotz*rotz)
        psi = math.asin(rotz/rotmag)/DTOR
        phi = arctan(rotx,roty)
        #print(' phi, psi, chi: ',phi,psi,chi)
        angles = [phi,psi,chi]
        return angles
#=======================================
def rot_vec(rmt,vec,inv=0):
  vec_rot = [0.,0.,0.]
  if(inv == 0):
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[i][j]*vec[j]
  else:
    for i in range(3):
      for j in range(3):
        vec_rot[i] += rmt[j][i]*vec[j]
  return vec_rot

def vdot(v1,v2):
  dot = 0.
  for k in range(3):
    dot += v1[k]*v2[k]
  return dot
#
def vcross(v1,v2):
  cross = [0.,0.,0.]
  cross[0] = v1[1]*v2[2] - v1[2]*v2[1]
  cross[1] = v1[2]*v2[0] - v1[0]*v2[2]
  cross[2] = v1[0]*v2[1] - v1[1]*v2[0]
  return cross
#
def vnorm(v1):
  dot = vdot(v1,v1)
  dot = math.sqrt(dot)
  if(dot > 0.):
    for k in range(3):
      v1[k] /= dot
  return dot
#
def vperp(v1):
  v2 = [v1[1],v1[2],v1[0]]
  v3 = vcross(v1,v2)
  return v3
#=======================================
