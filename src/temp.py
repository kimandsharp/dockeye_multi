def crd_flip(crds,natom):
   """ for any ligand conformer, generate 3 flipped versions around x,y,z 
       principal moment of inertia axes """
   #
   # find center, and set to 0,0,0
   gcen = np.zeros(3)
   for i in range(natom):
     for k in range(3):
       gcen[k] += crds[i][k]
   for k in range(3):
     gcen[k] /= natom
   #print('center: ',gcen)
   for i in range(natom):
     for k in range(3):
       crds[i][k] -= gcen[k]
   #
   # find moi
   moi = np.zeros((3,3))
   for i in range(natom):
     for k in range(3):
       for j in range(3):
         moi[k][j] += crds[i][k]*crds[i][j]
   #print(moi)
   #
   # find principal axes (eigenvectors of moi)
   # which will form rotation matrix to rotate into xyz frame
   # then flip around x, y, and then z
   e,ev = la.eigh(moi)
   #print(e)
   #print(ev)
   crds_rot = np.zeros((natom,3,3))
   xyz_rot = np.zeros(3)
   xyz_flip = np.zeros((3,3))
   for i in range(natom):
     for k in range(3):
       #
       # rotate original with major axes aligned along xyz
       xyz_rot[k] = 0.
       for j in range(3):
         xyz_rot[k] += crds[i][j]*ev[j][k]
       #
       # flip around x, y, and z
     # x rotation
     xyz_flip[0][0] = +xyz_rot[0]
     xyz_flip[1][0] = -xyz_rot[1]
     xyz_flip[2][0] = -xyz_rot[2]
     # y rotation
     xyz_flip[0][1] = -xyz_rot[0]
     xyz_flip[1][1] = +xyz_rot[1]
     xyz_flip[2][1] = -xyz_rot[2]
     # z rotation
     xyz_flip[0][2] = -xyz_rot[0]
     xyz_flip[1][2] = -xyz_rot[1]
     xyz_flip[2][2] = +xyz_rot[2]
     for l in range(3):
       #
       # rotate and translate back
       for k in range(3):
         crds_rot[i][k][l] = gcen[k]
         for j in range(3):
           crds_rot[i][k][l] += xyz_flip[j][l]*ev[k][j]
   return crds_rot
