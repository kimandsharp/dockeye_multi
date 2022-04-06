"""
functions to generate header, tail of  pymol file"""
def pymol_cgo_head(objectname,ctype):
    filename = objectname+'.py'
    pymol_file = open(filename,"w")
    string = "cmd.delete('"+objectname+"')\n"
    pymol_file.write(string)
    pymol_file.write('from pymol.cgo import *\n')
    pymol_file.write('from pymol import cmd\n')
    pymol_file.write('import math\n')
    pymol_file.write('obj = [\n')
    pymol_file.write('LINEWIDTH, 5.0,\n')
    string = '    BEGIN,' + ctype +',\n'
    pymol_file.write(string)
    return pymol_file

# tail of pymol file
def pymol_cgo_tail(pymol_file,objectname):
    pymol_file.write('   END\n')
    pymol_file.write('   ]\n')
    string = "cmd.load_cgo(obj,'"+objectname+"')\n"
    pymol_file.write(string)
    pymol_file.close()

def pymol_cgo_line(pymol_file,vbeg,vend,color):
    # put a line
    red = 1.
    green = 1.
    blue = 1.
    if((color >= 0.) and (color < 0.5)):
      red = 0.
      green = 2.*color
      blue = (1. - 2.*color)
    elif((color >= 0.5) and (color <= 1.0)):
      red = 2.*color  - 1.
      green = 2. - 2.*color
      blue = 0.
    pymol_file.write('COLOR, {:3.2f} , {:3.2f} , {:3.2f} ,'.format(red, green, blue))
    cgo_string = ' VERTEX, {:8.3f} , {:8.3f}, {:8.3f} ,\n'.format(vbeg[0],vbeg[1],vbeg[2])
    pymol_file.write(cgo_string)
    cgo_string = ' VERTEX, {:8.3f} , {:8.3f}, {:8.3f} ,\n'.format(vend[0],vend[1],vend[2])
    pymol_file.write(cgo_string)
