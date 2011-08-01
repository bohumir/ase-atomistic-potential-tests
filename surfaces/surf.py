#!/usr/bin/env python

import sys
import numpy as np

argc = len(sys.argv)
if argc < 8:
    print 'usage:', sys.argv[0], 'Al fcc 111 rep1 rep2 rep3 vac'
    sys.exit(1)

el1=sys.argv[1]
struct=sys.argv[2]
surf=sys.argv[3]
rep1=float(sys.argv[4])
rep2=float(sys.argv[5])
rep3=float(sys.argv[6])
mysize=(rep1,rep2,rep3)
vac=float(sys.argv[7])


def get_meam_lp_esub(infile="meamf"):
    from ase.parallel import paropen
    f = paropen(infile, 'r')
    while True:
        line = f.readline()
        if not line:
            break
        if el1+"S" in line:
            line = f.readline()
            # lp is sixth on the second line, esub seventh 
            lp = float(line.split()[5])
            esub = float(line.split()[6])
    f.close()
    return lp, esub

pair_style = "meam"
meamf = "meamf"
meamp = "meam.alsimgcufe"
subspec = [ el1 ]
pair_coeff = [ "* * meamf AlS SiS MgS CuS FeS "\
                   + meamp + " " + subspec[0] + "S " ]
parameters = { "pair_style" : pair_style, "pair_coeff" : pair_coeff }
files = [ meamf, meamp ]

from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=parameters, files=files, specorder=subspec)

if el1 != 'Mg':
    lp,esub = get_meam_lp_esub(infile=meamf)
else:
    lp = 3.2027793
    catoi = 0.991824332358
    esub = 1.51011430257

from ase.lattice.surface import *

print "el1: ", el1
print "struct: ", struct
print "surf: ", surf
print "size: ", mysize
print "lp: ", lp
print "vac: ", vac
if surf == "10m10":
    myorth = True
else:
    myorth = False

if el1 != 'Mg':
    atoms=surface(el1, struct, surf, size=mysize, a=lp, c=None, vacuum=vac, orthogonal=myorth)
else:
    c=lp*sqrt(8./3)*catoi
    atoms=surface(el1, struct, surf, size=mysize, a=lp, c=c, vacuum=vac, orthogonal=myorth)
    

atoms.set_calculator(calc)
# eV/A^2 to mJ/m^2
un=16021.765

# unrelaxed
#
print "orig"
n1=atoms.get_number_of_atoms()
ene1 = atoms.get_potential_energy()
ene1pa=ene1/n1
print "n1:", n1
print "ene1:", ene1
print "ene1pa:", ene1pa
print "ecoh_meam:", -esub

cell1=atoms.get_cell()
sarea1=np.linalg.norm(np.cross(cell1[0],cell1[1]))
esurff1=(ene1-n1*(-esub))/(2*sarea1)*un
print "sarea1: ", sarea1
print "esurff1: ", esurff1, " mJ/m^2"

# minimize pos
#
print "relaxed"
parameters["minimize"] = "1.0e-25 1.0e-25 10000 10000"
ene2m = atoms.get_potential_energy()
ene2mpa = ene2m/n1
print "ene2m:", ene2m
print "ene2mpa:", ene2mpa

# these are after step 2
atoms2 = calc.atoms
# this will get the wrong cell for bcc 110 and 111 (a1 twice loger),
# (looks like a lammps.py bug)
# cell2 = atoms2.get_cell()
# but we are not minimizing cell - atoms have original (and more exact) cell
cell2 = atoms.get_cell()

sarea2 = np.linalg.norm(np.cross(cell2[0],cell2[1]))
esurff2m = (ene2m-n1*(-esub))/(2*sarea2)*un
print "sarea2: ", sarea2
print "esurff2m: ", esurff2m, " mJ/m^2"

#from ase.visualize import view
#view(atoms)
