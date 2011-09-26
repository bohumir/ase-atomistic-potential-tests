#!/usr/bin/env python

import sys
import numpy as np

argc = len(sys.argv)
if argc < 10:
    print 'usage:', sys.argv[0], 'Al fcc lp EcohAl 111 rep1 rep2 rep3 vac'
    sys.exit(1)

el1 = sys.argv.pop(1)
struct = sys.argv.pop(1)
lp = float(sys.argv.pop(1))
if struct == 'hcp':
    if argc < 11:
        print >> sys.stderr, 'usage:',\
            sys.argv[0], 'Mg hcp lp catoi EcohMg 0001 rep1 rep2 rep3 vac'    
        sys.exit(1)
    else:
        catoi = float(sys.argv.pop(1))
esub = float(sys.argv.pop(1))
surf = sys.argv.pop(1)
rep1 = float(sys.argv.pop(1))
rep2 = float(sys.argv.pop(1))
rep3 = float(sys.argv.pop(1))
mysize = (rep1,rep2,rep3)
vac = float(sys.argv.pop(1))

print "el1: ", el1
print "struct: ", struct
print "lp: ", lp
if struct == 'hcp':
    print "catoi: ", catoi
print "ecoh: ", esub
print "surf: ", surf
print "size: ", mysize
print "vacuum: ", vac

import model
from model import pick_elements
species = [el1]
pick_elements(model, species)

from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=model.parameters, files=model.files)

from ase.lattice.surface import *

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
print "esub:", -esub

cell1=atoms.get_cell()
sarea1=np.linalg.norm(np.cross(cell1[0],cell1[1]))
esurff1=(ene1-n1*(-esub))/(2*sarea1)*un
print "sarea1: ", sarea1
print "esurff1: ", esurff1, " mJ/m^2"

# minimize pos
#
print "relaxed"
model.parameters["minimize"] = "1.0e-25 1.0e-25 10000 10000"
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
