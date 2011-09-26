#!/usr/bin/env python

import sys
import numpy as np

from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
class HCPOFactory(SimpleOrthorhombicFactory):
    "A factory for creating HCPO lattices."
    xtal_name = "hcpo"
    bravais_basis = [[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.5, 1.0/6.0, 0.5],
                     [0.0, 2.0/3.0, 0.5]]
    
HCPO = HCPOFactory()

argc = len(sys.argv)
if argc < 4:
    print 'usage:', sys.argv[0], 'Al fcc lp Si EcohSi'
    sys.exit(1)

el1 = sys.argv[1]
str = sys.argv[2]
lp1 = float(sys.argv[3])
el2 = sys.argv[4]
esub2 = float(sys.argv[5])
if str == 'hcp':
    if argc < 7:
        print 'usage:', sys.argv[0], 'Mg hcp lp Si EcohSi catoi'
        sys.exit(1)
    else:
        catoi = float(sys.argv[6])

species = [el1, el2]
import model
from model import pick_elements
pick_elements(model, species)

from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=model.parameters, files=model.files, specorder=species)

print "bulk: ", el1, str, lp1
if str == "hcp":
    print "catoi:", catoi
print "subst: ", el2, esub2

if str == "fcc" :
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic(symbol=el1,latticeconstant=lp1,size=(5,5,5))
elif str == "bcc" :
    from ase.lattice.cubic import BodyCenteredCubic
    atoms = BodyCenteredCubic(symbol=el1,latticeconstant=lp1,size=(5,5,5))
elif str == "dia" :
    from ase.lattice.cubic import Diamond
    atoms = Diamond(symbol=el1,latticeconstant=lp1,size=(3,3,3))
elif str == "hcp" :
    from math import sqrt
    a=lp1
    b=sqrt(3.)*lp1
    c=lp1*sqrt(8./3.)*catoi
    # this did not work
    # maybe wrong setting of LAMMPS volume optimization for triclinic cell
    # from ase.lattice.hexagonal import HexagonalClosedPacked
    # atoms = HexagonalClosedPacked(symbol=el1,latticeconstant=(lp1,c),size=(5,5,5))
    atoms = HCPO(symbol=el1,latticeconstant=(a,b,c),size=(8,4,4))

else :
    raise RuntimeError, "unknown structure \""+ str + "\""

atoms.set_calculator(calc)

v1=atoms.get_volume()
print "v1:", v1
n1=atoms.get_number_of_atoms()
print "n1:", n1
v1pa=v1/n1
print "v1pa:", v1pa

# initial
ene0 = atoms.get_potential_energy()
print "ene0:", ene0
print "ene0pa:", ene0/atoms.get_number_of_atoms()
ene0pa=ene0/atoms.get_number_of_atoms()

nm1 = n1-1
# esub in meamf is positive
subtr = nm1*ene0pa + 1.0*(-1.)*esub2

atoms[0].symbol = el2
print "subst n:", atoms[0].number

ene1nm = atoms.get_potential_energy()
print "ene1nm:", ene1nm
print "ene1nmf:", ene1nm - subtr

# minimize pos
model.parameters["minimize"] = "1.0e-25 1.0e-25 10000 10000"
#model.parameters["minimize"] = "1.0e-5 1.0e-5 10000 10000"
ene2m = atoms.get_potential_energy()
print "ene2m:", ene2m
print "ene2mf:", ene2m - subtr

atoms3 = calc.atoms

# minimize pos and cell
if str == "hcp":
    model.parameters["fix"] = "1 all box/relax aniso 0 couple xy vmax 0.000001"
else:
    model.parameters["fix"] = "1 all box/relax aniso 0 couple xyz vmax 0.000001"

atoms3.set_calculator(calc)
ene3v = atoms3.get_potential_energy()
print "ene3v:", ene3v
print "ene3vf:", ene3v - subtr
ene3vf = ene3v - subtr

esubf=ene3vf
print "esubf:", esubf

atoms4 = calc.atoms
v4 = atoms4.get_volume()
print "V4:", v4

vintf = v4-nm1*v1pa
print "vintf: ", vintf
print "vintfr: ", vintf/v1pa

#from ase.visualize import view
#view(atoms4)
