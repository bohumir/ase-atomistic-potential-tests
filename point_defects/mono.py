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
if argc < 3:
    print 'usage:', sys.argv[0], 'Al fcc lp [catoi]'
    sys.exit(1)

el1=sys.argv[1]
str=sys.argv[2]
lp = float(sys.argv[3])
if str == 'hcp':
    if argc < 4:
        print 'usage:', sys.argv[0], 'Mg hcp octa/tetr/dump lp catoi'    
    else:
        catoi = float(sys.argv[4])

from ase.calculators.lammps import LAMMPS

import model
calc = LAMMPS(parameters=model.parameters, files=model.files)
model.parameters["pair_coeff"][0]  += " " + el1 + model.ext

if str == "fcc" :
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic(symbol=el1,latticeconstant=lp,size=(5,5,5))
elif str == "bcc" :
    from ase.lattice.cubic import BodyCenteredCubic
    atoms = BodyCenteredCubic(symbol=el1,latticeconstant=lp,size=(5,5,5))
elif str == "dia" :
    from ase.lattice.cubic import Diamond
    atoms = Diamond(symbol=el1,latticeconstant=lp,size=(3,3,3))
elif str == "hcp" :
    from math import sqrt
    a=lp
    b=sqrt(3.)*lp
    c=lp*sqrt(8./3.)*catoi
    # this did not work
    # maybe wrong setting of LAMMPS volume optimization for triclinic cell
    # from ase.lattice.hexagonal import HexagonalClosedPacked
    # atoms = HexagonalClosedPacked(symbol=el1,latticeconstant=(lp,c),size=(5,5,5))
    atoms = HCPO(symbol=el1,latticeconstant=(a,b,c),size=(8,4,4))

else :
    raise RuntimeError, "unknown structure \""+ str + "\""

atoms.set_calculator(calc)

print "meamfvals:", el1, lp
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

# vacancy
atoms.pop(0)
nm1=n1-1
ene1nm = atoms.get_potential_energy()
print "ene1nm:", ene1nm
print "ene1nmpa:", ene1nm/nm1

# minimize pos
model.parameters["minimize"] = "1.0e-25 1.0e-25 10000 10000"
#model.parameters["minimize"] = "1.0e-5 1.0e-5 10000 10000"
ene2m = atoms.get_potential_energy()
print "ene2m:", ene2m
print "ene2mpa:", ene2m/nm1

atoms3 = calc.atoms

# minimize pos and cell
if str == "hcp":
    model.parameters["fix"] = "1 all box/relax aniso 0 couple xy vmax 0.000001"
else:
    model.parameters["fix"] = "1 all box/relax aniso 0 couple xyz vmax 0.000001"

atoms3.set_calculator(calc)
ene3v = atoms3.get_potential_energy()
print "ene3v:", ene3v
print "ene3vpa:", ene3v/nm1

evacf = ene3v-nm1*ene0pa
print "evacf:", evacf

atoms4 = calc.atoms
v4 = atoms4.get_volume()
print "V4:", v4

vvacf = v4-nm1*v1pa
print "vvacf: ", vvacf
print "vvacfr: ", vvacf/v1pa

#from ase.visualize import view
#view(atoms)
