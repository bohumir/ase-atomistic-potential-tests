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
    print 'usage:', sys.argv[0], 'Al fcc Si'
    sys.exit(1)

el1=sys.argv[1]
str=sys.argv[2]
el2=sys.argv[3]

def get_meam_lp_esub(elx, infile="meamf"):
    from ase.parallel import paropen
    f = paropen(infile, 'r')
    while True:
        line = f.readline()
        if not line:
            break
        if elx+"S" in line:
            line = f.readline()
            # lp is sixth on the second line, esub seventh 
            lp = float(line.split()[5])
            esub = float(line.split()[6])
    f.close()
    return lp, esub

pair_style = "meam"
meamf = "meamf"
meamp = "meam.alsimgcufe"
subspec = [ el1, el2 ]
pair_coeff = [ "* * meamf AlS SiS MgS CuS FeS "\
                   + meamp + " " + subspec[0] + "S " + subspec [1] + "S"]
parameters = { "pair_style" : pair_style, "pair_coeff" : pair_coeff }
files = [ meamf, meamp ]

from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=parameters, files=files, specorder=subspec)

lpmg = 3.2027793
catoimg = 0.991824332358
esubmg = 1.51011430257

if el1 != 'Mg':
    lp1,esub1 = get_meam_lp_esub(el1,infile=meamf)
else:
    lp1 = lpmg
    catoi = catoimg
    esub1 = esubmg

if el2 != 'Mg':
    lp2,esub2 = get_meam_lp_esub(el2,infile=meamf)
else:
    lp2 = lpmg
    esub2 = esubmg

print "meamfvals1:", el1, lp1, esub1
print "meamfvals2:", el2, lp2, esub2

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
parameters["minimize"] = "1.0e-25 1.0e-25 10000 10000"
#parameters["minimize"] = "1.0e-5 1.0e-5 10000 10000"
ene2m = atoms.get_potential_energy()
print "ene2m:", ene2m
print "ene2mf:", ene2m - subtr

atoms3 = calc.atoms

# minimize pos and cell
if str == "hcp":
    parameters["fix"] = "1 all box/relax aniso 0 couple xy vmax 0.000001"
else:
    parameters["fix"] = "1 all box/relax aniso 0 couple xyz vmax 0.000001"

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
