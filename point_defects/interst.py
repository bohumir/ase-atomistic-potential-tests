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
    print 'usage:', sys.argv[0], 'Al fcc octa/tetr/dump'
    sys.exit(1)

el1=sys.argv[1]
str=sys.argv[2]
dfc=sys.argv[3]

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

if str != 'hcp':
    lp,esub = get_meam_lp_esub(infile=meamf)
else:
    lp = 3.2027793
    catoi = 0.991824332358
    esub = -1.51011430257

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

print "meamfvals:", el1, lp, esub
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

from ase.atom import Atom
# insert
if dfc == "octa":
    if str == "fcc":
        atoms.append(Atom(el1, position=(0,0,0.5*lp)))
    elif str == "hcp":
        atoms.append(Atom(el1, position=(0,lp*2./3.*sqrt(3.)/2.,lp*1./2.*0.5*sqrt(8./3.)*catoi)))
    elif str == "bcc":
        atoms.append(Atom(el1, position=(1./2.*lp,0.,1./2.*lp)))
elif dfc == "tetr":
    if str == "fcc":
        atoms.append(Atom(el1, position=(1./4.*lp,1./4.*lp,1./4.*lp)))
    elif str == "hcp":
        atoms.append(Atom(el1, position=(lp*0.5,lp*1./3.*sqrt(3.)/2.,lp*1./3.*0.5*sqrt(8./3.)*catoi)))
    elif str == "bcc":
        atoms.append(Atom(el1, position=(1./4.*lp,0.,1./2.*lp)))
elif dfc == "dumb":
    if str == "fcc":
        atoms[0].position = [0, 0, 0.2*lp]
        atoms.append(Atom(el1, position=(0,0,-0.2*lp)))
    elif str == "bcc":
        atoms[0].position = [0.25*lp, 0.25*lp, 0.]
        atoms.append(Atom(el1, position=(-0.25*lp, -0.25*lp, 0.)))
    elif str == "hcp":
        atoms[0].position = [0., 0., 0.4*lp]
        atoms.append(Atom(el1, position=(0., 0., -0.4*lp)))
    elif str == "dia":
        atoms[0].position = [-0.1*lp, 0.1*lp, 0]
        atoms.append(Atom(el1, position=(0.1*lp, -0.1*lp, 0)))
elif dfc == "dum111":
    if str == "bcc":
        off=0.175
        atoms[0].position = [off*lp, off*lp, off*lp]
        atoms.append(Atom(el1, position=(-off*lp, -off*lp, -off*lp)))
elif dfc == "dum100":
    if str == "bcc":
        off=0.3
        atoms[0].position = [off*lp, 0., 0.]
        atoms.append(Atom(el1, position=(-off*lp, 0., 0.)))
else:
    raise RuntimeError, "unknown structure/defect pair: \""+ str + "/" + dfc + "\""

np1=n1+1
ene1nm = atoms.get_potential_energy()
print "ene1nm:", ene1nm
print "ene1nmpa:", ene1nm/np1

#from ase.visualize import view
#view(atoms)

# minimize pos
parameters["minimize"] = "1.0e-10 1.0e-10 10000 10000"
#parameters["minimize"] = "1.0e-5 1.0e-5 10000 10000"
ene2m = atoms.get_potential_energy()
print "ene2m:", ene2m
print "ene2mpa:", ene2m/np1

atoms3 = calc.atoms

# minimize pos and cell
if str == "hcp":
    parameters["fix"] = "1 all box/relax aniso 0 couple xy vmax 0.000001"
else:
    parameters["fix"] = "1 all box/relax aniso 0 couple xyz vmax 0.000001"

atoms3.set_calculator(calc)
ene3v = atoms3.get_potential_energy()
print "ene3v:", ene3v
print "ene3vpa:", ene3v/np1

eintf = ene3v-np1*ene0pa
print "eintf:", eintf

atoms4 = calc.atoms
v4 = atoms4.get_volume()
print "V4:", v4

vintf = v4-np1*v1pa
print "vintf: ", vintf
print "vintfr: ", vintf/v1pa

#from ase.visualize import view
#view(atoms4)
