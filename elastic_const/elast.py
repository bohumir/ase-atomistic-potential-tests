#!/usr/bin/env python

import sys
import numpy as np
from ase.lattice.spacegroup import crystal

argc = len(sys.argv)
if argc < 5:
    print 'usage:', sys.argv[0], 'Al Mg nacl lp'
    sys.exit(1)

el1 = sys.argv[1]
el2 = sys.argv[2]
str = sys.argv[3]
lp = float(sys.argv[4])

print "el1:", el1, "el2", el2, "str", str, "lp:", lp

import model
from model import pick_elements
species = [el1, el2]
pick_elements(model, species)

from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=model.parameters, files=model.files, specorder=species)

if str == "fcc" or str == "nacl" or str == "cu2mg" or str == "mgcu2" or str == "zns" or str == "caf2" or str == "f2ca" or str == "alfe3" or str == "fe3al":
    primc = np.array([[0.0, 0.5, 0.5],
                      [0.5, 0.0, 0.5],
                      [0.5, 0.5, 0.0]])
    if str == "fcc":
        elems = [el1]
        poss = [(0, 0, 0)];
    elif str == "nacl":
        elems = [el1, el2]
        poss = [(0, 0, 0),
                (0.5, 0.5, 0.5)]
    elif str == "cu2mg":
        elems = [el1,el1,el1,el1,el2,el2]
        poss = [(0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875)]
    elif str == "zns":
        elems = [el1,el2]
        poss = [(0.0, 0.0, 0.0),
                (0.25, 0.25, 0.25)]
    elif str == "caf2":
        elems = [el1,el2,el2]
        poss = [(0.0, 0.0, 0.0),
                (0.25, 0.25, 0.25),
                (-0.25, -0.25, -0.25)]
    elif str == "f2ca":
        elems = [el1,el1,el2]
        poss = [(0.25, 0.25, 0.25),
                (-0.25, -0.25, -0.25),
                (0.0, 0.0, 0.0)]
    elif str == "alfe3":
        elems = [el1,el2,el2,el2]
        poss = [(0.0, 0.0, 0.0),
                (-0.5, 0.5, 0.5),
                (-0.25, -0.25, -0.25),
                (0.25, 0.25, 0.25)]
    elif str == "fe3al":
        elems = [el1,el1,el1,el2]
        poss = [(-0.5, 0.5, 0.5),
                (-0.25, -0.25, -0.25),
                (0.25, 0.25, 0.25),
                (0.0, 0.0, 0.0)]
    elif str == "cu2mg":
        elems = [el1,el1,el1,el1,el2,el2]
        poss = [(0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875)]
    elif str == "mgcu2":
        elems = [el1,el1,el2,el2,el2,el2]
        poss = [(0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875),
                (0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0)]

elif str == "sc" or str == "aucu3" or str == "cu3au" or str == "cscl" or str == "cr3si" or str == "sicr3":
    primc = np.array([[1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])
    if str == "sc":
        elems = [el1]
        poss = [(0,0,0)]
    elif str == "aucu3":
        elems = [el1,el2,el2,el2]
        poss = [(0.0, 0.0, 0.0),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0)]
    elif str == "cu3au":
        elems = [el1,el1,el1,el2]
        poss = [(0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.0, 0.0, 0.0)]
    elif str == "cscl":
        elems = [el1,el2]
        poss = [(0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5)]
    elif str == "sicr3":
        elems = [el1,el1,el2,el2,el2,el2,el2,el2]
        poss = [(0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5),
                (0.25, 0.50, 0.00),
                (0.75, 0.50, 0.00),
                (0.00, 0.25, 0.50),
                (0.00, 0.75, 0.50),
                (0.50, 0.00, 0.25),
                (0.50, 0.00, 0.75)]
    elif str == "cr3si":
        elems = [el1,el1,el1,el1,el1,el1,el2,el2]
        poss = [(0.25, 0.50, 0.00),
                (0.75, 0.50, 0.00),
                (0.00, 0.25, 0.50),
                (0.00, 0.75, 0.50),
                (0.50, 0.00, 0.25),
                (0.50, 0.00, 0.75),
                (0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5)]

# per cell
vpc = abs(np.linalg.det(primc))

lat0 = lp*primc
mys = crystal(elems,
              poss,
              cell=lat0)

mys.set_calculator(calc)
epa0 = mys.get_potential_energy()/mys.get_number_of_atoms()
vpa0 = mys.get_volume()/mys.get_number_of_atoms()
print "epa0:", epa0
print "vpa0:", vpa0, "\n"

from ase.units import kJ
units = 1.6021765e-19*1e30 # eV/angstrom^3
from ase.utils.eos import EquationOfState

volumes = []
energies = []
elstr = el1+el2+'-'+str
for x in np.linspace(0.98, 1.02, 5):
    mys.set_cell(lat0 * x, scale_atoms = True)
    volumes.append(mys.get_volume()/mys.get_number_of_atoms())
    energies.append(mys.get_potential_energy()/mys.get_number_of_atoms())
print "Volumespa:", volumes
print "Energiespa:", energies

eos = EquationOfState(volumes, energies)
vpao, epaoe, B1 = eos.fit()
lpopt1 = pow(vpao*mys.get_number_of_atoms()/vpc, 1/3.)

print 'vpao:', vpao, 'A^3'
print 'epaoe:', epaoe, 'eV'
print 'B1:', B1 / kJ * 1.0e24, 'GPa'
print 'lpopt1:', lpopt1
#eos.plot(elstr+'-eos.pdf')#,show = True)

latopt = lpopt1*primc
mys.set_cell(latopt, scale_atoms = True)
epaoa = mys.get_potential_energy()/mys.get_number_of_atoms()

perc = 0.001
unit = ([1,0,0],
        [0,1,0],
        [0,0,1]);

# c44
defmrh = np.array([[0, perc, perc], [perc, 0, perc], [perc, perc, 0]])
lat1 = np.dot(latopt,defmrh+unit)
mys.set_cell(lat1,scale_atoms = True)

ene1 = mys.get_potential_energy()/mys.get_number_of_atoms()
dele = (ene1-epaoa)/vpao*units
c44 = dele/6/perc/perc/1e9 # GPa
print "ene1:", ene1*mys.get_number_of_atoms()
print "epaoa:", epaoa
print "del:", ene1-epaoa
print "c44:", c44, "GPa\n"

# c44 Dynamo
defmdy = np.array([[0, perc, 0], [perc, 0, 0], [0, 0, 0]])
lat3 = np.dot(latopt,defmdy+unit)
mys.set_cell(lat3,scale_atoms = True)

ene3 = mys.get_potential_energy()/mys.get_number_of_atoms()
dele = (ene3-epaoa)/vpao*units
c44dy = dele/2/perc/perc/1e9 # GPa
print "ene3:", ene3*mys.get_number_of_atoms()
print "epaoa:", epaoa
print "del:", ene3-epaoa
print "c44dy:", c44dy, "GPa\n"

# (c11-c12)/2
defmsi = np.array([[perc, 0, 0], [0, 1/(1+perc)-1, 0], [0, 0, 0]])
lat2 = np.dot(latopt,defmsi+unit)
mys.set_cell(lat2,scale_atoms = True)

ene2 = mys.get_potential_energy()/mys.get_number_of_atoms()
dele = (ene2-epaoa)/vpao*units
gammap = dele/2/perc/perc/1e9 # GPa
tmp =  ene2*mys.get_number_of_atoms()
print 'ene2: %0.9f' % tmp
print "epaoa:", epaoa
print "del:", ene2-epaoa
print "(c11-c12)/2:", gammap, "GPa\n"
