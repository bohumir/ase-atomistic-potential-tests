#!/usr/bin/env python

import sys
import numpy as np
from ase.lattice.spacegroup import crystal

argc = len(sys.argv)
if argc < 6 and argc != 4:
    print 'usage:', sys.argv[0], 'Al Si nacl 3.353 4.63'
    print 'usage:', sys.argv[0], 'Al bcc 3.353'
    sys.exit(1)

if argc == 4:
    el1 = sys.argv[1]
    str = sys.argv[2]
    esub1 = sys.argv[3]
    print "el1:", el1, "esub1:", esub1
    binary = 0
else:
    el1 = sys.argv[1]
    el2 = sys.argv[2]
    str = sys.argv[3]
    esub1 =  float(sys.argv[4])
    esub2 =  float(sys.argv[5])
    print "el1:", el1, "el2:", el2, "esub1:", esub1, "esub2:", esub2
    binary = 1

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
        elems = [el1,el2]
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

n1=elems.count(el1)
print "n1:", n1
if binary == 1:
    n2 = elems.count(el2)
    print "n2:", n2

est_nnd = 2.5
# assume cubic packing to estimate lp from nnd, down 0.9
vpa = 0.95*est_nnd*est_nnd*est_nnd

# per cell
vpc = abs(np.linalg.det(primc))
numatpc = len(elems)

lp = pow(vpa*numatpc/vpc, 1/3.)
print "initlp:", lp
print "initvpa:", vpa

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
from ase.utils.eos import EquationOfState

lestim = lat0
volumes = []
energies = []
elstr=el1+el2+'-'+str
for x in np.linspace(0.98, 1.02, 5):
    mys.set_cell(lestim * x, scale_atoms=True)
    volumes.append(mys.get_volume()/mys.get_number_of_atoms())
    energies.append(mys.get_potential_energy()/mys.get_number_of_atoms())
print "Volumespa:", volumes
print "Energiespa:", energies

eos = EquationOfState(volumes, energies)
vpa1, epa1, B1 = eos.fit()
lpopt1 = pow(vpa1*(n1+n2)/vpc, 1/3.)

print 'vpa1:', vpa1, 'A^3'
print 'epa1:', epa1, 'eV'
print 'B1:', B1 / kJ * 1.0e24, 'GPa'
print 'lpopt1:', lpopt1
#eos.plot(elstr+'-eos.pdf')#,show=True)

if binary == 1:
    hof = epa1+esub1*n1/(n1+n2)+esub2*n2/(n1+n2)
else:
    hof = epa1/n1-esub1
print "hof1:", hof*1000, "meV/atom\n"

# optimize again
lestim = lpopt1*primc
volumes = []
energies = []
elstr=el1+el2+'-'+str
for x in np.linspace(0.98, 1.02, 5):
    mys.set_cell(lestim * x, scale_atoms=True)
    volumes.append(mys.get_volume()/mys.get_number_of_atoms())
    energies.append(mys.get_potential_energy()/mys.get_number_of_atoms())
print "Volumespa:", volumes
print "Energiespa:", energies

eos = EquationOfState(volumes, energies)
vpa1, epa1, B1 = eos.fit()
lpopt1 = pow(vpa1*(n1+n2)/vpc, 1/3.)

print 'vpa1:', vpa1, 'A^3'
print 'epa1:', epa1, 'eV'
print 'B1:', B1 / kJ * 1.0e24, 'GPa'
print 'lpopt1:', lpopt1

if binary == 1:
    hof = epa1+esub1*n1/(n1+n2)+esub2*n2/(n1+n2)
else:
    hof = epa1/n1-esub1
print "hof1:", hof*1000, "meV/atom"

#eos.plot(elstr+'-eos.pdf')#,show=True)
