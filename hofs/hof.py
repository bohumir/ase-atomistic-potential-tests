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
else:
    el1 = sys.argv[1]
    el2 = sys.argv[2]
    str = sys.argv[3]
    esub1 =  float(sys.argv[4])
    esub2 =  float(sys.argv[5])

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

est_nnd = 2.5
# assume cubic packing to estimate lp from nnd, down 0.9
volpa = 0.95*est_nnd*est_nnd*est_nnd

# per cell
volpc = abs(np.linalg.det(primc))
numatpc = len(elems)

lp = pow(volpa*numatpc/volpc, 1/3.)
print "initlp:", lp, "initvolpa:", volpa

lat0 = lp*primc
mys = crystal(elems,
              poss,
              cell=lat0)

mys.set_calculator(calc)
model.parameters["minimize"] = "1.0e-25 1.0e-25 100000 100000"
model.parameters["fix"] = "1 all box/relax aniso 0 couple xyz vmax  0.00001"

ene0 = mys.get_potential_energy()
atomsmin=calc.atoms
vol0 = atomsmin.get_volume()
print "ene0:", ene0
print "vol0:", vol0
print "vol0pa:", vol0/mys.get_number_of_atoms()

n1=elems.count(el1)
print "esub1:", esub1, "n1:", n1

if argc != 4:
    n2 = elems.count(el2)
    print "esub2:", esub2, "n2:", n2
    hof = (ene0+n1*esub1+n2*esub2)/(n1+n2)
else:
    hof = ene0/n1-esub1

print "hof:", hof*1000, "meV/atom"
