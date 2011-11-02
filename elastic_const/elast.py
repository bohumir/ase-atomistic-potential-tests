#!/usr/bin/env python

from ase.units import kJ, _e

# obtain species, structure, and lattice parameter from command line
#
import sys
argc = len(sys.argv)
if argc < 5:
    print 'usage:', sys.argv[0], 'Al Mg nacl lp'
    sys.exit(1)

el1 = sys.argv[1]
el2 = sys.argv[2]
struct = sys.argv[3]
lp = float(sys.argv[4])
print "el1:", el1, "el2:", el2, "str:", struct, "lp:", lp

# read model specification from ./model.py file
# pick elements from the model
#
import model
from model import pick_elements
species = [el1, el2]
pick_elements(model, species)

# initialize LAMMPS calculator
#
from ase.calculators.lammps import LAMMPS
calc = LAMMPS(parameters=model.parameters, files=model.files,
              specorder=species)

# setup structure
#
import numpy as np
if struct in ["fcc", "nacl", "cu2mg", "mgcu2",
              "zns", "caf2", "f2ca", "alfe3", "fe3al"]:
    # reference cell
    refcell = np.array([[0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5],
                        [0.5, 0.5, 0.0]])
    if struct == "fcc":
        elems = [el1]
        poss = [(0, 0, 0)]
    elif struct == "nacl":
        elems = [el1, el2]
        poss = [(0, 0, 0),
                (0.5, 0.5, 0.5)]
    elif struct == "cu2mg":
        elems = [el1, el1, el1, el1, el2, el2]
        poss = [(0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875)]
    elif struct == "zns":
        elems = [el1, el2]
        poss = [(0.0, 0.0, 0.0),
                (0.25, 0.25, 0.25)]
    elif struct == "caf2":
        elems = [el1, el2, el2]
        poss = [(0.0, 0.0, 0.0),
                (0.25, 0.25, 0.25),
                (-0.25, -0.25, -0.25)]
    elif struct == "f2ca":
        elems = [el1, el1, el2]
        poss = [(0.25, 0.25, 0.25),
                (-0.25, -0.25, -0.25),
                (0.0, 0.0, 0.0)]
    elif struct == "alfe3":
        elems = [el1, el2, el2, el2]
        poss = [(0.0, 0.0, 0.0),
                (-0.5, 0.5, 0.5),
                (-0.25, -0.25, -0.25),
                (0.25, 0.25, 0.25)]
    elif struct == "fe3al":
        elems = [el1, el1, el1, el2]
        poss = [(-0.5, 0.5, 0.5),
                (-0.25, -0.25, -0.25),
                (0.25, 0.25, 0.25),
                (0.0, 0.0, 0.0)]
    elif struct == "cu2mg":
        elems = [el1, el1, el1, el1, el2, el2]
        poss = [(0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875)]
    elif struct == "mgcu2":
        elems = [el1, el1, el2, el2, el2, el2]
        poss = [(0.125, 0.125, 0.125),
                (0.875, 0.875, 0.875),
                (0.5, 0.5, 0.5),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0)]

elif struct in ["sc", "aucu3", "cu3au", "cscl", "cr3si", "sicr3"]:
    # reference cell
    refcell = np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])
    if struct == "sc":
        elems = [el1]
        poss = [(0.0, 0.0, 0.0)]
    elif struct == "aucu3":
        elems = [el1, el2, el2, el2]
        poss = [(0.0, 0.0, 0.0),
                (0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0)]
    elif struct == "cu3au":
        elems = [el1, el1, el1, el2]
        poss = [(0.0, 0.5, 0.5),
                (0.5, 0.0, 0.5),
                (0.5, 0.5, 0.0),
                (0.0, 0.0, 0.0)]
    elif struct == "cscl":
        elems = [el1, el2]
        poss = [(0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5)]
    elif struct == "sicr3":
        elems = [el1, el1, el2, el2, el2, el2, el2, el2]
        poss = [(0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5),
                (0.25, 0.50, 0.00),
                (0.75, 0.50, 0.00),
                (0.00, 0.25, 0.50),
                (0.00, 0.75, 0.50),
                (0.50, 0.00, 0.25),
                (0.50, 0.00, 0.75)]
    elif struct == "cr3si":
        elems = [el1, el1, el1, el1, el1, el1, el2, el2]
        poss = [(0.25, 0.50, 0.00),
                (0.75, 0.50, 0.00),
                (0.00, 0.25, 0.50),
                (0.00, 0.75, 0.50),
                (0.50, 0.00, 0.25),
                (0.50, 0.00, 0.75),
                (0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5)]

# create structure from elements, positions, and initial cell
#
from ase.lattice.spacegroup import crystal
init_cell = lp * refcell
atoms = crystal(elems, poss, cell=init_cell)

# assign calculator, get energy and volume per atom
#
atoms.set_calculator(calc)
epa0 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
vpa0 = atoms.get_volume() / atoms.get_number_of_atoms()
print "epa0:", epa0
print "vpa0:", vpa0, "\n"

# reoptimize/check volume
#
volumes = []
energies = []
for x in np.linspace(0.98, 1.02, 5):
    atoms.set_cell(init_cell * x, scale_atoms=True)
    volumes.append(atoms.get_volume() / atoms.get_number_of_atoms())
    energies.append(atoms.get_potential_energy() / atoms.get_number_of_atoms())
print "per atom volumes:", volumes
print "per atom energies:", energies

# fit EOS
#
from ase.utils.eos import EquationOfState
eos = EquationOfState(volumes, energies)
vpaf, epaf, B1 = eos.fit()
print "vpaf:", vpaf, "A^3"
print "epaf:", epaf, "eV"
print "B1:", B1 / kJ * 1.0e24, "GPa"

# get optimal lattice parameter from optimal volume
#
volrc = abs(np.linalg.det(refcell))
optlp = pow(vpaf * atoms.get_number_of_atoms() / volrc, 1. / 3.)
print "optlp:", optlp, "\n"

# get actual energy at optimal volume
#
opt_cell = optlp * refcell
atoms.set_cell(opt_cell, scale_atoms=True)
epao = atoms.get_potential_energy() / atoms.get_number_of_atoms()

strain = 0.001
diag = ([1, 0, 0],
        [0, 1, 0],
        [0, 0, 1])

# c44
#
defm1 = np.array([[0, strain, strain],
                  [strain, 0, strain],
                  [strain, strain, 0]])
cell1 = np.dot(opt_cell, defm1 + diag)
atoms.set_cell(cell1, scale_atoms=True)

ene1 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
dele = (ene1 - epao) / vpaf * _e / 1.0e-30  # eV/angstrom^3
c44 = dele / 6 / strain / strain / 1e9  # GPa
print "epao:", epao
print "ene1:", ene1
print "del:", ene1 - epao
print "c44:", c44, "GPa\n"

# c44 other way
#
defm2 = np.array([[0, strain, 0],
                  [strain, 0, 0],
                  [0, 0, 0]])
cell2 = np.dot(opt_cell, defm2 + diag)
atoms.set_cell(cell2, scale_atoms=True)

ene2 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
dele = (ene2 - epao) / vpaf * _e / 1.0e-30  # eV/angstrom^3
c44o = dele / 2 / strain / strain / 1e9  # GPa
print "epao:", epao
print "ene2:", ene2
print "del:", ene2 - epao
print "c44o:", c44o, "GPa\n"

# (c11-c12)/2
#
defm3 = np.array([[strain, 0, 0],
                  [0, 1 / (1 + strain) - 1, 0],
                  [0, 0, 0]])
cell3 = np.dot(opt_cell, defm3 + diag)
atoms.set_cell(cell3, scale_atoms=True)

ene3 = atoms.get_potential_energy() / atoms.get_number_of_atoms()
dele = (ene3 - epao) / vpaf * _e / 1.0e-30  # eV/angstrom^3
gammap = dele / 2 / strain / strain / 1e9  # GPa
print "epao:", epao
print "ene3:", ene3
print "del:", ene3 - epao
print "(c11-c12)/2:", gammap, "GPa\n"
