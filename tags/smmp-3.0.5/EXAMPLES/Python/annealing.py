#!/usr/bin/env python
# annealing.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
"""Minimize the energy of met-enkaphalin starting from a given structure using
simulated annealing. In simulated annealing the temperature of a Monte Carlo
simulation is slowly reduced from Tmax to Tmin.
"""
# Adds the source directory to Python's search path.
import sys
sys.path.append('../..')
import smmp, universe, protein

# Initialize the Universe to T=300K with the ECEPP/3 force field, no solvent 
# term (st = 0) and the sub directory SMMP/ as library path. Except for the
# solvent term, these are the default values. Alternatively, we could have 
# written
# myUniverse = universe.Universe(st=0)
# to get the same result.
myUniverse = universe.Universe(T=300, ff = 'ecepp2', st = 0, libdir ='SMMP/')
# Create a new protein object from the sequence file ../enkefa.seq and
# set the dihedral angles according to the values given in ../enkefa.var.
p = protein.Protein('../enkefa.seq', '../enkefa.ann')
# Make myUniverse aware of p.
myUniverse.add(p)
seed = 81236
smmp.sgrnd(seed)

Tmin = 100
Tmax = 1000
equilibrationSweeps = 100
sweeps = 100000
measurementInterval = 1000
randomStart = 1
smmp.anneal(equilibrationSweeps, sweeps, measurementInterval, Tmax, Tmin, randomStart)
