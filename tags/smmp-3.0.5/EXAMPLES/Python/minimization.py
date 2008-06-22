#!/usr/bin/env python
# minimization.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
"""Minimize the energy of met-enkaphalin starting from a given structure using
conjugate gradient. 
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
myUniverse = universe.Universe(T=300, ff = 'ecepp3', st = 0, libdir ='SMMP/')
# Create a new protein object from the sequence file ../enkefa.seq and
# set the dihedral angles according to the values given in ../enkefa.var.
p = protein.Protein('../enkefa.seq', '../enkefa.var')
# Make myUniverse aware of p.
myUniverse.add(p)
typeOfMinimization = 1
maxNumberOfIterations = 10000
desiredPrecision = 1.0e-9

smmp.minim(typeOfMinimization, maxNumberOfIterations, desiredPrecision)
