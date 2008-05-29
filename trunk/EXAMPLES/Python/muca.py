#!/usr/bin/env python
# muca.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
import sys

sys.path.append('../..')

import smmp, universe, protein

smmp.epar_l.flex = 0
smmp.epar_l.sh2 = 0
smmp.epar_l.epsd = 0
smmp.epar_l.ientyp = 0
smmp.isolty.itysol = 0
smmp.init_energy('./SMMP/')
smmp.mol_i.ntlml=0

smmp.sgrnd(31433)
smmp.updchois.upchswitch = 0
smmp.updchois.rndord = 0
smmp.updchois.bgsprob = 0.3
smmp.init_lund()

p = protein.Protein('../enkefa.seq', '../enkefa.var')

smmp.multicanonical.mulcan_par(100000, 500, 1000, -12, 20, 1.0, 0)
#smmp.multicanonical.mulcan_sim(100, 100000, 10, 1000, -12, 20, 1.0, 0)

