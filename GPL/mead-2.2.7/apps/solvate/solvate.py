#! /usr/bin/env python

# solvate.py
# ID: $Id: solvate.py,v 1.5 2000/08/02 16:32:54 ttnttn Exp $
# SOURCE: $Source: /cvs-repository/bashford/cvsroot/mead/apps/solvate/solvate.py,v $

import sys
from MEAD import *

# these should go into option processing
molname = sys.argv[1]
print 'molname = ', molname

epsin = 4.0
istrength = 0.4
eext = 80.0
rad = 1.4
exclus_radius = 2.0
interesting = Coord (0, 0, 0)

# this is what we need to interface to from MMTK/PMV
# one approach is to extend AtomSet
a = AtomSet ()
a.read (molname + ".pqr")

rho = ChargeDist (AtomChargeSet (a))
eps = DielectricEnvironment (TwoValueDielectricByAtoms (a, epsin, eext, rad))
ely = ElectrolyteEnvironment (ElectrolyteByAtoms (a, istrength, exclus_radius))

fdm = FinDiffMethod ()
fdm.read (molname + ".ogm")
fdm.resolve (a.geom_cent(), interesting);

phi = ElstatPot (fdm, eps, rho, ely)
phi.solve()

prod_sol = energy_of (rho, phi)         # do NOT reverse args!
print "prod_sol = ", prod_sol


vac_istrength = 0.0
vac_eext = 1.0
econv = 331.842

elyvac = ElectrolyteEnvironment (ElectrolyteByAtoms (a, vac_istrength,
                                                     exclus_radius))
vac_eps = DielectricEnvironment (TwoValueDielectricByAtoms (a, epsin,
                                                            vac_eext, rad))

vac_phi = ElstatPot (fdm, vac_eps, rho, elyvac)
vac_phi.solve()

prod_vac = energy_of (rho, vac_phi)
print "prod_vac = ", prod_vac

print "solvation_energy = ", (prod_sol - prod_vac) / 2 * econv

# solvate.py ends here
