from PyMead import MEAD
econv = 331.842 # conversion factor for elem_crg^2/Angstroms to kcal/mole

def solven(atset, fdm, epsin, eps_sol):
    
    istrength = 0.0
    eext = eps_sol
    rad = 1.4
    exclus_radius = 2.0

    rho = MEAD.AtomChargeSet (atset)
    print "Constructing new TwoValueDielectricByAtoms for solvent"
    eps = MEAD.TwoValueDielectricByAtoms (atset, epsin, eext, rad)
    print "Constructing new ElectrolyteByAtoms for solvent"
    ely = MEAD.ElectrolyteByAtoms (atset, istrength, exclus_radius)

    print "Constructing new FinDiffElstatPot for solvent"
    phi_sol = MEAD.FinDiffElstatPot (fdm, eps, rho, ely)
    phi_sol.solve()

    prod_sol = rho * phi_sol

    vac_istrength = 0.0
    vac_eext = 1.0

    print "Constructing new ElectrolyteByAtoms for vacuum"
    elyvac = MEAD.ElectrolyteByAtoms (atset, vac_istrength, exclus_radius)
    print "Constructing new TwoValueDielectricByAtoms for vacuum"
    vac_eps = MEAD.TwoValueDielectricByAtoms (atset, epsin, vac_eext, rad)

    print "Constructing new FinDiffElstatPot for vacuum"
    vac_phi = MEAD.FinDiffElstatPot (fdm, vac_eps, rho, elyvac)
    vac_phi.solve()

    prod_vac = rho * vac_phi
    print 'prod_vac = ', prod_vac

    return (prod_sol - prod_vac) / 2 * econv

# Now a very simple example of using to compute solvation of an ion

# Specify things about the physical system
crd = MEAD.Coord(0.0,0.0,0.0)
chrg = 1.0
rad = 2.0
at = MEAD.Atom()
at.atname = "NA"
at.resname = "ION"
at.resnum = 1
at.coord = crd
at.rad = rad
at.charge = chrg
ats = MEAD.AtomSet("Sodium Ion")
ats.insert(at)
epssol = 80.0
eps_in = 2.0

# Describe the set of lattices to be used in the finite difference method
fdm = MEAD.FinDiffMethod ()
fdm.add_level(41, 1.0, MEAD.CenteringStyle.ON_ORIGIN)  # the coarsest level
fdm.add_level(41, 0.25, MEAD.CenteringStyle.ON_ORIGIN) # a finer one in a smaller region.
fdm.resolve(crd,crd)

for i in range(1, 2):
   print "calling solven at iteration ", i
   result = solven(ats, fdm, eps_in, epssol)
   print "the answer computed by solven is:", result

   print "For comparison, the analytical result computed using"
   print "the Born formula is:",
   print - econv / (2*2.0) * (1.0 - 1.0/epssol)
