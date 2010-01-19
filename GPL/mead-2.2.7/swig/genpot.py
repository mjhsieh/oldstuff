from PyMead import MEAD
import Numeric

def calcpot(alist):
    """Calculate potential grid of molecule and return as Numeric array.
    Argument is a list of Atoms (no duplicates)."""

    epsin = 4.0
    istrength = 0.4
    eext = 80.0
    rad = 1.4
    exclus_radius = 2.0

    a = MEAD.AtomSet("Calcpot Atoms")
    for at in alist:
      a.insert(at)

    at = alist[0]
    print 'The atom used is ', at.atname

    rho = MEAD.AtomChargeSet (a)
    eps = MEAD.TwoValueDielectricByAtoms (a, epsin, eext, rad)
    ely = MEAD.ElectrolyteByAtoms (a, istrength, exclus_radius)

    fdm = MEAD.FinDiffMethod ()
    fdm.add_level(41, 1.0, a.geom_cent())
    phi = MEAD.FinDiffElstatPot (fdm, eps, rho, ely)
    phi.solve()

    print "Solving for phi done."
    cls = phi.coarse_lattice_spec()
    print "grid dimension is", cls.get_grid_dim(), "cubed"
    print "grid spacing is", cls.get_spacing(), "Angstroms"
    c = cls.get_center()
    print "grid center is point", c.x, c.y, c.z


    cf = phi.get_cuberep (cls) # returns Numeric array
    print "Shape of array is", cf.shape
    print "maximum value is", Numeric.maximum.reduce(cf.flat)
    print "minimum value is", Numeric.minimum.reduce(cf.flat)

    return cf

MEAD.Blab.level = 2

crd = MEAD.Coord(0.0, 0.0, 0.0)
at = MEAD.Atom()
at.atname = "NA"
at.resname = "ION"
at.resnum = 1
at.coord = crd
at.rad = 2.0
at.charge = 1.0

for i in range(1,2):
   print 'Calling calcpot at iteration ', i
   cf = calcpot( [at] )
   cf.shape
