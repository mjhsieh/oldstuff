from MEAD import *

a = AtomSet()
crd = Coord(0.0, 0.0, 0.0)
a.build_from_vectors(["ION"], ["NA"], [1], [crd], [2.0], [1.0])

epsin = 4.0
istrength = 0.4
eext = 80.0
rad = 1.4
exclus_radius = 2.0
interesting = Coord (0, 0, 0)

rho = AtomChargeSet (a)
eps = TwoValueDielectricByAtoms (a, epsin, eext, rad)
ely = ElectrolyteByAtoms (a, istrength, exclus_radius)

fdm = FinDiffMethod ()
fdm.add_level(41, 1.0, a.geom_cent())
fdm.add_level(41, 0.25, a.geom_cent())
#fdm.read ("born.ogm")
#fdm.resolve (a.geom_cent(), interesting);

phi = FinDiffElstatPot (fdm, eps, rho, ely)
phi.solve()

print "The potential in solvent is now calculated"
print "getting info about course lattice"
cls = phi.coarse_lattice_spec()
print "grid dimension is", cls.get_grid_dim(), "cubed"
print "grid spacing is", cls.get_spacing(), "Angstroms"
c = cls.get_center()
print "grid center is point", c.x, c.y, c.z

print "now pull out the course field"
cf = phi.get_cuberep (cls)
print "Got it.  It's shape is", cf.shape
import Numeric
_max = Numeric.maximum.reduce
print "maximum value is", _max(_max(_max(cf)))
print "maximum value is", Numeric.maximum.reduce(cf.flat)
print "minimum value is", Numeric.minimum.reduce(cf.flat)

prod_sol = rho * phi 
print "prod_sol = ", prod_sol


vac_istrength = 0.0
vac_eext = 1.0
econv = 331.842

elyvac = ElectrolyteByAtoms (a, vac_istrength, exclus_radius)
vac_eps = TwoValueDielectricByAtoms (a, epsin, vac_eext, rad)

vac_phi = FinDiffElstatPot (fdm, vac_eps, rho, elyvac)
vac_phi.solve()

prod_vac = rho * vac_phi
print "prod_vac = ", prod_vac

print "solvation_energy = ", (prod_sol - prod_vac) / 2 * econv
