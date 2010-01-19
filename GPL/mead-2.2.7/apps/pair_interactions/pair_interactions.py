import sys
#
# FOR TESTING
#
#sys.path.insert(0,'../../../../swig')

import getopt
import types
import copy
import MEAD

# Set all charges in AtomSet ats to zero

def zero_charges(ats):
  if not isinstance(ats, MEAD.AtomSet):
    raise TypeError, 'zero_charges: argument must be an AtomSet'
  for atomID in ats.keys():
     ats[atomID].charge = 0.0


# Merge two AtomSets together

def merge_atoms(a1, a2):
  if not isinstance(a1, MEAD.AtomSet):
    raise TypeError, 'merge_atoms: first argument must be an AtomSet'
  if not isinstance(a2, MEAD.AtomSet):
    raise TypeError, 'merge_atoms: second argument must be an AtomSet'
  ats = copy.deepcopy(a1)
  ats.update(a2)
  return ats


# Print out the cubic lattic levels in the finite difference method

def print_fdm_levels(fdm):
  if not isinstance(fdm, MEAD.FinDiffMethod):
    raise TypeError, 'print_fdm_levels: first argument must be a FinDiffMethod'
  print 'Using finite difference method with lattice levels:'
  print ' Level             Center               Grid Dimension    Spacing'
  l=1
  for level in fdm.levels():
    cntr = level.get_center()
    print '%4d %s %0.3f %s %0.3f %s %0.3f %s %4d %s %0.3f' % (l, '     (', cntr.x, ',', cntr.y, ',', cntr.z, ')       ', level.get_grid_dim(), '          ', level.get_spacing())
    l = l + 1

# Compute interactions between molecules using two calculations of
# electrostatic potential for each molecules' charges zeroed out

def pair_interactions(a1, a2, fdm, epsin, ionicstr):
  if not isinstance(a1, MEAD.AtomSet):
    raise TypeError, 'pair_interactions: first argument must be an AtomSet'
  if not isinstance(a2, MEAD.AtomSet):
    raise TypeError, 'pair_interactions: second argument must be an AtomSet'
  if not isinstance(fdm, MEAD.FinDiffMethod):
    raise TypeError, 'pair_interactions: third argument must be a FinDiffMethod'

# Define a new atom set with zero charges
  a1zero = MEAD.AtomSet(a1)
  zero_charges(a1zero)

# Define a new atom set with zero charges
  a2zero = MEAD.AtomSet(a2)
  zero_charges(a2zero)

# Resolve the lattice levels to the geometric center of the two molecules
  if (not fdm.is_resolved()):
    geom_cent1 = a1.geom_cent()
    geom_cent2 = a2.geom_cent()
    geom_cent = geom_cent1 + (geom_cent2 - geom_cent1) / 2.0
    center_of_interest = geom_cent
    fdm.resolve(geom_cent, center_of_interest)

# print_fdm_levels(fdm)
  fdm.write()

# Electrolyte set to whatever the user passes in, but probably should be zero
  if ionicstr != 0:
    print '%s %0.3f %s' % ('pair_interactions: WARNING ionic strength of', ionicstr, 'in Uniform Electrolyte')

  ely = MEAD.UniformElectrolyte(ionicstr)

# Solve the electrostatic potential with the first molecule charges zeroed out
  ats = merge_atoms(a1zero, a2)
  eps = MEAD.TwoValueDielectricByAtoms(ats, epsin)
  acs = MEAD.AtomChargeSet(ats)
  phi1 = MEAD.FinDiffElstatPot(fdm, eps, acs, ely)
  phi1.solve()
  print 'Calculation with molecule 1 charges zeroed is done'

# Solve the electrostatic potential with the second molecule charges zeroed out
  ats.clear()
  ats = merge_atoms(a1, a2zero)
  eps = MEAD.TwoValueDielectricByAtoms(ats, epsin)
  acs = MEAD.AtomChargeSet(ats)
  phi2 = MEAD.FinDiffElstatPot(fdm, eps, acs, ely)
  phi2.solve()
  print 'Calculation with molecule 2 charges zeroed is done'

  return phi1, phi2

if __name__ == '__main__':
  #
  # MAIN

  # Process the command line arguments
  #
  opts, pargs = getopt.getopt(sys.argv[1:], '', ['epsin=', 'epsext=', 'solrad=', 'sterln=', 'ionicstr=', 'blab1', 'blab2', 'blab3'])

  for opt in opts:
    if opt[0] == '--epsin':
      try:
        epsin = float(opt[1])
      except:
        raise ValueError, 'epsin must be a number'
    elif opt[0] == '--epsext':
      try:
        MEAD.PhysCond.epsext = float(opt[1])
      except:
        raise ValueError, 'epsext must be a number'
    elif opt[0] == '--solrad':
      try:
        MEAD.PhysCond.solrad = float(opt[1])
      except:
        raise ValueError, 'solrad must be a number'
    elif opt[0] == '--sterln':
      try:
        MEAD.PhysCond.sterln = float(opt[1])
      except:
        raise ValueError, 'sterln must be a number'
    elif opt[0] == '--ionicstr':
      try:
        MEAD.PhysCond.ionicstr = float(opt[1])
      except:
        raise ValueError, 'ionicstr must be a number'
    elif opt[0] == '--blab1':
      MEAD.Blab.level = 1
    elif opt[0] == '--blab2':
      MEAD.Blab.level = 2
    elif opt[0] == '--blab3':
      MEAD.Blab.level = 3
    else:
      print 'Unknown option:', opt[0], 'is ignored'

  # You must supply the epsin value on the command line
  try:
    epsinid = id(epsin)
  except:
    raise ValueError, 'epsin must be defined using --epsin <value>'

  # epsin must be a number
  assert (isinstance(epsin, types.IntType) or isinstance(epsin, types.FloatType)), 'epsin must be a Number type'

  # Better have two molecule names on the command line
  if len(pargs) != 2:
    raise RuntimeError, 'must supply two molecule names as input'

  molname1 = pargs[0]
  molname2 = pargs[1]

  # Ionic strength will be zero by default unless changed by the user
  ionicstr = MEAD.PhysCond.ionicstr

  econv = MEAD.PhysCond.econv

  # Prepare output files with names derived from the two molecule names
  m1m2 = molname1 + '_' + molname2
  interaction_1_2_filename = m1m2 + '_interaction_1_2.dat'
  interaction_2_1_filename = m1m2 + '_interaction_2_1.dat'

  interaction_1_2_file = open(interaction_1_2_filename, 'w')
  interaction_2_1_file = open(interaction_2_1_filename, 'w')

  # OK, let's begin
  print 'Starting', sys.argv[0], 'for pair of molecules named', molname1, 'and', molname2
  print 'with interior dielectric constant =', epsin
  print 'using the following physical conditions:'
  MEAD.PhysCond.write()
  if MEAD.Blab.level > 0:
    print 'Blab level set to', MEAD.Blab.level
  else:
    print 'No blab level set (so no blabbing)'

  # Read in the first molecule
  a1 = MEAD.AtomSet(molname1)
  a1.read()

  # Read in the second molecule
  a2 = MEAD.AtomSet(molname2)
  a2.read()

  # Read in the finite difference method lattice levels
  fdm = MEAD.FinDiffMethod()
  ogm_filename = molname1 + '.ogm'
  fdm.read(ogm_filename)

  # Compute the interactions between the two molecules
  phi1, phi2 = pair_interactions(a1, a2, fdm, epsin, ionicstr)

  # Write the interaction energies for each atom in molecule 1 with molecule 2
  for atomID in a1.keys():
    a = a1[atomID]
    line = a.atname + ' ' + a.resname + ' ' + str(a.resnum) + ' ' + str(a.charge * phi1.value(a.coord) * econv) + '\n'
    interaction_1_2_file.write(line)

  # Write the interaction energies for each atom in molecule 2 with molecule 1
  for atomID in a2.keys():
    a = a2[atomID]
    line = a.atname + ' ' + a.resname + ' ' + str(a.resnum) + ' ' + str(a.charge * phi2.value(a.coord) * econv) + '\n'
    interaction_2_1_file.write(line)

  # All done!
  interaction_1_2_file.close()
  interaction_2_1_file.close()
