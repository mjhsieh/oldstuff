"""Incomplete example.  For illustrative purposes only."""

import MEAD

def pair_analysis(phi1_sol, phi1_vac, phi1_alone_sol, phi1_alone_vac,
                  phi2_sol, phi2_vac, phi2_alone_vac, phi2_alone_sol,
                  chrg1, chrg2):
    """The phi should be ElstatPot, and the chrg ChargeDist.
    Their underlying specific type (e.g. AtomChargeDist) doesn't matter"""

    DGsol1 = (MEAD.energy_of(chrg1, phi1_alone_sol)
              - MEAD.energy_of(chrg1, phi1_alone_vac)) / 2.0
    DGsol2 = (MEAD.energy_of(chrg2, phi2_alone_sol)
              - MEAD.energy_of(chrg2, phi2_alone_vac)) / 2.0

    DGsol1_diel2 = (MEAD.energy_of(chrg1, phi1_sol)
                    - MEAD.energy_of(chrg1, phi1_vac)) / 2.0
    DGsol2_diel1 = (MEAD.energy_of(chrg2, phi2_sol)
                    - MEAD.energy_of(chrg2, phi2_vac)) / 2.0

    DGdesol1 = DGsol1 - DGsol1_diel2
    DGdesol2 = DGsol2 - DGsol2_diel1
    print "desolvation of molecules 1 and 2 are:", DGdesol2, DGdesol2


def use_pair_analysis(atset1, atset2, both_atsets):
    # create phi_sol, etc. as as in solven

    # but for reference "vacuum" state
    epsin = 4.0
    ely_ref = MEAD.ElectrolyteEnvironment \
                 (MEAD.UniformElectrolyte (0.0))
    eps_ref = MEAD.DielectricEnvironment \
                 (MEAD.UniformDielectric (epsin))
    phi1_ref = MEAD.ElstatPot (eps_ref, chrg1, ely_ref)
    # ph1_ref is now a simple Coulombic potential
    phi1_alone_ref = phi_1
    
    # .... call to pair_analysis using phi1_ref, etc.
