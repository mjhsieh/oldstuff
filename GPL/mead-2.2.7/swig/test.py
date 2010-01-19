r'''>>> from PyMead.MEAD import *
    >>> from copy import *
    >>> a = Atom()
    >>> c = Coord(1,2,3)
    >>> a.coord = c
    >>> a.atname = "AT"
    >>> print a.coord.x, a.coord.y, a.coord.z
    1.0 2.0 3.0
    >>> print a.atname
    AT
    >>> b = deepcopy(a)
    >>> d = Coord(2,3,4)
    >>> at3 = AtomID(1, "AT", "CID")
    >>> at2 = AtomID(1, "AT")
    >>> at2.resnum
    1
    >>> a.atname = "A"
    >>> b = Atom()
    >>> b.atname = "B"
    >>> al = [a,b]
    >>> ats = AtomSet(al)
    >>> id = AtomID(a)
    >>> at5 = ats[id]
    >>> a2 = Atom()
    >>> a2.atname = "C"
    >>> id = AtomID(a2)
    >>> acs = AtomChargeSet(ats)
    >>> lk = acs.keys()
    >>> lv = acs.values()
    >>> lp = acs.items()
    >>> lpt = acs.pointcharges()
    >>> acs = AtomChargeSet(ats)
    >>> acs2 = AtomChargeSet(acs)
    >>> pc1 = PointCharge()
    >>> pc1.charge=1
    >>> pc2 = PointCharge()
    >>> pc2.charge=2
    >>> mpc = ManyPointCharge([pc1,pc2])
    >>> mpc2 = acs + mpc
    >>> opc = OnePointCharge(1.,d)
    >>> opc2 = OnePointCharge(2.,d)
    >>> mpc3 = opc + opc2
    >>> pc3 = mpc3.pointcharges()
    >>> pc3[0].charge
    3.0
    >>> opc3 = OnePointCharge(3,c)
    >>> mpc4 = opc + opc3
    >>> pc4 = mpc4.pointcharges()
    >>> pc4[0].charge
    1.0
    >>> pc4[1].charge
    3.0
    >>> pc3[0].coord.x
    2.0
    >>> ud = UniformDielectric(80.)
    >>> dsp = DielectricSphere(4.0, 80.0, 10.0, c)
    >>> ue = UniformElectrolyte(1.0)
    >>> esp = ElySphere(1.0, c, 10.0)
    >>> cls = CubeLatSpec(3, 5.0, c)
    >>> ecr = esp.get_cuberep(cls)
    >>> fdm = FinDiffMethod()
    >>> fdm.add_level(3, 5.0, c)
    >>> asp = AnalySphere(dsp, mpc, esp)
    >>> fde = FinDiffElstatPot(fdm, ud, acs, ue)
    >>> asp.solve()
    >>> fde.solve()
    >>> fde.value(c)
    0.0
    >>> asp.value(c)
    0.126329064369
    >>> asp.value(d)
    0.0695440098643
    >>> fde.value(d)
    0.0
    >>> epc = asp + fde
    >>> epc.value(d)
    0.0695440098643
    >>> fepc = epc.field(c)
    >>> fasp = asp.field(c)
    >>> epc4 = asp * 4.
    >>> epc4 = 4. * asp
    >>> epc4.value(c)
    0.505316257477
    >>> epc4.value(c) == asp.value(c) * 4.
    1
'''

def run(args = None):
    if args is not None:
        import sys
        sys.argv = args
    import doctest, test
    print doctest.testmod(test, verbose=1)

if __name__ == '__main__':
    run()
