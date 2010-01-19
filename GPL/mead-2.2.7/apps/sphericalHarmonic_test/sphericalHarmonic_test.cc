#include "MEAD/SphericalHarmonic.h"
#include "MEAD/globals.h"
#include "MEAD/MEADexcept.h"
#include "MEAD/AtomChargeSet.h"
#include <iostream>

#include <fstream>
#include <sstream>
#include <math.h>

#include <typeinfo>
using std::type_info;

using std::vector;

template<class FunType> void plotfun (const FunType& funob, char * datname)
{
  std::ofstream ofs(datname);
  for (double x = -1.0; x < 1.0; x += 0.01)
    ofs << x << "   " << funob(x) << endl;
  ofs << 1.0 << "   " << funob(1.0) << endl;
}  

inline double norm(double x) {return x*x;}

// Trapazoid-rule integration of funob**2 from -1.0 to 1.0 using ninter points
template<class FunType> double normal (const FunType& funob,
				       unsigned ninter=2000)
{
  double dx = 2.0/ninter;
  double prenorm = dx*norm(funob(-1.0))/2.0;
  //  for (double x = -1.0; x < 1.0; x += dx) {
  for (int i=1; i<ninter; ++i) {
    double x = -1.0 + i*dx;
    prenorm += norm(funob(x))*dx;
  }
  prenorm += dx*norm(funob(1.0))/2.0;
  prenorm += dx*funob(1.0)*funob(1.0)/2.0;
  return prenorm;
}  

// like normal, above, but for a function defined on unit sphere.
template<class FunType> double sphenorm (const FunType& funob,
					 unsigned ninter=100)
{
  int theta_inters = ninter;
  int phi_inters = ninter;
  double dtheta = M_PI/theta_inters;
  double dphi = 2.0*M_PI/phi_inters;
  // Trapezoid formula over theta, but using fact that integrand
  // fun*sin(theta) is formally zero at poles.
  double prenorm = 0.0;
  for (int itheta=1; itheta < theta_inters; ++itheta) {
    double theta = itheta*dtheta;
    // for the circular phi integration the usual factors of 1/2
    // at the ends in the trapezoid rule don't apply.
    double phisum = 0.0;
    for (int iphi=0; iphi < phi_inters; ++iphi) {
      double phi = iphi*dphi;
      phisum += norm(funob(theta, phi)) * dphi;
    }
    prenorm += phisum*sin(theta)*dtheta;
  }
  return prenorm;
}

inline double sphunitfun(double theta, double phi) {return 1.0;}

#include <typeinfo>

template<class FunType>
void do_normtest(int max_ell,
		 double expected_norm = 1.0,
		 double (*normfunctional)(const FunType&, unsigned) = normal<FunType>)
{
  const type_info* ti;
  ti = &typeid(FunType);
  const char* tn = ti->name();
  cout << "entering do_normtest for FunType, " << tn << endl;
  cout << "max_ell = " << max_ell << endl;
  for (int ell=0; ell <= max_ell; ++ell) {
    cout << "ell =" << ell << endl;
    for (int m=0; m <= ell; ++m) {
      FunType fun(ell,m);
      int numnodes = ell - m; // I think that's how many nodes func has
      //      // Use 2000 integration points, or 200 per node, whichever is more.
      //      int ninter = numnodes > 10 ? 200*numnodes : 2000;
      // Use 100 integration points, or 10 per node, whichever is more.
      int ninter = numnodes > 10 ? 10*numnodes : 100;
      double norm = normfunctional(fun, ninter);
      double dev = expected_norm - norm;
      if (fabs(dev) > 0.01) {
	cerr << "deviation of " << dev << " for ell = " << ell
	     << " m = " << m << endl;
	return;
      }
    }
  }
}


void PreSpheHarm_normtest(int max_ell)
{
  cout << "entering PreSpheHarm_normtest with max_ell = " << max_ell << endl;
  cout << "start building the series" << endl;
  vector< vector < PreSpheHarm > > vvs;
  make_PreSpheHarm_series(max_ell, vvs);
  cout << "finished building the series" << endl;
  for (int ell=0; ell <= max_ell; ++ell) {
    cout << "ell =" << ell << endl;
    for (int m=0; m <= ell; ++m) {
      int numnodes = ell - m; // I think that's how many nodes func has
      // Use 2000 integration points, or 200 per node, whichever is more.
      int ninter = numnodes > 10 ? 200*numnodes : 2000;
      double norm = normal(vvs[ell][m], ninter);
      double dev = 1.0 - norm*2*M_PI;
      if (fabs(dev) > 0.01) {
	cerr << "deviation of " << dev << " for ell = " << ell
	     << " m = " << m << endl;
      }
    }
  }
}

ostream& angfracprin(int num, int den, ostream& ost)
{
  if (num < 0 && den > 0) {ost << "-"; num = -num;}
  if (num > 0 && den < 0) {ost << "-"; den = -den;}
  if (num < 0 && den < 0) {den = -den; num = -num;}
  if (num == den)
    ost << "pi";
  else {
    ost << num;
    if (num != 0) ost << "pi/" << den;
  }
  return ost;
}

#include <MEAD/MomentAnalysis.h>
#include <MEAD/MomentAnalysis_tmplts.h>

struct NormY {
  NormY(const SphericalHarmonic &y)
    : _y(y) {}
  complex<double> operator()(double theta, double phi) const
  {
    complex<double> yval = _y(theta, phi);
    return conj(yval) * yval;
  }
  SphericalHarmonic _y;
};



#include "MEAD/ManyPointCharge.h"
#include "MEAD/UniformDielectric.h"
#include "MEAD/DielectricSphere.h"
#include "MEAD/UniformElectrolyte.h"
#include "MEAD/ElySphere.h"
#include "MEAD/ElstatPot.h"
#include "MEAD/Debye.h"
#include "MEAD/AnalySphere.h"


void testAnal(ManyPointCharge* rhop)
{
  const int maxell=3;
  ChargeDist rho(rhop);
  cout << "The charge distribution has spherical moments" << endl;
  Moments rhomoms = momentsOfChargeDist(rho, maxell);
  print_moments(rhomoms);

  cout << "The potential of the charge distribution in vaccum ..." << endl;
  ElstatPot vacphi(DielectricEnvironment(new UniformDielectric(1.0)),
		   ChargeDist(rhop),
		   ElectrolyteEnvironment(new UniformElectrolyte(0.0)));
  vacphi.solve();
  cout << "has effective spherical moments" << endl;
  Moments vacmoms = momentsOfElstatPot(vacphi, maxell, 5.0);
  print_moments(vacmoms);

  float dielrad = 4.0;
  float epsin = 8.0;
  float epsext = 1.0;
  cout << "its potential inside sphere of radius " << dielrad
       << ", inner dielectric " << epsin << " (outside vacuum)" << endl;
  
  ElstatPot
    sphphi(new AnalySphere(new DielectricSphere(epsin, epsext, dielrad,
						Coord(0,0,0)),
			   rhop,
			   new ElySphere(0.0, Coord(0,0,0), dielrad),
			   5));
  sphphi.solve();
  cout << "has effective spherical moments" << endl;
  Moments sphmoms = momentsOfElstatPot(sphphi, maxell, 5.0);
  print_moments(sphmoms);

  compare_moments(vacmoms, rhomoms, "vacuum_pot", "charge_dist");
  compare_moments(sphmoms, rhomoms, "sphere_pot", "charge_dist");
  cout << "for a sphere with epsin = " << epsin << ", epsext = " << epsext
       << ", expected ratios are:" << endl;
  for (int ell=0; ell <= maxell; ++ell) {
    cout << ell << ":  " <<
      (1.0 + (epsin - epsext)*(ell+1)/(ell*epsin + epsext*(ell+1)))/epsin
	 << endl;
  }
}


main()
{
  
  try {
    ManyPointCharge* combo = new ManyPointCharge();

    cout << "SINGLE UNIT CHARGE AT ORIGIN" << endl;
    ManyPointCharge* rhop = new ManyPointCharge();
    rhop->push_back(PointCharge(Coord(0,0,0), 1.0));
    combo->insert(combo->end(), rhop->begin(), rhop->end());
    testAnal(rhop);

    cout << "DIPOLE OF 2 ALONG Z AXIS" << endl;
    rhop = new ManyPointCharge();
    rhop->push_back(PointCharge(Coord(0,0,1.0), 1.0));
    rhop->push_back(PointCharge(Coord(0,0,-1.0), -1.0));
    combo->insert(combo->end(), rhop->begin(), rhop->end());
    testAnal(rhop);

    cout << "DIPOLE OF 2 ALONG X AXIS" << endl;
    rhop = new ManyPointCharge();
    rhop->push_back(PointCharge(Coord( 1.0,  0.0,  0),  1.0));
    rhop->push_back(PointCharge(Coord(-1.0,  0.0,  0), -1.0));
    combo->insert(combo->end(), rhop->begin(), rhop->end());
    testAnal(rhop);

    cout << "DIPOLE OF 2 ALONG Y AXIS" << endl;
    rhop = new ManyPointCharge();
    rhop->push_back(PointCharge(Coord( 0.0,  1.0,  0),  1.0));
    rhop->push_back(PointCharge(Coord( 0.0, -1.0,  0), -1.0));
    combo->insert(combo->end(), rhop->begin(), rhop->end());
    testAnal(rhop);

    cout << "THE COMBINATION OF ALL THE ABOVE" << endl;
    testAnal(combo);

    // What's below might be irrelevant.

    cout << "the area of the unit sphere is " << sphenorm(sphunitfun) << endl;
    SphericalHarmonic Y21(2,1);
    cout << "Y21: " << Y21 << endl;
    cout << "norm of Y21: " << sphenorm(Y21) << endl;
    Moments moms = sphericalComponents(Y21, 2);
    for (int ell=0; ell <= 2; ++ell) {
      for (int m=-ell; m <= ell; ++m) {
	Moments::momtype mom = moms(ell, m);
	cout << ell << "   " << m << "   ";
	if (norm(mom) < 1.0e-15)
	  cout << "zero";
	else
	  cout << mom;
	cout << endl;
      }
    }




    cout << "\nSYSTEMATIC TEST OF NORMS OF SphericalHarmonics" << endl;
    cout << "Sorry the test of norms is broken" << endl;
    do_normtest(50, 1.0, sphenorm<SphericalHarmonic>);


    cout << "test of spherint with sphunitfun" << endl;
    cout << spherint<double>(sphunitfun) << endl;

    cout << "test of spherint with NormY(Y21): " << endl;
    NormY n21(Y21);
    cout << spherint< complex<double> > (n21) << endl;

    PreSpheHarm spy30(3,0);
    cout << "spy30: " << spy30  << endl;
    cout << "      normal*2*pi = " << normal(spy30)*2*M_PI << endl;
    plotfun(spy30, "spy30.dat");
    cout << endl;

    PreSpheHarm spy31(3,1);
    cout << "spy31: " << spy31  << endl;
    cout << "      normal*2*pi = " << normal(spy31)*2*M_PI << endl;
    plotfun(spy31, "spy31.dat");
    cout << endl;

    PreSpheHarm spy32(3,2);
    cout << "spy32: " << spy32  << endl;
    cout << "      normal*2*pi = " << normal(spy32)*2*M_PI << endl;
    plotfun(spy32, "spy32.dat");
    cout << endl;

    PreSpheHarm spy33(3,3);
    cout << "spy33: " << spy33  << endl;
    cout << "      normal*2*pi = " << normal(spy33)*2*M_PI << endl;
    plotfun(spy33, "spy33.dat");
    cout << endl;

    Legendre P_5(5);
    cout << "P_5:           " << P_5 << endl;

    PreSpheHarm spy50(5,0);
    cout << "spy50: " << spy50  << endl;
    cout << "      normal*2*pi = " << normal(spy50)*2*M_PI << endl;
    plotfun(spy50, "spy50.dat");
    cout << endl;

    PreSpheHarm spy53(5,3);
    cout << "spy53: " << spy53  << endl;
    cout << "      normal*2*pi = " << normal(spy53)*2*M_PI << endl;
    plotfun(spy53, "spy53.dat");
    cout << endl;

    PreSpheHarm spy54(5,4);
    cout << "spy54: " << spy54  << endl;
    cout << "      normal*2*pi = " << normal(spy54)*2*M_PI << endl;
    plotfun(spy54, "spy54.dat");
    cout << endl;

    PreSpheHarm spy25_21(25,21);
    cout << "spy25_21: " << spy25_21  << endl;
    cout << "      normal*2*pi = " << normal(spy25_21)*2*M_PI << endl;
    plotfun(spy25_21, "spy25_21.dat");

    PreSpheHarm spy39_0(39,0);
    cout << "spy39_0: " << spy39_0 << endl;
    plotfun(spy39_0, "spy39_0");

    PreSpheHarm spy47_0(47,0);
    cout << "spy47_0: " << spy47_0 << endl;
    plotfun(spy47_0, "spy47_0");


    cout << "BEGINNING QUICK TESTING OF PRESPHEHARM" << endl;
    PreSpheHarm_normtest(50);
    cout << "SYSTEMATIC NORMALIZATION TESTS OF PreSpheHarm" << endl;
    do_normtest<PreSpheHarm>(50);
    cout << "done with systematic tests" << endl;

  }
  catch(MEADexcept e) {
    cerr << e.get_error1() << e.get_error2() << e.get_error3() << endl;
    exit(1);
  }

  cout << "try the old stuff" << endl;
  AssocLegendre P00;
  AssocLegendre P20(2,0);
  AssocLegendre P21(2,1);
  AssocLegendre P22(2,2);
  cout << P00 << endl;
  cout << P20 << endl;
  cout << P21 << endl;
  cout << P22 << endl;
  plotfun(P00, "P00.dat");
  plotfun(P20, "P20.dat");
  plotfun(P21, "P21.dat");
  plotfun(P22, "P22.dat");
  AssocLegendre P25_21(25,20);
  cout << "P25_21 : " << P25_21 << endl;
  plotfun(P25_21, "P25_21.dat");


  Poly p(1.0, -1.0);
  cout << p << endl;
  cout << "size " << p.size() << endl;
  cout << "p*p    :       " << p*p << endl;
  cout << "p + p  :       " << p + p << endl;
  Legendre P_0;
  cout << "P_0:           " << P_0 << endl;
  Legendre P_0_also(0);
  cout << "P_0_also:      " << P_0_also << endl;
  Legendre P_1 = order_one_Legendre();
  cout << "P_1:           " << P_1 << endl;
  Legendre P_1_also(1);
  cout << "P_1_also:      " << P_1_also << endl;

  cout << "P_1*P_1:       " << P_1*P_1 << endl;
  cout << "P_0 + P_1*P_1: " << P_0 + P_1*P_1 << endl;

  Legendre P_2(2);
  cout << "P_2:           " << P_2 << endl;
  Legendre P_3(3);
  cout << "P_3:           " << P_3 << endl;
  Legendre P_4(4);
  cout << "P_4:           " << P_4 << endl;

  int ellmax = 5;
  cout << "MAKE THE SERIES UP TO " << ellmax << endl;
  std::vector<Legendre> lvec;
  make_legendre_series(ellmax, &lvec);
  for (int i=0; i<lvec.size(); ++i) {
    cout << i << ":   " << lvec[i] << endl;
    std::ostringstream ostr;
    ostr << "P" << i << ".dat" << std::ends;
    std::ofstream ofs(ostr.str().c_str());
    for (double x = -1.0; x < 1.0; x += 0.1)
      ofs << x << "   " << lvec[i](x) << endl;
  }
  // Stuff adapted from Dave L's kirkwood_test
  {
    float eps_in = 8.0;
    float eps_ext = 1.0;
    double r_b = 1.0;
    double r_a = 2.0;

    Coord origin(0,0,0);
    Atom a;
    a.coord = origin;
    a.rad = 1.0;
    a.charge = 1.0;
    a.atname = "NA";
    a.resnum = 1;
    AtomSet aset;
    aset.insert(a);
    ChargeDist arho(new AtomChargeSet(aset));

    ManyPointCharge* mrhop = new ManyPointCharge();
    mrhop->push_back(PointCharge(Coord(0,0,0), 1.0));
    ChargeDist mrho(mrhop);

    ElectrolyteEnvironment ely(new ElySphere()); //(0.0, origin, r_a));
    DielectricEnvironment spheps(new DielectricSphere
				 (eps_in, eps_ext, r_b , origin));
    ElstatPot phi_sph(spheps, mrho, ely);
    phi_sph.solve();

    const double h = 1.0e-5; // used for checking boundary conditions
    const double epsilon = 1.0e-5; // used for checking boundary conditions
    Coord p1a(0,0,r_b - epsilon - h/2.0);
    Coord p1b(0,0,r_b - epsilon + h/2.0);
    double E1 =  -(phi_sph.value(p1b) -   phi_sph.value(p1a)) /h;
    double D1 = eps_in *  E1;
    Coord p2a(0,0,r_b + epsilon - h/2.0);
    Coord p2b(0,0,r_b + epsilon + h/2.0);
    double E2 = -(phi_sph.value(p2b) -   phi_sph.value(p2a)) /h;
    double D2 = eps_ext *  E2;

    phi_sph.displacement(p1b);
    phi_sph.displacement(p2a);
    cout << "E1 = "<< E1 << endl;
    cout << "E2 = "<< E2 << endl;
    cout << "D1 = "<< D1 << endl;
    cout << "D2 = "<< D2 << endl;
  }
}


