/**
 *  @file    vpmg-private.h
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @brief   Class Vpmg private method declaration
 *  @version $Id: vpmg-private.h 1251 2008-03-30 22:32:00Z sobolevnrm $
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2008, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2008.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _VPMG_PRIVATE_H_
#define _VPMG_PRIVATE_H_

#include "apbscfg.h"
#include "apbs/vpmg.h"

/* ///////////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////////// */

/** 
 * @brief  Evaluate a cubic B-spline 
 * @author  Nathan Baker
 * @return  Cubic B-spline value
 */
VPRIVATE double bspline2(
        double x  /** Position */
        );

/** 
 * @brief  Evaluate a cubic B-spline derivative 
 * @author  Nathan Baker
 * @return  Cubic B-spline derivative
 */
VPRIVATE double dbspline2(
        double x  /** Position */
        );

/**
 * @brief   Return 2.5 plus difference of i - f
 * @author  Michael Schnieders
 * @return  (2.5+((double)(i)-(f)))
 */
VPRIVATE double VFCHI4(
		int i,
		double f
		);

/**
 * @brief   Evaluate a 5th Order B-Spline (4th order polynomial)
 * @author: Michael Schnieders
 * @return  5th Order B-Spline
 */
VPRIVATE double bspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate a 5th Order B-Spline derivative (4th order polynomial)
 * @author: Michael Schnieders
 * @return  5th Order B-Spline derivative
 */
VPRIVATE double dbspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate the 2nd derivative of a 5th Order B-Spline 
 * @author: Michael Schnieders
 * @return  2nd derivative of a 5th Order B-Spline
 */
VPRIVATE double d2bspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate the 3rd derivative of a 5th Order B-Spline 
 * @author: Michael Schnieders
 * @return  3rd derivative of a 5th Order B-Spline
 */
VPRIVATE double d3bspline4(
         double x /** Position */
         );

/**
 * @brief  Determines energy from polarizeable charge and interaction with 
 *         fixed charges according to Rocchia et al.
 * @author  Nathan Baker
 * @return  Energy in kT
 */           
VPRIVATE double Vpmg_polarizEnergy(
         Vpmg *thee, 
         int extFlag  /** If 1, add external energy contributions to 
                       result */
         ); 
/** 
 * @brief  Calculates charge-potential energy using summation over delta
 *         function positions (i.e. something like an Linf norm)
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyPoint(
        Vpmg *thee, 
        int extFlag  /** If 1, add external energy contributions to 
                       result */
        );

/** 
 * @brief  Calculates charge-potential energy as integral over a volume
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyVolume(
        Vpmg *thee, 
        int extFlag  /** If 1, add external energy contributions to 
                       result */
        );

/**
* @brief Selects a spline based surface method from either VSM_SPLINE,
 *        VSM_SPLINE5 or VSM_SPLINE7
 * @author David Gohara
 */
VPRIVATE void Vpmg_splineSelect(
		int srfm,		/** Surface method, currently VSM_SPLINE, 
		VSM_SPLINE5, or VSM_SPLINE7 */
		Vacc *acc,		/** Accessibility object */
		double *gpos,	/** Position array -> array[3] */
		double win,		/** Spline window */
		double infrad,	/** Inflation radius */
		Vatom *atom,	/** Atom object */
		double *force	/** Force array -> array[3] */
		);

/**
 * @brief  For focusing, fill in the boundaries of the new mesh based on the
 * potential values in the old mesh
 * @author  Nathan Baker
 */
VPRIVATE void focusFillBound(
        Vpmg *thee,  /** New PMG object (the one just created) */
        Vpmg *pmg  /** Old PMG object */
        );

/**
 * @brief  Increment all boundary points by
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) to add the
 *         effect of the Debye-Huckel potential due to a single charge
 * @author  Nathan Baker
 */
VPRIVATE void bcfl1(
        double size,  /** Size of the ion */
        double *apos,  /** Position of the ion */
        double charge,  /** Charge of the ion */
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */
        double *gxcf,  /** Set to x-boundary values */
        double *gycf,  /** Set to y-boundary values */
        double *gzcf,  /** Set to z-boundary values */
        double *xf,  /** Boundary point x-coordinates */
        double *yf,  /** Boundary point y-coordinates */
        double *zf,  /** Boundary point z-coordinates */
        int nx,  /** Number of grid points in x-direction */
        int ny,  /** Number of grid points in y-direction */
        int nz /** Number of grid points in y-direction */
        );

/**
 * @brief  Increment all boundary points to include the Debye-Huckel 
 *         potential due to a single multipole site. (truncated at quadrupole)
 * @author Michael Schnieders
 */
VPRIVATE void bcfl2(
        double size,  /** Size of the ion */
        double *apos,  /** Position of the ion */
        double charge,  /** Charge of the ion */
        double *dipole, /** Dipole of the ion */
        double *quad,   /** Traceless Quadrupole of the ion */
        double xkappa,  /** Exponential screening factor */
        double eps_p,   /** Solute dielectric */
        double eps_w,   /** Solvent dielectric */
        double T,       /** Temperature */
        double *gxcf,  /** Set to x-boundary values */
        double *gycf,  /** Set to y-boundary values */
        double *gzcf,  /** Set to z-boundary values */
        double *xf,  /** Boundary point x-coordinates */
        double *yf,  /** Boundary point y-coordinates */
        double *zf,  /** Boundary point z-coordinates */
        int nx,  /** Number of grid points in x-direction */
        int ny,  /** Number of grid points in y-direction */
        int nz /** Number of grid points in y-direction */
        );  

/**
 * @brief  This routine serves bcfl2. It returns (in tsr) the contraction 
 *         independent portion of the Debye-Huckel potential tensor 
 *         for a spherical ion with a central charge, dipole and quadrupole.
 *         See the code for an in depth description. 
 * 
 * @author Michael Schnieders
 */
VPRIVATE void multipolebc(
        double r,      /** Distance to the boundary */
        double kappa,  /** Exponential screening factor */ 
        double eps_p,  /** Solute dielectric */
        double eps_w,  /** Solvent dielectric */
        double rad,    /** Radius of the sphere */
        double tsr[3]  /** Contraction-independent portion of each tensor */
        );

/**
 * @brief  Calculate
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) due to a
 *         specific ion at a specific point
 * @author  Nathan Baker
 * @returns  Value of above function in arbitrary units (dependent on
 *           pre-factor)
 */
VPRIVATE double bcfl1sp(
        double size,  /** Atom size */
        double *apos,  /** Atom position */
        double charge,  /** Atom charge */ 
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */  
        double *pos  /** Function evaluation position */
        );

/** 
 * @brief  Fill boundary condition arrays
 * @author  Nathan Baker
 */
VPRIVATE void bcCalc(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill all operator coefficient arrays 
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoef(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from pre-calculated maps
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMap(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from a molecular surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMol(
        Vpmg *thee
        );

/** 
 * @brief  Fill ion (nonlinear) operator coefficient array from a molecular
 * surface calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolIon(
        Vpmg *thee
        );

/** 
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolDiel(
        Vpmg *thee
        );

/** 
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation without smoothing
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolDielNoSmooth(
        Vpmg *thee
        );

/** 
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation with smoothing.
 *
 *         Molecular surface, dielectric smoothing following an implementation
 *         of Bruccoleri, et al.  J Comput Chem 18 268-276 (1997).
 *
 *         This algorithm uses a 9 point harmonic smoothing technique - the point 
 *         in question and all grid points 1/sqrt(2) grid spacings away.
 *
 * @note   This uses thee->a1cf, thee->a2cf, thee->a3cf as temporary storage.
 * @author  Todd Dolinsky
 */
VPRIVATE void fillcoCoefMolDielSmooth(
        Vpmg *thee
        );

/** 
 * @brief  Fill operator coefficient arrays from a spline-based surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefSpline(
        Vpmg *thee
        );

/** 
* @brief  Fill operator coefficient arrays from a 5th order polynomial 
*         based surface calculation
* @author  Michael Schnieders 
*/
VPRIVATE void fillcoCoefSpline3(
		Vpmg *thee
		);

/** 
 * @brief  Fill operator coefficient arrays from a 7th order polynomial 
 *         based surface calculation
 * @author  Michael Schnieders 
 */
VPRIVATE void fillcoCoefSpline4(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill source term charge array
 * @returns  Success/failure status
 * @author  Nathan Baker
 */
VPRIVATE Vrc_Codes fillcoCharge(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from a pre-calculated map
 * @returns  Success/failure status
 * @author  Nathan Baker
 */
VPRIVATE Vrc_Codes fillcoChargeMap(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from linear interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline1(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from cubic spline interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline2(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array for the use of permanent multipoles
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoPermanentMultipole(
        Vpmg *thee
        );
        
/**
 * @brief  Fill source term charge array for use of induced dipoles 
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoInducedDipole(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array for non-local induced
 * dipoles
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoNLInducedDipole(
        Vpmg *thee
        );

/**
 * @brief  For focusing, set external energy data members in new Vpmg object
 *         based on energy calculations on old Vpmg object from regions
 *         outside the indicated partition.
 * @author  Nathan Baker, Todd Dolinsky
 */
VPRIVATE void extEnergy(
        Vpmg *thee,  /** Newly created PMG manager */
        Vpmg *pmgOLD,  /** Old PMG manager, source of energies */ 
        PBEparm_calcEnergy extFlag,  /** Energy calculation flag */
        double partMin[3],  /** Partition lower corner */
        double partMax[3],  /** Partition upper corner */ 
        int bflags[6]  /** Which boundaries to include the calculation */
        );

/** 
 * @brief  Charge-field force due to a linear spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline1(
        Vpmg *thee, 
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

/** 
 * @brief  Charge-field force due to a cubic spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline2(
        Vpmg *thee, 
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

/** 
* @brief  Charge-field force due to a quintic spline charge function
* @author  Michael Schnieders
*/
VPRIVATE void qfForceSpline4(
							 Vpmg *thee, 
							 double *force,  /** Set to force */
							 int atomID  /** Valist atom ID */
							 );


/**
 * @brief  Calculate the solution to Poisson's equation with a simple
 *         Laplacian operator and zero-valued Dirichlet boundary conditions.
 *         Store the solution in thee->u.
 * @author  Nathan Baker
 * @note  Vpmg_fillco must be called first
 */
VPRIVATE void zlapSolve(
        Vpmg *thee,
        double **solution,  /** Solution term vector */
        double **source,  /** Source term vector */
        double **work1  /** Work vector */
        );

/** 
 * @brief  Mark the grid points inside a sphere with a particular value.  This
 *         marks by resetting the the grid points inside the sphere to the
 *         specified value.
 * @author  Nathan Baker
 */
VPRIVATE void markSphere(
        double rtot,  /** Sphere radius */
        double *tpos,  /** Sphere position */
        int nx,  /** Number of grid points */
        int ny,  /** Number of grid points */
        int nz,  /** Number of grid points */
        double hx,  /** Grid spacing */
        double hy,  /** Grid spacing */
        double hzed,  /** Grid spacing */
        double xmin,  /** Grid lower corner */
        double ymin,  /** Grid lower corner */
        double zmin,  /** Grid lower corner */
        double *array,  /** Grid values */
        double markVal  /** Value to mark with */
        );

/**
 * @brief Vpmg_qmEnergy for SMPBE
 * @author Vincent Chu
 */
VPRIVATE double Vpmg_qmEnergySMPBE(Vpmg *thee, int extFlag);
VPRIVATE double Vpmg_qmEnergyNONLIN(Vpmg *thee, int extFlag);

/* Added by Vincent Chu 9/13/06 for SMPB */
#define VCUB(x)            ((x)*(x)*(x))
#define VLOG(x)            (log(x))

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define VFCHI(iint,iflt) (1.5+((double)(iint)-(iflt)))

/* ///////////////////////////////////////////////////////////////////////////
// External FORTRAN routines 
/////////////////////////////////////////////////////////////////////////// */
#define F77BCOLCOMP VF77_MANGLE(bcolcomp, BCOLCOMP)
VEXTERNC void F77BCOLCOMP(int *iparm, double *rparm, int *iwork, 
  double *rwork, double *nzval, int *rowind, int *colptr, int *flag);

#define F77PCOLCOMP VF77_MANGLE(pcolcomp, PCOLCOMP)
VEXTERNC void F77PCOLCOMP(int *nrow, int *ncol, int *nonz, 
  double *nzval, int *rowind, int *colptr, 
  char *path, char *title, char *mxtype);

#define F77MGSZ VF77_MANGLE(mgsz, MGSZ)
VEXTERNC void F77MGSZ(int *mgcoar, int *mgdisc, int *mgsolv, int *nx, int *ny,
  int *nz, int *nlev, int *nxc, int *nyc, int *nyz, int *nf, int *nc, 
  int *narr, int *narrc, int *n_rpc, int *n_iz, int *n_ipc, int *nrwk, 
  int *niwk);

#define F77MGSZAQUA VF77_MANGLE(mgszaqua, MGSZAQUA)
VEXTERNC void F77MGSZAQUA(int *mgcoar, int *mgdisc, int *mgsolv, int *nx, int *ny,
					  int *nz, int *nlev, int *nxc, int *nyc, int *nyz, int *nf, int *nc, 
					  int *narr, int *narrc, int *n_rpc, int *n_iz, int *n_ipc, int *nrwk, 
					  int *niwk);

#define F77PACKMG VF77_MANGLE(packmg, PACKMG)
VEXTERNC void F77PACKMG(int *iparm, double *rparm, int *nrwk, int *niwk,
  int *nx, int *ny, int *nz, int *nlev, int *nu1, int *nu2, int *mgkey, 
  int *itmax, int *istop, int *ipcon, int *nonlin, int *mgsmoo, int *mgprol, 
  int *mgcoar, int *mgsolv, int *mgdisc, int *iinfo, double *errtol,
  int *ipkey, double *omegal, double *omegan, int *irite, int *iperf);

#define F77CGMGDRIV VF77_MANGLE(cgmgdriv, CGMGDRIV)
VEXTERNC void F77CGMGDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NEWDRIV    VF77_MANGLE(newdriv, NEWDRIV)
VEXTERNC void F77NEWDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77MGDRIV     VF77_MANGLE(mgdriv, MGDRIV)
VEXTERNC void F77MGDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NCGHSDRIV  VF77_MANGLE(ncghsdriv, NCGHSDRIV)
VEXTERNC void F77NCGHSDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NSORDRIV   VF77_MANGLE(nsordriv, NSORDRIV)
VEXTERNC void F77NSORDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NGSRBDRIV  VF77_MANGLE(ngsrbdriv, NGSRBDRIV)
VEXTERNC void F77NGSRBDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NWJACDRIV  VF77_MANGLE(nwjacdriv, NWJACDRIV)
VEXTERNC void F77NWJACDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

#define F77NRICHDRIV  VF77_MANGLE(nrichdriv, NRICHDRIV)
VEXTERNC void F77NRICHDRIV(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf, double *tcf);

/* TEMPORARY USEAQUA */
#define F77CGMGDRIVAQUA VF77_MANGLE(cgmgdrivaqua, CGMGDRIVAQUA)
VEXTERNC void F77CGMGDRIVAQUA(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf);

/* TEMPORARY USEAQUA */
#define F77NEWDRIVAQUA    VF77_MANGLE(newdrivaqua, NEWDRIVAQUA)
VEXTERNC void F77NEWDRIVAQUA(int *iparm, double *rparm, int *iwork, double *rwork,
  double *u, double *xf, double *yf, double *zf, double *gxcf, double *gycf,
  double *gzcf, double *a1cf, double *a2cf, double *a3cf, double *ccf,
  double *fcf);

#define F77TSECND     VF77_MANGLE(tsecnd, TSECND)
#define F77VPMGANORM  VF77_MANGLE(vpmganorm, VPMGANORM)
#define F77VPMGABAND  VF77_MANGLE(vpmgaband, VPMGABAND)
#define F77DPBFA      VF77_MANGLE(dpbfa, DPBFA)
#define F77DPBDI      VF77_MANGLE(dpbdi, DPBDI)
#define F77EIGDRIV    VF77_MANGLE(eigdriv, EIGDRIV)
#define F77ANORMDRIV  VF77_MANGLE(anormdriv, ANORMDRIV)

#define F77MYPDEFINITLPBE VF77_MANGLE(mypdefinitlpbe, MYPDEFINITLPBE)
VEXTERNC void F77MYPDEFINITLPBE(int *nion, double *ionQ, double *ionConc);

#define F77MYPDEFINITNPBE VF77_MANGLE(mypdefinitnpbe, MYPDEFINITNPBE)
VEXTERNC void F77MYPDEFINITNPBE(int *nion, double *ionQ, double *ionConc);

#define F77MYPDEFINITSMPBE VF77_MANGLE(mypdefinitsmpbe, MYPDEFINITSMPBE)
VEXTERNC void F77MYPDEFINITSMPBE(int *nion, double *ionQ, double *ionConc,
									double *smvolume,double *smsize);

#define F77MYPDEFCLEAR VF77_MANGLE(mypdefclear, MYPDEFCLEAR)
VEXTERNC void F77MYPDEFCLEAR();

#endif

