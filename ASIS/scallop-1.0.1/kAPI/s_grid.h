#ifndef _s_gridNDIM_h
#define _s_gridNDIM_h

/*

s_grid is a single processor/serial grid, an extension of the standard
KeLP Grid, with the kinds of extras I use a lot

*/

#include "config.h"
#include "GridNDIM.h"

#define DIMNDIMD

class f_regionNDIM {
 public:
#if defined(DIM2D)
    int l0, l1, u0, u1;
    f_region2(const Region2& r): l0(r.lower(0)), l1(r.lower(1)),
	u0(r.upper(0)), u1(r.upper(1)) {}
#elif defined(DIM3D)
    int l0, l1, l2, u0, u1, u2;
    f_region3(const Region3& r): l0(r.lower(0)), l1(r.lower(1)),
	l2(r.lower(2)), u0(r.upper(0)), u1(r.upper(1)), u2(r.upper(2)) {}
#else
#error Only 2D and 3D s_grids are defined.
#endif
};    
 
class s_gridNDIM: public GridNDIM<double> {
 public:
    f_regionNDIM domain;

#if defined(DIM2D)
    s_grid2() : Grid2<double>(), domain(Region2(0, 0, 0, 0)) {}
#elif defined(DIM3D)
    s_grid3() : Grid3<double>(), domain(Region3(0, 0, 0, 0, 0, 0)) {}
#else
#error Only 2D and 3D s_grids are defined.
#endif

    s_gridNDIM(const RegionNDIM& region, const int processor) 
	: GridNDIM<double>(region, processor), domain(region) { 
	cout << region.upper(1) << endl;
    }

    void resize(const RegionNDIM& R, const int alloc=TRUE) {
	GridNDIM<double>::resize(R, alloc);
	domain = f_regionNDIM(R);
    }
    void fill(const double val, const RegionNDIM& R);
    void fill(const double val) {fill(val,this->region());	}
   void negate();
   void mult(double s);
   double norm(const int p);
};

inline void DefaultPack(s_gridNDIM& G, const RegionNDIM& R, double *buf) {  
    G.PackRegion(buf,R);
}

inline void DefaultUnpack(s_gridNDIM& G, const RegionNDIM& R, double *buf) {
    G.UnPackRegion(buf,R);
}

inline void DefaultCopy(s_gridNDIM& FG, const RegionNDIM& FR, 
			s_gridNDIM& TG, const RegionNDIM& TR) {
    TG.CopyRegion(TR,FG,FR);
}

#endif
