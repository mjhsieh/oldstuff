#ifndef helpers_H
#define helpers_H

#include "kelp.h"
#include "def3.h"
#include "p_domain3.h"

void refine(P_DOMAIN& fine, const P_DOMAIN& coarse, int nref);

void coarsen(P_DOMAIN& coarse, const P_DOMAIN& fine, int nref);

P_DOMAIN coarsen(const P_DOMAIN& fine, int nref);

P_DOMAIN ncoarsen(const P_DOMAIN& fine, int nref);

P_DOMAIN* ncoarsen(const P_DOMAIN* fine, int nref);

P_DOMAIN* nrefine(const P_DOMAIN* coarse, int nref);

P_DOMAIN nrefine(const P_DOMAIN& coarse, int nref);

REGION coarsen(const REGION& fine, int nref);

REGION refine(const REGION& fine, int nref);

REGION ncoarsen(const REGION& fine, int nref);

REGION nrefine(const REGION& fine, int nref);

REGION grow_lo(const REGION& r, int direction, int n);

REGION grow_hi(const REGION& r, int direction, int n);

REGION shift(const REGION& r, int direction, int n);

POINT coarsen(const POINT& p, int nref);

POINT shift(const POINT& p, int direction, int n);

#endif // ifndef helpers_H       
    
