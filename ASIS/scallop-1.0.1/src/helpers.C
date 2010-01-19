#include "helpers.h"

void refine(P_DOMAIN& fine, const P_DOMAIN& coarse, int nref) {
    fine.resize(coarse.size());
    for_1(i, coarse) {
	fine.setowner(i, coarse(i).owner());
    } end_for;
    POINT lo;
    POINT hi;
    for_1(i, coarse) {
	lo = coarse(i).lower();
	hi = coarse(i).upper();
	lo *= nref;
	hi = (hi+1)*nref - 1;
	fine.setlower(i, lo);
	fine.setupper(i, hi);
    } end_for;
}

void coarsen(P_DOMAIN& coarse, const P_DOMAIN& fine, int nref) {
    coarse.resize(fine.size());
    for_1(i, fine) {
	coarse.setowner(i, fine(i).owner());
    } end_for;
    POINT lo;
    POINT hi;
    for_1(i, fine) {
	lo = fine(i).lower();
	hi = fine(i).upper();
	lo(0) = ( (lo(0) < 0) ? ( (lo(0) + 1)/nref - 1 ) : lo(0)/nref );
	lo(1) = ( (lo(1) < 0) ? ( (lo(1) + 1)/nref - 1 ) : lo(1)/nref );
	lo(2) = ( (lo(2) < 0) ? ( (lo(2) + 1)/nref - 1 ) : lo(2)/nref );
	hi(0) = ( (hi(0) < 0) ? ( (hi(0) + 1)/nref - 1 ) : hi(0)/nref );
	hi(1) = ( (hi(1) < 0) ? ( (hi(1) + 1)/nref - 1 ) : hi(1)/nref );
	hi(2) = ( (hi(2) < 0) ? ( (hi(2) + 1)/nref - 1 ) : hi(2)/nref );
	coarse.setlower(i, lo);
	coarse.setupper(i, hi);
    } end_for;
}

P_DOMAIN coarsen(const P_DOMAIN& fine, int nref) {
    P_DOMAIN tmp(fine);
    coarsen(tmp, fine, nref);
    return P_DOMAIN(tmp);
}

       
P_DOMAIN ncoarsen(const P_DOMAIN& fine, int nref) {
    P_DOMAIN tmp(fine);
    POINT lo;
    POINT hi;
    for_1(i, fine) {
	lo = fine(i).lower();
	hi = fine(i).upper();
	lo(0) = ( (lo(0) < 0) ? ( (lo(0) + 1)/nref - 1) : lo(0)/nref );
	lo(1) = ( (lo(1) < 0) ? ( (lo(1) + 1)/nref - 1) : lo(1)/nref );
	lo(2) = ( (lo(2) < 0) ? ( (lo(2) + 1)/nref - 1) : lo(2)/nref );
	hi(0) = ( (hi(0) < 0) ? hi(0)/nref : (hi(0) - 1)/nref + 1 );
	hi(1) = ( (hi(1) < 0) ? hi(1)/nref : (hi(1) - 1)/nref + 1 );
	hi(2) = ( (hi(2) < 0) ? hi(2)/nref : (hi(2) - 1)/nref + 1 );
	tmp.setlower(i, lo);
	tmp.setupper(i, hi);
    } end_for;
    return P_DOMAIN(tmp);
}

P_DOMAIN* ncoarsen(const P_DOMAIN* fine, int nref) {
    P_DOMAIN *tmp = new P_DOMAIN(*fine);
    POINT lo;
    POINT hi;
    for_1(i, *fine) {
	lo = (*fine)(i).lower();
	hi = (*fine)(i).upper();
	lo(0) = ( (lo(0) < 0) ? ( (lo(0) + 1)/nref - 1) : lo(0)/nref );
	lo(1) = ( (lo(1) < 0) ? ( (lo(1) + 1)/nref - 1) : lo(1)/nref );
	lo(2) = ( (lo(2) < 0) ? ( (lo(2) + 1)/nref - 1) : lo(2)/nref );
	hi(0) = ( (hi(0) < 0) ? hi(0)/nref : (hi(0) - 1)/nref + 1 );
	hi(1) = ( (hi(1) < 0) ? hi(1)/nref : (hi(1) - 1)/nref + 1 );
	hi(2) = ( (hi(2) < 0) ? hi(2)/nref : (hi(2) - 1)/nref + 1 );
	tmp->setlower(i, lo);
	tmp->setupper(i, hi);
    } end_for;
    return tmp;
}


P_DOMAIN *nrefine(const P_DOMAIN* coarse, int nref) {
    P_DOMAIN *tmp = new P_DOMAIN(*coarse);
    POINT lo;
    POINT hi;
    for_1(i, *coarse) {
	lo = (*coarse)(i).lower();
	hi = (*coarse)(i).upper();
	lo *= nref;
	hi *= nref;
	tmp->setlower(i, lo);
	tmp->setupper(i, hi);
    } end_for;
    return tmp;
}

P_DOMAIN nrefine(const P_DOMAIN& coarse, int nref) {
    P_DOMAIN tmp(coarse);
    POINT lo;
    POINT hi;
    for_1(i, coarse) {
	lo = coarse(i).lower();
	hi = coarse(i).upper();
	lo *= nref;
	hi *= nref;
	tmp.setlower(i, lo);
	tmp.setupper(i, hi);
    } end_for;
    return P_DOMAIN(tmp);
}

REGION coarsen(const REGION& fine, int nref) {
    REGION tmp(fine);
    POINT lo;
    POINT hi;
    lo = fine.lower();
    hi = fine.upper();
    lo(0) = ( (lo(0) < 0) ? ( (lo(0) + 1)/nref - 1 ) : lo(0)/nref );
    lo(1) = ( (lo(1) < 0) ? ( (lo(1) + 1)/nref - 1 ) : lo(1)/nref );
    lo(2) = ( (lo(2) < 0) ? ( (lo(2) + 1)/nref - 1 ) : lo(2)/nref );
    hi(0) = ( (hi(0) < 0) ? ( (hi(0) + 1)/nref - 1 ) : hi(0)/nref );
    hi(1) = ( (hi(1) < 0) ? ( (hi(1) + 1)/nref - 1 ) : hi(1)/nref );
    hi(2) = ( (hi(2) < 0) ? ( (hi(2) + 1)/nref - 1 ) : hi(2)/nref );
    tmp.setlower(lo);
    tmp.setupper(hi);
    return REGION(tmp);
}


REGION refine(const REGION& coarse, int nref) {
    REGION tmp(coarse);
    POINT lo;
    POINT hi;
    lo = coarse.lower();
    hi = coarse.upper();
    lo *= nref;
    hi = (hi+1)*nref - 1;
    tmp.setlower(lo);
    tmp.setupper(hi);
    return REGION(tmp);
}
    

REGION ncoarsen(const REGION& fine, int nref) {
    REGION tmp(fine);
    POINT lo;
    POINT hi;
    lo = fine.lower();
    hi = fine.upper();
    lo(0) = ( (lo(0) < 0) ? ( (lo(0) + 1)/nref - 1) : lo(0)/nref );
    lo(1) = ( (lo(1) < 0) ? ( (lo(1) + 1)/nref - 1) : lo(1)/nref );
    lo(2) = ( (lo(2) < 0) ? ( (lo(2) + 1)/nref - 1) : lo(2)/nref );
    hi(0) = ( (hi(0) < 0) ? hi(0)/nref : (hi(0) - 1)/nref + 1 );
    hi(1) = ( (hi(1) < 0) ? hi(1)/nref : (hi(1) - 1)/nref + 1 );
    hi(2) = ( (hi(2) < 0) ? hi(2)/nref : (hi(2) - 1)/nref + 1 );
    tmp.setlower(lo);
    tmp.setupper(hi);
    return REGION(tmp);
}

REGION nrefine(const REGION& coarse, int nref) {
    REGION tmp(coarse);
    POINT lo;
    POINT hi;
    lo = coarse.lower();
    hi = coarse.upper();
    lo *= nref;
    hi *= nref;
    tmp.setlower(lo);
    tmp.setupper(hi);
    return REGION(tmp);
}

REGION grow_lo(const REGION& r, int direction, int n) {
    REGION tmp(r);
    if (direction == 0) {
	tmp.setlower(r.lower(0) - n, r.lower(1), r.lower(2));
    }
    else if (direction == 1) {
	tmp.setlower(r.lower(0), r.lower(1) - n, r.lower(2));
    }
    else if (direction == 2) {
	tmp.setlower(r.lower(0), r.lower(1), r.lower(2) - n);
    }
    return REGION(tmp);
}


REGION grow_hi(const REGION& r, int direction, int n) {
    REGION tmp(r);
    if (direction == 0) {
	tmp.setupper(r.upper(0) + n, r.upper(1), r.upper(2));
    }
    else if (direction == 1) {
	tmp.setupper(r.upper(0), r.upper(1) + n, r.upper(2));
    }
    else if (direction == 2) {
	tmp.setupper(r.upper(0), r.upper(1), r.upper(2) + n);
    }
    return REGION(tmp);
}


REGION shift(const REGION& r, int direction, int n) {
    REGION tmp(r);
    if (direction == 0) {
	tmp.setlower(r.lower(0) + n, r.lower(1), r.lower(2));
	tmp.setupper(r.upper(0) + n, r.upper(1), r.upper(2));
    }
    else if (direction == 1) {
	tmp.setlower(r.lower(0), r.lower(1) + n, r.lower(2));
	tmp.setupper(r.upper(0), r.upper(1) + n, r.upper(2));
    }
    else if (direction == 2) {
	tmp.setlower(r.lower(0), r.lower(1), r.lower(2) + n);
	tmp.setupper(r.upper(0), r.upper(1), r.upper(2) + n);
    }
    return REGION(tmp);
}


POINT coarsen(const POINT& p, int nref) {
    POINT tmp(p);
    tmp(0) = ((tmp(0) < 0) ? (-abs(tmp(0) + 1)/nref) - 1 : tmp(0)/nref);
    tmp(1) = ((tmp(1) < 0) ? (-abs(tmp(1) + 1)/nref) - 1 : tmp(1)/nref);
    tmp(2) = ((tmp(2) < 0) ? (-abs(tmp(2) + 1)/nref) - 1 : tmp(2)/nref);
    return POINT(tmp);
}


POINT shift(const POINT& p, int direction, int n) {
    POINT tmp(p);
    if (direction == 0) {
	tmp(0) += n;
    }
    else if (direction == 1) {
	tmp(1) += n;
    }
    else if (direction == 2) {
	tmp(2) += n;
    }
    return POINT(tmp);
}
