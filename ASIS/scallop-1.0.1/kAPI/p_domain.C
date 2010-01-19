#include "p_domainNDIM.h"

p_domainNDIM* grow(const p_domainNDIM* b, int dim) {
    p_domainNDIM *tmp = new p_domainNDIM(*b);
    for_1(i, *b) 
	tmp->grow(i, dim);
//	tmp(i).RegionNDIM::grow(dim);
    end_for;
    return tmp;
}

p_domainNDIM grow(const p_domainNDIM& b, int dim) {
    p_domainNDIM tmp = p_domainNDIM(b);
    for_1(i, b) 
	tmp.grow(i, dim);
//	tmp(i).RegionNDIM::grow(dim);
    end_for;
    return p_domainNDIM(tmp);
}
