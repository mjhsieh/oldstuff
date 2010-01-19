#include "stencil.h"

#define FORT_SMLAP F77_FUNC(smlap, SMLAP)

extern "C" {
    void FORT_SMLAP(FAA_TYPE, double*, double*);
};

stencil::stencil(P_DOMAIN *input_coarse_domain, P_DOMAIN *input_fine_domain) {
    coarse_domain = input_coarse_domain;
    fine_domain = input_fine_domain;
    coarse_data = new P_GRID(*coarse_domain);
    fine_data = new P_GRID(*fine_domain);
}

// This is just for debugging, but I need to implement the FORT_SMLAP
// routine first.
double stencil::laplacian()
{

    double Result = 0.;
//    FORT_SMLAP(FORT_ARRAY_ARG_PTR(coarse_data), &CoarseH, &Result);

    return Result;
}

void stencil::reset()
{
    coarse_data->instantiate(*coarse_domain);
    coarse_data->fill(0.);
}    
