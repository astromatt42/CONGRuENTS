#ifndef MSP_funcs_h
#define MSP_funcs_h

#include <stdio.h>
#include <math.h>

#include <gsl_sf_gamma.h>

#include "physical_constants.h"


double L_prompt__ergsm1( double M_star__Msol )
{
    double nu_sdp = 0.1;
    double norm_SDP_ergsm1Msolm1 = pow(10., 28.5);
    return nu_sdp * norm_SDP_ergsm1Msolm1 * M_star__Msol;
}

double dNdE_MSP_prompt___GeVm1( double E__GeV, double L_prompt__ergsm1, )
{
    double alpha_prompt = 2.0;
    double E_cut_prompt__GeV = 1.e0;
    double E_min__GeV = 1.e-2;
    
    return L_prompt__ergsm1 * erg__GeV/pow(E_cut_prompt__GeV, 2) * pow( E__GeV/E_cut_prompt__GeV, -1.*alpha ) * exp( -1.* E__GeV/E_cut_prompt__GeV )/gsl_sf_gamma_inc_Q( 2.+alpha, E_min__GeV/E_cut_prompt__GeV );
}



#endif
