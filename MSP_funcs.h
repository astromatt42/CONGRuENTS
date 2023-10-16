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

double d2NdEdt_MSP_prompt___GeVm1( double E__GeV, double L_prompt__ergsm1 )
{
    double alpha_prompt = 0.039;
    double E_cut_prompt__GeV = 1.e0;
    double E_min__GeV = 1.e-2;
    
    return L_prompt__ergsm1 * erg__GeV/pow(E_cut_prompt__GeV, 2) * pow( E__GeV/E_cut_prompt__GeV, alpha_prompt ) * exp( -1.* E__GeV/E_cut_prompt__GeV )/gsl_sf_gamma_inc_Q( 2.+alpha_prompt, E_min__GeV/E_cut_prompt__GeV );
}

double d2NdEdt_MSP_e___GeVm1( double E__GeV, double L_prompt__ergsm1 )
{
    double alpha_prompt = 0.039;
    double gamma_MSP = 3. * alpha_prompt + 1.;
    double E_cut_prompt__GeV = 1.e0;
    double rho_c__km = 30.;
    double E_cut_e__GeV = m_e__GeV * pow( 2. * rho_c__km * E_cut_prompt__GeV * 2. * M_PI/(3.*h__GeVs*c__cmsm1/1.e5), 1./3. );
    double E_min__GeV = 1.e-2;
    double L_IC__egsm1 = 9. * L_prompt__ergsm1;
    
    return L_IC__egsm1 * erg__GeV/pow(E_cut_e__GeV, 2) * pow( E__GeV/E_cut_e__GeV, gamma_MSP ) * exp( -1.* E__GeV/E_cut_e__GeV )/gsl_sf_gamma_inc_Q( 2.+gamma_MSP, E_min__GeV/E_cut_e__GeV );
}

#endif
