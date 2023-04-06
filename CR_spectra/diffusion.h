#ifndef diffusion_h
#define diffusion_h

#include <math.h>
#include "physical_constants.h"
#include "gsl_decs.h"

/*
 * simple assumption of diffusion according to MW estimate in MK2020+ for now
 * 
 */

double dEdtm1_diff__GeVsm1( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 )
{
//    double D = 3.e27 * pow(E_e__GeV, 0.5); //[cm^2 s^{-1}]
    double D = gsl_spline_eval( gso_1D_D__cm2sm1.spline, E_e__GeV, gso_1D_D__cm2sm1.acc ); //[cm^2 s^{-1}]
    return -1. * E_e__GeV * D/pow(h__pc * pc__cm,2);
}


double deldelEm1dEdtm1_diff__sm1( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 )
{
//    double D = 3.e27 * pow(E_e__GeV, 0.5); //[cm^2 s^{-1}]
    double D = gsl_spline_eval( gso_1D_D__cm2sm1.spline, E_e__GeV, gso_1D_D__cm2sm1.acc ); //[cm^2 s^{-1}]
    double deldelED = gsl_spline_eval_deriv( gso_1D_D__cm2sm1.spline, E_e__GeV, gso_1D_D__cm2sm1.acc ); //[cm^2 s^{-1} GeV^{-1}]
//    return -1. * 3./2. * D/pow(h__pc * pc__cm,2);
    return -1./pow(h__pc * pc__cm,2) * ( D + E_e__GeV * deldelED );
}


double tau_diff__s( double E_e__GeV, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1 )
{
    return -1.* E_e__GeV/dEdtm1_diff__GeVsm1( E_e__GeV, h__pc, gso_1D_D__cm2sm1 );
}


#endif
