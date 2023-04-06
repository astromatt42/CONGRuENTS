#ifndef freefreeabsorption_h
#define freefreeabsorption_h

//IO for testing, don't need this after
//#include <stdio.h>
#include <math.h>
//Draine 2010
double gaunt_ff( double Z_i, double nu__Hz, double T__K )
{
return 4.691 * ( 1. - 0.118 * log( Z_i/10. * nu__Hz/1.e9 * pow( T__K/1.e4, -3./2.) ) );
}

double kap_nu__cmm1( double n_e__cmm3, double n_i__cmm3, double Z_i, double nu__Hz, double T__K )
{
return 8.31e-26 * ( 1. - 0.118 * log( Z_i/10. * nu__Hz/1.e9 * pow( T__K/1.e4, -3./2.) ) ) * n_e__cmm3 * n_i__cmm3 * pow( Z_i, 2 ) * pow( T__K/1.e4, -3./2.) * pow( nu__Hz/1.e9, -2);
}

#endif
