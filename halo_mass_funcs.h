#ifndef halo_mass_funcs_h
#define halo_mass_funcs_h

#include <math.h>
#include "gsl_decs.h"

struct halo_mass_obj
{
    gsl_spline_object_1D logM1_hm_gso1D;
    gsl_spline_object_1D logMstar0_hm_gso1D; 
    gsl_spline_object_1D beta_hm_gso1D;
    gsl_spline_object_1D delta_hm_gso1D;
    gsl_spline_object_1D gamma_hm_gso1D;
};


struct halo_mass_obj halo_mass_init( void )
{
    size_t n = 16;
    //data from Rodriguez-Puebla et al. (2017): Constraining the Galaxy-Halo Connection Over The Last 13.3 Gyrs... for values
    //parametrisation Behroozi et al. (2010): A COMPREHENSIVE ANALYSIS OF UNCERTAINTIES AFFECTING THE STELLAR MASS–HALO MASS RELATION FOR 0 < z < 4
    double z[16] = { 0., 0.10, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00  };
    double logM1[16] = { 12.56, 12.58, 12.61, 12.68, 12.77, 12.89, 13.01, 13.15, 13.33, 13.51, 14.02, 14.97, 14.86, 17.43, 17.27, 16.79 };
    double logMstar0[16] = { 10.88, 10.90, 10.93, 10.99, 11.08, 11.19, 11.31, 11.47, 11.73, 12.14, 12.73, 14.31, 14.52, 16.69, 20.24, 21.89 };
    double beta[16] = { 0.48, 0.48, 0.48, 0.48, 0.50, 0.51, 0.53, 0.54, 0.55, 0.55, 0.59, 0.60, 0.58, 0.55, 0.52, 0.43 };
    double delta[16] = { 0.30, 0.29 , 0.27, 0.23, 0.18, 0.12, 0.03, -0.10, -0.34, -0.44, -0.44, -0.44, -0.44, -0.44, -0.44, -0.44 };
    double gamma[16] = { 1.56, 1.52, 1.46, 1.39, 1.33, 1.27, 1.22, 1.17, 1.16, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92 };

    struct halo_mass_obj hm_obj;

    hm_obj.logM1_hm_gso1D = gsl_so1D( n, z, logM1 );
    hm_obj.logMstar0_hm_gso1D = gsl_so1D( n, z, logMstar0 );
    hm_obj.beta_hm_gso1D = gsl_so1D( n, z, beta );
    hm_obj.delta_hm_gso1D = gsl_so1D( n, z, delta );
    hm_obj.gamma_hm_gso1D = gsl_so1D( n, z, gamma );

    return hm_obj;
}


double halo_mass__Msol( struct halo_mass_obj hm_obj, double M_star__Msol, double z )
{
    double logM1 = gsl_so1D_eval( hm_obj.logM1_hm_gso1D, z );
    double Mstar0 = pow( 10., gsl_so1D_eval( hm_obj.logMstar0_hm_gso1D, z ) );
    double beta = gsl_so1D_eval( hm_obj.beta_hm_gso1D, z );
    double delta = gsl_so1D_eval( hm_obj.delta_hm_gso1D, z );
    double gamma = gsl_so1D_eval( hm_obj.gamma_hm_gso1D, z );
    return pow( 10., logM1 + beta * log10( M_star__Msol/Mstar0 ) + pow( M_star__Msol/Mstar0, delta )/( 1. + pow( M_star__Msol/Mstar0, -gamma ) ) - 0.5);
}


//Convert Reffective (half-light) to half mass, using result from Suess et al. (2019) extrapolated to z = 0: HALF-MASS RADII FOR ∼ 7, 000 GALAXIES AT 1.0 ≤ z ≤ 2.5: MOST OF THE EVOLUTION IN THE MASS-SIZE RELATION IS DUE TO COLOR GRADIENTS
//Below z= 1 use Suess (2019) Half-mass Radii of Quiescent and Star-forming Galaxies Evolve Slowly from: Implications for Galaxy Assembly Histories
double R_half_mass__kpc( double R_e__kpc, double z )
{
    double s;
    double b;


    if ( z < 1.)
    {
        return 0.7 * R_e__kpc;
    }
    else if ( z < 1.5 )
    {
        s = -0.325;
        b = -0.180;
    }
    else if ( z < 2.0 )
    {
        s = -0.100;
        b = -0.081;
    }
    else
    {
        s = -0.034;
        b = -0.029;
    }
    return R_e__kpc * pow(10., s * (log10(R_e__kpc) - 1.) + b);
}

//Rohr et al. (2022): The galaxy–halo size relation of low-mass galaxies in FIRE
//invert the relationship
double R_vir__kpc( double R_half_mass__kpc )
{
    return 35. * pow( 18.87 * R_half_mass__kpc, 1.07 );

}






#endif
