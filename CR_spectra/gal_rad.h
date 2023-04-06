#ifndef gal_rad_h
#define gal_rad_h

#include <math.h>
#include "gsl_decs.h"
#include "file_io.h"

/**********************************************************************************************************************************/
/**********************************************************************************************************************************/

//Magnelli 2018: Tdust - sSFR relation - Eqn. 5
double Tdust__K( double z, double SFR__Msolyrm1, double M_star__Msol )
{
    return 98.*pow(1.+z,-0.065) + 6.9*log10(SFR__Msolyrm1/M_star__Msol);
}

/**********************************************************************************************************************************/
/**********************************************************************************************************************************/

/*
 * These are fits to Dustpedia galaxies, data taken from http://dustpedia.astro.noa.gr
 */

//Use this for FIR component, reradiated and get dust Temperature from Magnelli
double Labs__Lsol( double SFR__Msolyrm1 )
  {
  return pow( 10., 1.096548 * log10(SFR__Msolyrm1/0.56) + 9.710084 );
  }

double LFIR__Lsol( double SFR__Msolyrm1 )
  {
  return Labs__Lsol( SFR__Msolyrm1 );
  }

//Use this for low T (3500K) merged Draine field
double Lobs_old__Lsol( double M_star__Msol )
  {
  return pow( 10., 0.8480565 * log10(M_star__Msol/0.56) + 1.521623 );
  }

double L3000K__Lsol( double M_star__Msol )
  {
  return 0.574 * Lobs_old__Lsol( M_star__Msol );
  }

double L4000K__Lsol( double M_star__Msol )
  {
  return 0.426 * Lobs_old__Lsol( M_star__Msol );
  }

//Use this for high T (7500K) Draine field
double Lobs_young__Lsol( double SFR__Msolyrm1 )
  {
  if (log10( SFR__Msolyrm1 ) > -2.6 ){return pow( 10., 0.7969616 * log10(SFR__Msolyrm1/0.56) + 9.007323 );}
  else {return pow( 10., 0.9867204 * log10(SFR__Msolyrm1/0.56) + 9.476592 );}
  }

//Use this for high T (7500K) Draine field scaled down
double L7500K__Lsol( double SFR__Msolyrm1 )
  {
  return 0.763 * Lobs_young__Lsol( SFR__Msolyrm1 );
  }

//Use this for the UV Draine field scaled down
double LUV__Lsol( double SFR__Msolyrm1 )
  {
  return 0.237 * Lobs_young__Lsol( SFR__Msolyrm1 );
  }

/**********************************************************************************************************************************/
/**********************************************************************************************************************************/

/*
 * Planck's law number density
 * params[0] is T__K
 */
double dndEphot_BB__cmm3GeVm1( double *params, double E_phot__GeV )
{
    return (8.*M_PI*pow(E_phot__GeV,2))/pow(h__GeVs * c__cmsm1, 3) * 1./( exp( E_phot__GeV/(k_B__GeVKm1 * params[0]) ) - 1. );
}

/*
 * modified BB/Planck's law number density: Persic 2008, Yun & Cailli 2002 - FIR
 */
double dndEphot_modBB__cmm3GeVm1( double *params, double E_phot__GeV )
{
//    return dndEphot_BB__cmm3GeVm1( params, E_phot__GeV );
    double E_0 = 2e12 * h__GeVs;
//    double sigma = 1.;
    return dndEphot_BB__cmm3GeVm1( params, E_phot__GeV ) * E_phot__GeV/E_0;
//    return dndEphot_BB__cmm3GeVm1( params, E_phot__GeV ) * pow( E_phot__GeV/E_0, sigma );
}

/*
 * Draine: Mattis field normalised to 7500K BB field (W3)
 */
double dndEphot_UVMattis__cmm3GeVm1(  double *params, double E_phot__GeV )
{
    double lambda_micron = c__cmsm1 * h__GeVs/E_phot__GeV * 1.e4;
    if ( lambda_micron > 0.1340 && lambda_micron <= 0.2460 )
    {
        return 2.373/pow(E_phot__GeV,2) * erg__GeV * pow( lambda_micron, -0.6678 );
    }
    else if ( lambda_micron > 0.1100 && lambda_micron <= 0.1340 )
    {
        return 6.825e1/pow(E_phot__GeV,2) * erg__GeV * lambda_micron;
    }
    else if ( lambda_micron > 0.0912 && lambda_micron <= 0.1100 )
    {
        return 1.287e5/pow(E_phot__GeV,2) * erg__GeV * pow( lambda_micron, 4.4172 );
    }
    else
    {
        return 0.;
    }
}












double dndE_BB__cmm3GeVm1( double T__K, double E_phot__GeV )
{
    return (8.*M_PI*pow(E_phot__GeV,2))/pow(h__GeVs * c__cmsm1, 3) * 1./( exp( E_phot__GeV/(k_B__GeVKm1 * T__K) ) - 1. );
}


double dndE_modBB__cmm3GeVm1( double T__K, double E_phot__GeV )
{
//    return dndE_BB__cmm3GeVm1( T__K, E_phot__GeV );
    double E_0__GeV = 2e12 * h__GeVs;
//    double sigma = 1.;
    return dndE_BB__cmm3GeVm1( T__K, E_phot__GeV ) * E_phot__GeV/E_0__GeV;
//    return dndE_BB__cmm3GeVm1( T__K, E_phot__GeV ) * pow( E_phot__GeV/E_0__GeV, sigma );
}


/*
 * Radiation energy density of a BB 
 */
double u_rad_BB__GeVcmm3( double T__K )
{
    return arad__GeVcmm3Km4 * pow( T__K, 4 );
}

/* 
 * Radiation energy density of modified BB
 * This currently assumes sigma = 1, check notebook for full soln of integral need polylog and gamma functions
 */
double u_rad_modBB__GeVcmm3( double T__K )
{
//    return u_rad_BB__GeVcmm3( T__K );
    double E_0__GeV = 2e12 * h__GeVs;
    return 24.8863 * 8. * M_PI/( pow( h__GeVs * c__cmsm1, 3 ) * E_0__GeV ) * pow( k_B__GeVKm1 * T__K, 5 );
}

double u_rad_UVMattis__GeVcmm3( void )
{
    return 4450.1668;
}



//Peaks in energy of black body spectra

double E_BB_peak__GeV( double T__K )
{
    return 1.59362 * k_B__GeVKm1 * T__K;
}

double E_modBB_peak_GeV( double T__K )
{
//    return E_BB_peak__GeV( T__K );
    return 2.82144 * k_B__GeVKm1 * T__K;
}








/* Determine the dilution factor */
double C_dil( double u_rad__GeVcmm3, double L_obs__Lsol, double re__kpc, double h__pc )
  {
  return L_obs__Lsol * Lsol__GeVsm1/( u_rad__GeVcmm3 * 4. * M_PI * pow(re__kpc * 1e3 * pc__cm,2) * c__cmsm1 );
  }




/*
 * diluted black body params[ T__K, u_rad__GeVcmm3, L_obs__Lsol, re__kpc, h__pc ]
 */

/*
double dndEphot_dBB__cmm3GeVm1( double *params, double E_phot__GeV )
  {
//  return C_dil( u_rad_BB__GeVcmm3( params[0] ), params[1], params[2], params[3] ) * dndEphot_BB__cmm3GeVm1( params[0], E_phot__GeV );
  return dndEphot_BB__cmm3GeVm1( params, E_phot__GeV );
  }

double dndEphot_dmodBB__cmm3GeVm1( double *params, double E_phot__GeV )
  {
//  return C_dil( u_rad_modBB__GeVcmm3( params[0] ), params[1], params[2], params[3] ) * dndEphot_modBB__cmm3GeVm1(  params[0] , E_phot__GeV );
  return dndEphot_modBB__cmm3GeVm1( params , E_phot__GeV );
  }

double dndEphot_dUVMattis__cmm3GeVm1( double *params, double E_phot__GeV )
  {
//  return C_dil( u_rad_BB__GeVcmm3( params[0] ), params[1], params[2], params[3] ) * dndEphot_UVMattis__cmm3GeVm1( E_phot__GeV );
  return dndEphot_UVMattis__cmm3GeVm1( params, E_phot__GeV );
  }
*/

/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/*
 * To improve computational efficiency, we combine the 3000K, 4000K and 7500K black bodies into a single function, we can use one
 * interpolation object for these, integrating over the minimum of energy of the 3000K BB and the max of 7500K
 * Likewise, we combine the CMB black body and the FIR black body which vary from galaxy to galaxy
 * ACTUALLY WE CANT DO THIS, BECAUSE Cdil DEPENDS ON EACH GALAXY
 */

// params 0 : double T_CMB__K, 1 : double T_dust__K, 2 : double Mstar__Msol, 3 : double SFR_Msolyrm1, 4 : double re_light_kpc, 5 : double h__pc


int Cdils( double *params, double Cdils[5] )
{
    //FIR field from dust
    Cdils[0] = C_dil( u_rad_modBB__GeVcmm3( params[1] ),  LFIR__Lsol( params[3] ), params[4], params[5] );
    //Draine 3000K
    Cdils[1] = C_dil( u_rad_BB__GeVcmm3( 3000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    //Draine 4000K
    Cdils[2] = C_dil( u_rad_BB__GeVcmm3( 4000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    //Draine 7500K
    Cdils[3] = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( params[3] ), params[4], params[5] );
    //Draine UV Mattis field
    Cdils[4] = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( params[3] ), params[4], params[5] );

    return 0;
}

double dndEphot_total__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_3000, Cdil_4000, Cdil_7500, Cdil_UV, Cdil_FIR;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //FIR field from dust
    Cdil_FIR = C_dil( u_rad_modBB__GeVcmm3( params[1] ),  LFIR__Lsol( params[3] ), params[4], params[5] );
    //Draine 3000K
//    Cdil_3000 = C_dil( u_rad_BB__GeVcmm3( 3000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    Cdil_3000 = C_dil( u_rad_BB__GeVcmm3( 3000. ), L3000K__Lsol( params[2] ), params[4], params[5] );
    //Draine 4000K
//    Cdil_4000 = C_dil( u_rad_BB__GeVcmm3( 4000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    Cdil_4000 = C_dil( u_rad_BB__GeVcmm3( 4000. ),  L4000K__Lsol( params[2] ), params[4], params[5] );
    //Draine 7500K
    Cdil_7500 = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( params[3] ), params[4], params[5] );
    //Draine UV Mattis field
    Cdil_UV = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( params[3] ), params[4], params[5] );

    dndEphot_total__cmm3GeVm1 = dndE_BB__cmm3GeVm1( params[0], E_phot__GeV ) + Cdil_FIR * dndE_modBB__cmm3GeVm1( params[1], E_phot__GeV ) +
    Cdil_3000 * dndE_BB__cmm3GeVm1( 3000., E_phot__GeV ) + Cdil_4000 * dndE_BB__cmm3GeVm1( 4000., E_phot__GeV ) +
    Cdil_7500 * dndE_BB__cmm3GeVm1( 7500., E_phot__GeV ) + Cdil_UV * dndEphot_UVMattis__cmm3GeVm1( nphot_params, E_phot__GeV );

    return dndEphot_total__cmm3GeVm1;
}


double dndEphot_CMB__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    nphot_params[0] = params[0];
    dndEphot_total__cmm3GeVm1 += dndEphot_BB__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}

double dndEphot_FIR__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_FIR;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //FIR field from dust
    Cdil_FIR = C_dil( u_rad_modBB__GeVcmm3( params[1] ),  LFIR__Lsol( params[3] ), params[4], params[5] );
    nphot_params[0] = params[1];
    dndEphot_total__cmm3GeVm1 += Cdil_FIR * dndEphot_modBB__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}

double dndEphot_3000__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_3000;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //Draine 3000K
    Cdil_3000 = C_dil( u_rad_BB__GeVcmm3( 3000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    nphot_params[0] = 3000.;
    dndEphot_total__cmm3GeVm1 += Cdil_3000 * dndEphot_BB__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}

double dndEphot_4000__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_4000;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //Draine 4000K
    Cdil_4000 = C_dil( u_rad_BB__GeVcmm3( 4000. ), Lobs_old__Lsol( params[2] ), params[4], params[5] );
    nphot_params[0] = 4000.;
    dndEphot_total__cmm3GeVm1 += Cdil_4000 * dndEphot_BB__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}

double dndEphot_7500__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_7500;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //Draine 7500K
    Cdil_7500 = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( params[3] ), params[4], params[5] );
    nphot_params[0] = 7500.;
    dndEphot_total__cmm3GeVm1 += Cdil_7500 * dndEphot_BB__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}

double dndEphot_UV__cmm3GeVm1( double *params, double E_phot__GeV )
{
    double Cdil_UV;
    double nphot_params[1];
    double dndEphot_total__cmm3GeVm1 = 0.;

    //Draine UV Mattis field
    Cdil_UV = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( params[3] ), params[4], params[5] );
    //nphot_params is only a dummy here
    dndEphot_total__cmm3GeVm1 += Cdil_UV * dndEphot_UVMattis__cmm3GeVm1( nphot_params, E_phot__GeV );
    return dndEphot_total__cmm3GeVm1;
}


void write_radspecs( unsigned long int n_gal, double *z, double *T_dust__K, double *M_star__Msol, double *SFR__Msolyrm1, double *Re__kpc, double *h__pc, char *outfp )
{
    unsigned long int i;

    double nphot_params[6];
    double Cdil_3000, Cdil_4000, Cdil_7500, Cdil_UV, Cdil_FIR;



    unsigned int j;
    unsigned int n_E = 1000;
    double E__GeV[n_E];

    logspace_array( n_E, 1.e-13, 2.e-8, E__GeV );

    FILE *rad_specs = fopen( string_cat(outfp, "/rad_specs.txt"), "w+" );

    fprintf( rad_specs, "dnde_rad_CMB__cmm3GeVm1 dnde_rad_FIR__cmm3GeVm1 dnde_rad_3000__cmm3GeVm1 dnde_rad_4000__cmm3GeVm1 dnde_rad_7500__cmm3GeVm1 dnde_rad_UV__cmm3GeVm1 dnde_total__cmm3GeVm1\n" );

    for (j = 0; j < n_E; j++)
    {
        fprintf( rad_specs, "%e ", E__GeV[j] );
    }
    fprintf( rad_specs, "\n" );


    for (i = 0; i < n_gal; i++)
    {

        nphot_params[0] = T_0_CMB__K * (1.+z[i]);
        nphot_params[1] = T_dust__K[i];
        nphot_params[2] = M_star__Msol[i];
        nphot_params[3] = SFR__Msolyrm1[i];
        nphot_params[4] = Re__kpc[i];
        nphot_params[5] = h__pc[i];

        //FIR field from dust
        Cdil_FIR = C_dil( u_rad_modBB__GeVcmm3( nphot_params[1] ),  LFIR__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
        //Draine 3000K
        Cdil_3000 = C_dil( u_rad_BB__GeVcmm3( 3000. ), L3000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
        //Draine 4000K
        Cdil_4000 = C_dil( u_rad_BB__GeVcmm3( 4000. ), L4000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
        //Draine 7500K
        Cdil_7500 = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
        //Draine UV Mattis field
        Cdil_UV = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );


        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", dndE_BB__cmm3GeVm1( nphot_params[0], E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", Cdil_FIR * dndE_modBB__cmm3GeVm1( nphot_params[1], E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", Cdil_3000 * dndE_BB__cmm3GeVm1( 3000., E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", Cdil_4000 * dndE_BB__cmm3GeVm1( 4000., E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", Cdil_7500 * dndE_BB__cmm3GeVm1( 7500., E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", Cdil_UV * dndEphot_UVMattis__cmm3GeVm1( nphot_params, E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );
        for (j = 0; j < n_E; j++)
        {
            fprintf( rad_specs, "%e ", dndEphot_total__cmm3GeVm1( nphot_params, E__GeV[j] ) );
        }
        fprintf( rad_specs, "\n" );

    }
    fclose(rad_specs);

    return;

}






#endif
