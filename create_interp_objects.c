#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl_interp2d.h>
#include <gsl_integration.h>
#include <gsl_sf_hyperg.h>
#include <gsl_spline.h>
#include <gsl_spline2d.h>
#include "gen_funcs.h"
#include "spectra_funcs.h"
//#include "EBL_funcs.h"
#include "math_funcs.h"
//#include "cosmo_funcs.h"
#include "physical_constants.h"
#include "gal_rad.h"
#include "inverse_Compton.h"
//#include "spec_integrate.h"
//#include "CR_steadystate_3.h"

//#include "CRe_steadystate.h"

//#include "halo_mass_funcs.h"

#include <omp.h>

#include <time.h>

#include "data_objects.h"
#include "data_calc.h"

#include "file_io.h"

double q_p_inject = 2.2;
double q_e_inject = 2.2;

double E_gam_lims__GeV[2] = { 1.e-16, 1.e8 };
double T_CR_lims__GeV[2] = { 1.e-3, 1.e8 };
double E_CRe_lims__GeV[2] = { 1.e-3 + m_e__GeV, 1.e8 + m_e__GeV };






int main( int argc, char *argv[] )
{

    unsigned long int i;


    /* Read in the gals */

    char infile[strlen(argv[1]) + 1];
    snprintf(infile, strlen(argv[1]) + 1, "%s", argv[1]);
    char datadir[strlen(argv[2]) + 1];
    snprintf(datadir, strlen(argv[2]) + 1, "%s", argv[2]);


    FILE *gals_in = fopen( infile , "r" );
    fscanf( gals_in , "%*[^\n]\n");
    unsigned long int n_gal;
    fscanf( gals_in , "%lu\n" , &n_gal);
    printf("Reading %lu galaxies from file %s...\n" , n_gal, infile);
    fscanf( gals_in , "%*[^\n]\n");

    double z[n_gal], M_star__Msol[n_gal], Re__kpc[n_gal], SFR__Msolyrm1[n_gal];
    for (i = 0; i < n_gal; i++)
    {
        fscanf( gals_in , "%le %le %le %le\n" , &z[i], &M_star__Msol[i], &Re__kpc[i], &SFR__Msolyrm1[i] ); 
    }
    fclose(gals_in);

/*################################################################################################################################*/

/* Calculate calorimetry fraction */

    double Sig_star__Msolpcm2[n_gal], Sig_SFR__Msolyrm1pcm2[n_gal], Sig_gas__Msolpcm2[n_gal];
    double A_Re__pc2[n_gal], sig_gas__kmsm1[n_gal], h__pc[n_gal], n_H__cmm3[n_gal], n__cmm3[n_gal], T_dust__K[n_gal];

    #pragma omp parallel for schedule(dynamic)
    for (i = 0; i < n_gal; i++)
    {
        A_Re__pc2[i] = M_PI * pow( Re__kpc[i] * 1e3, 2 );
        Sig_star__Msolpcm2[i] = M_star__Msol[i] / ( 2. * A_Re__pc2[i] );
        Sig_SFR__Msolyrm1pcm2[i] = SFR__Msolyrm1[i] / ( 2. * A_Re__pc2[i] );
        Sig_gas__Msolpcm2[i] = Sigma_gas_Shi_iKS__Msolpcm2( Sig_SFR__Msolyrm1pcm2[i], Sig_star__Msolpcm2[i] );
        sig_gas__kmsm1[i] = sigma_gas_Yu__kmsm1( SFR__Msolyrm1[i] );
        T_dust__K[i] = Tdust__K( z[i], SFR__Msolyrm1[i], M_star__Msol[i] );
    }

    double chi = 1e-4, M_A = 2.0;
    double mu_H = 1.4, mu_p = 1.17;

    double G_h = 4.302e-3;
    double u_LA, v_Ai;
    double L_A, B__G[n_gal], B_halo__G[n_gal];

    #pragma omp parallel for schedule(dynamic) private(u_LA, v_Ai, L_A )

    for (i = 0; i < n_gal; i++)
    {

        h__pc[i] = pow( sig_gas__kmsm1[i], 2 )/( M_PI * G_h * ( Sig_gas__Msolpcm2[i] + 
                   sig_gas__kmsm1[i]/sigma_star_Bezanson__kmsm1( M_star__Msol[i], Re__kpc[i] ) * Sig_star__Msolpcm2[i] ) );

        n_H__cmm3[i] = Sig_gas__Msolpcm2[i]/( mu_H * m_H__kg * 2. * h__pc[i] ) * Msol__kg/pow( pc__cm, 3 );

        n__cmm3[i] = n_H__cmm3[i] * mu_H/mu_p;

        u_LA = sig_gas__kmsm1[i]/sqrt(2.);
        v_Ai = 1000. * ( u_LA/10. )/( sqrt(chi/1e-4) * M_A );
        L_A = h__pc[i]/pow(M_A,3);

        B__G[i] = sqrt(4.* M_PI * chi * n__cmm3[i] * mu_p * m_H__kg * 1e3) * v_Ai * 1e5;

        if (log10(SFR__Msolyrm1[i]/M_star__Msol[i]) > -10.)
        {
            B_halo__G[i] = B__G[i]/3.;
        }
        else
        {
            B_halo__G[i] = B__G[i]/1.5;
        }

    }


/*################################################################################################################################*/
/*################################################################################################################################*/


    double E_phot_lims__GeV[2] = { E_BB_peak__GeV( T_0_CMB__K ) * 1e-4, 1.e-7 };


    size_t n_pts[2] = { 1000, 1000 };

    double x_SY_lims[2];

    x_SY_lims[0] = (2.*pow(M_PI,2)*pow(m_e__g,2)*pow(c__cmsm1,3)) / (3.*e__esu*maxval( n_gal, B__G )*h__ergs) * 
                   (E_gam_lims__GeV[0] * m_e__GeV)/(pow(E_CRe_lims__GeV[1],2)) * 0.1;

    x_SY_lims[1] = (2.*pow(M_PI,2)*pow(m_e__g,2)*pow(c__cmsm1,3)) / (3.*e__esu*minval( n_gal, B_halo__G )/10.*h__ergs) * 
                   (E_gam_lims__GeV[1] * (1+maxval( n_gal, z )) * m_e__GeV)/(pow(E_CRe_lims__GeV[0],2)) * 10.;


    IC_object ICo = load_IC_do_files( n_pts, E_gam_lims__GeV, E_CRe_lims__GeV, E_phot_lims__GeV,
                                      (maxval( n_gal, z ) + 1.) * T_0_CMB__K, 0.5, 
                                      minval( n_gal, T_dust__K ), maxval( n_gal, T_dust__K ), 5., datadir );


    data_object_2D do_2D_BS = load_BS_do_files( n_pts, E_gam_lims__GeV, E_CRe_lims__GeV, datadir );
    gsl_spline_object_2D gso2D_BS = do2D_to_gso2D( do_2D_BS );

    data_object_1D do_1D_SY = load_SY_do_files( &(n_pts[0]), x_SY_lims, datadir );
    gsl_spline_object_1D gso1D_SY = do1D_to_gso1D( do_1D_SY );

    printf("Finished creating interpolation objects!\n");

    return 0.;

}
