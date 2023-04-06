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
#include "EBL_funcs.h"
#include "math_funcs.h"
#include "cosmo_funcs.h"
#include "physical_constants.h"
#include "gal_rad.h"
#include "inverse_Compton.h"
#include "spec_integrate.h"
//#include "CR_steadystate_3.h"

#include "CRe_steadystate.h"

#include "halo_mass_funcs.h"

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



    //Initiate the halo mass routine
//    struct halo_mass_obj hm_obj = halo_mass_init();

    unsigned long int i, j, k;

    time_t seconds;






/*################################################################################################################################*/

/* Read tau EBL from file Dominguez data */

    size_t nx_E, ny_z;
//    FILE *tau_in = fopen( "input/tau_Eg_z_Dominguez.txt", "r" );
    FILE *tau_in = fopen( "input/tau_Eg_z_Franceschini.txt", "r" );
//    FILE *tau_in = fopen( "input/tau_Eg_z_Gilmore.txt", "r" );
    fscanf( tau_in , "%*[^\n]\n");
    fscanf( tau_in , "%lu \n", &ny_z );
    fscanf( tau_in , "%*[^\n]\n");
    fscanf( tau_in , "%lu \n", &nx_E );
    fscanf( tau_in , "%*[^\n]\n");

    double xa_E[nx_E];
    double ya_z[ny_z]; 
    double za_tau[nx_E * ny_z];

    /* File IO */
    /* Careful with \n characters in za_tau read-in */
    for (i=0;i<ny_z;i++){fscanf( tau_in ,"%lf", &ya_z[i] );}
    fscanf( tau_in , "\n");
    fscanf( tau_in , "%*[^\n]\n");
    for (i=0;i<nx_E;i++){fscanf( tau_in ,"%le", &xa_E[i] );}
    fscanf( tau_in , "\n");
    fscanf( tau_in , "%*[^\n]\n");
    for ( i = 0 ; i < (nx_E * ny_z) ; i++){fscanf( tau_in ,"%le", &za_tau[i] );}
    fclose(tau_in);


    /* Assign fdata for interpolation for data & error - energy in eV */
    fd_in fdata_in;
    fdata_in.nx = nx_E;
    fdata_in.ny = ny_z;
    fdata_in.xa = (double*) malloc(nx_E * sizeof(double));
    fdata_in.ya = (double*) malloc(ny_z * sizeof(double));
    fdata_in.za = (double*) malloc( (nx_E * ny_z) * sizeof(double));
    for (i=0; i < nx_E; i++) fdata_in.xa[i] = xa_E[i];
    for (i=0; i < ny_z; i++) fdata_in.ya[i] = ya_z[i];
    for (i=0; i < (ny_z * nx_E); i++) fdata_in.za[i] = za_tau[i];
    fdata_in.T = gsl_interp2d_bilinear;
    fdata_in.interp = gsl_interp2d_alloc(fdata_in.T, fdata_in.nx, fdata_in.ny);
    fdata_in.xacc = gsl_interp_accel_alloc();
    fdata_in.yacc = gsl_interp_accel_alloc();
    gsl_interp2d_init(fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, fdata_in.nx, fdata_in.ny);

/*
    fd_in fdata_in_err;
    fdata_in_err.nx = nx_E;
    fdata_in_err.ny = ny_z;
    fdata_in_err.xa = (double*) malloc(nx_E * sizeof(double));
    fdata_in_err.ya = (double*) malloc(ny_z * sizeof(double));
    fdata_in_err.za = (double*) malloc( (nx_E * ny_z) * sizeof(double));
    for (i=0; i < nx_E; i++) fdata_in_err.xa[i] = xa_E[i];
    for (i=0; i < ny_z; i++) fdata_in_err.ya[i] = ya_z[i];
    for (i=0; i < (ny_z * nx_E); i++) fdata_in_err.za[i] = za_tau_err[i];
    fdata_in_err.T = gsl_interp2d_bilinear;
    fdata_in_err.interp = gsl_interp2d_alloc(fdata_in.T, fdata_in.nx, fdata_in.ny);
    fdata_in_err.xacc = gsl_interp_accel_alloc();
    fdata_in_err.yacc = gsl_interp_accel_alloc();
    gsl_interp2d_init(fdata_in_err.interp, fdata_in_err.xa, fdata_in_err.ya, fdata_in_err.za, fdata_in_err.nx, fdata_in_err.ny);
*/

/*################################################################################################################################*/



    /* Read in the gals */

    char infile[strlen(argv[1]) + 1];
    snprintf(infile, strlen(argv[1]) + 1, "%s", argv[1]);
    char datadir[strlen(argv[2]) + 1];
    snprintf(datadir, strlen(argv[2]) + 1, "%s", argv[2]);
    char outfp[strlen(argv[3]) + 1];
    snprintf(outfp, strlen(argv[3]) + 1, "%s", argv[3]);


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

    double chi = 1e-4, M_A = 2.0, beta = 0.25;
    double sigma_pp_cm2 = 40e-27, mu_H = 1.4, mu_p = 1.17;

    /* Number of stars that go SN per solar mass of stars formed - C2003 IMF */
    double n_SN_Msolm1 = 1.321680e-2;
    /* fraction of energy that goes into CRs */
    double f_EtoCR = 0.1;

    /* fraction of CRp energy that goes into primary CRes */
    double f_CRe_CRp = 0.2;
    /* Energy of each SN */
    double E_SN_erg = 1e51;
    /* Pre-factor for the ion Alfven speed, default is 1. */
    double f_vAi = 1.;


    double T_p_cutoff__GeV = 1e8;
    double T_e_cutoff__GeV = 1e5;


    unsigned int n_T_CR = 1000;

    double T_CR__GeV[n_T_CR];
    logspace_array( n_T_CR, T_CR_lims__GeV[0], T_CR_lims__GeV[1], T_CR__GeV );
    write_1D_file( n_T_CR, T_CR__GeV, "T_CR__GeV", string_cat(outfp, "/T_CR.txt") );

    double E_CRe__GeV[n_T_CR];
    logspace_array( n_T_CR, E_CRe_lims__GeV[0], E_CRe_lims__GeV[1], E_CRe__GeV );

    int n_E_gam = 500;
    double E_gam__GeV[n_E_gam];
    logspace_array( n_E_gam, E_gam_lims__GeV[0], E_gam_lims__GeV[1], E_gam__GeV );
    write_1D_file( n_E_gam, E_gam__GeV, "E_gam__GeV", string_cat(outfp, "/E_gam.txt") );

  
    double **f_cal = malloc(sizeof *f_cal * n_gal);
    if (f_cal){for (i = 0; i < n_gal; i++){f_cal[i] = malloc(sizeof *f_cal[i] * n_T_CR);}}
    double **D__cm2sm1 = malloc(sizeof *D__cm2sm1 * n_gal);
    if (D__cm2sm1){for (i = 0; i < n_gal; i++){D__cm2sm1[i] = malloc(sizeof *D__cm2sm1[i] * n_T_CR);}}
    double **D_e__cm2sm1 = malloc(sizeof *D_e__cm2sm1 * n_gal);
    if (D_e__cm2sm1){for (i = 0; i < n_gal; i++){D_e__cm2sm1[i] = malloc(sizeof *D_e__cm2sm1[i] * n_T_CR);}}
    double **D_e_z2__cm2sm1 = malloc(sizeof *D_e_z2__cm2sm1 * n_gal);
    if (D_e_z2__cm2sm1){for (i = 0; i < n_gal; i++){D_e_z2__cm2sm1[i] = malloc(sizeof *D_e_z2__cm2sm1[i] * n_T_CR);}}





    double G_h = 4.302e-3, eta_pp = 0.5;
    double u_LA, v_Ai, t_loss_s, C[n_gal], Ce_Esm1[n_gal];
    double v_st, v_ste, L_A, Gam_0, D0, tau_eff, B__G[n_gal], B_halo__G[n_gal];


    double CnormE[n_gal];


    #pragma omp parallel for schedule(dynamic) private(j, u_LA, v_Ai, L_A, D0, t_loss_s, v_st, v_ste, Gam_0, seconds, tau_eff )

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



    //    CnormE[i] = C_norm_E(E_cut[i]);
        CnormE[i] = C_norm_E( q_p_inject, m_p__GeV, T_p_cutoff__GeV);

        D0 = v_Ai * L_A * 1e5 * pc__cm;

        t_loss_s = 1./(1./(1./( n__cmm3[i] * sigma_pp_cm2 * eta_pp * c__cmsm1 )) + 1./(pow( h__pc[i] * pc__cm,2)/D0) );

        C[i] = SFR__Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg__GeV/yr__s * t_loss_s/( CnormE[i] * 2. * A_Re__pc2[i] * 2. * h__pc[i] * pow(pc__cm, 3) );

        Ce_Esm1[i] = f_CRe_CRp * SFR__Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg__GeV/( yr__s * C_norm_E( q_e_inject, m_e__GeV, T_e_cutoff__GeV ) );

        for (j = 0; j < n_T_CR; j++)
        {
            v_st = fmin( f_vAi * v_Ai * (1. + 2.3e-3 * pow( sqrt(pow(T_CR__GeV[j],2) + 2. * m_p__GeV * T_CR__GeV[j]) , q_p_inject-1.) * 
                   pow(n_H__cmm3[i]/1e3, 1.5) * (chi/1e-4) * M_A/( u_LA/10. * C[i]/2e-7 )), c__cmsm1/1e5);
            D__cm2sm1[i][j] = v_st * L_A * 1e5 * pc__cm;


            tau_eff = 9.9 * Sig_gas__Msolpcm2[i]/1e3 * h__pc[i]/1e2 * 1e27/D__cm2sm1[i][j];
            Gam_0 = 41.2 * h__pc[i]/1e2 * v_st/1e3 * 1e27/D__cm2sm1[i][j];
            f_cal[i][j] = 1. - 1./( gsl_sf_hyperg_0F1( beta/(beta+1.) , tau_eff/pow(beta+1.,2) ) + 
                       tau_eff/Gam_0 * gsl_sf_hyperg_0F1( (beta+2.)/(beta+1.) , tau_eff/pow(beta+1., 2)) );

if (i == 10){f_cal[i][j] = f_cal[i][j] * 0.1;}

        }

        for (j = 0; j < n_T_CR; j++)
        {
            v_ste = fmin( f_vAi * v_Ai * (1. + 2.3e-3 * pow( sqrt(pow(T_CR__GeV[j],2) + 2. * m_e__GeV * T_CR__GeV[j]) , q_p_inject-1.) * 
                    pow(n_H__cmm3[i]/1e3, 1.5) * (chi/1e-4) * M_A/( u_LA/10. * C[i]/2e-7 )), c__cmsm1/1e5);
            D_e__cm2sm1[i][j] = v_ste * L_A * 1e5 * pc__cm;
            D_e_z2__cm2sm1[i][j] = fmin( v_Ai * (1. + 2.3e-3 * pow( sqrt(pow(T_CR__GeV[j],2) + 2. * m_e__GeV * T_CR__GeV[j]) , q_p_inject-1.) * 
                                pow(n_H__cmm3[i]/1e3/1e3, 1.5) * (1./1e-4) * M_A/( u_LA/10. * ((1.-f_cal[i][0]) * C[i])/2e-7 )),
                                c__cmsm1/1e5) * L_A * 1e5 * pc__cm;
        }

    }

    write_2D_file( n_gal, n_T_CR, f_cal, "", string_cat(outfp, "/fcal.txt") );

    FILE *galdata_out = fopen( string_cat(outfp, "/gal_data.txt"), "w+" );
    fprintf( galdata_out, "h__pc n_H__cmm3 B__G sigmag__kmsm1 Are__pc2 Sigmagas__Msolpcm2 SigmaSFR__Msolyrm1pcm2 Sigmastar__Msolpcm2 Tdust__K\n" );
    for (i = 0; i < n_gal; i++)
    {
        fprintf( galdata_out, "%e %e %e %e %e %e %e %e %e\n", h__pc[i], n_H__cmm3[i], B__G[i], sig_gas__kmsm1[i], A_Re__pc2[i], 
                 Sig_gas__Msolpcm2[i], Sig_SFR__Msolyrm1pcm2[i], Sig_star__Msolpcm2[i], T_dust__K[i] );
    }
    fclose(galdata_out);



/*################################################################################################################################*/
/*################################################################################################################################*/

/* Calculate spectra */


    double **tau_gg = malloc(sizeof *tau_gg * n_gal);
    if (tau_gg){for (i = 0; i < n_gal; i++){tau_gg[i] = malloc(sizeof *tau_gg[i] * n_E_gam);}}
    double **tau_EBL = malloc(sizeof *tau_EBL * n_gal);
    if (tau_EBL){for (i = 0; i < n_gal; i++){tau_EBL[i] = malloc(sizeof *tau_EBL[i] * n_E_gam);}}






//  double **specs_casc_obs = malloc(sizeof *specs_casc_obs * n_gal);
//  if (specs_casc_obs){for (i = 0; i < n_gal; i++){specs_casc_obs[i] = malloc(sizeof *specs_casc_obs[i] * n_E_gam);}}
//  double **specs_casc = malloc(sizeof *specs_casc * n_gal);
//  if (specs_casc){for (i = 0; i < n_gal; i++){specs_casc[i] = malloc(sizeof *specs_casc[i] * n_E_gam);}}



//  double **specs_obs = malloc(sizeof *specs_obs * n_gal);
//  if (specs_obs){for (i = 0; i < n_gal; i++){specs_obs[i] = malloc(sizeof *specs_obs[i] * n_E_gam);}}

//  double **specs_L_emit = malloc(sizeof *specs_L_emit * n_gal);
//  if (specs_L_emit){for (i = 0; i < n_gal; i++){specs_L_emit[i] = malloc(sizeof *specs_L_emit[i] * n_E_gam);}}





//    double C_gam;
//  double spec_emit[n_Esteps];



    double nphot_params[7];





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
    data_object_2D_free( do_2D_BS );


    data_object_1D do_1D_SY = load_SY_do_files( &(n_pts[0]), x_SY_lims, datadir );
    gsl_spline_object_1D gso1D_SY = do1D_to_gso1D( do_1D_SY ); //init_gso_1D_sync( x_SY_lims, 1000, "log" );
    data_object_1D_free( do_1D_SY );



    FILE *fileout = fopen( string_cat(outfp, "/Urad_Ub.txt"), "w+" );
    fprintf( fileout, "u_rad__eVcmm3 u_B__eVcmm3 u_rad_CMB__eVcmm3 u_rad_FIR__eVcmm3 u_rad_3000__eVcmm3 u_rad_4000__eVcmm3 u_rad_7500__eVcmm3 u_rad_UV__eVcmm3\n" );// P_kin__eVcmm3 P_CR__eVcmm3, P_grav__eVcmm3 f_edd\n" );
    for (i = 0; i < n_gal; i++)
    {

        nphot_params[0] = T_0_CMB__K * (1.+z[i]);
        nphot_params[1] = T_dust__K[i];
        nphot_params[2] = M_star__Msol[i];
        nphot_params[3] = SFR__Msolyrm1[i];
        nphot_params[4] = Re__kpc[i];
        nphot_params[5] = h__pc[i];

        fprintf( fileout, "%e %e %e %e %e %e %e %e\n", 
            ISRF_integrate__GeVcmm3( dndEphot_total__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9, //u_rad
            pow(B__G[i],2)/(8.*M_PI)/GeV__erg * 1.e9, //u_B
            ISRF_integrate__GeVcmm3( dndEphot_CMB__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9,
            ISRF_integrate__GeVcmm3( dndEphot_FIR__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9,
            ISRF_integrate__GeVcmm3( dndEphot_3000__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9,
            ISRF_integrate__GeVcmm3( dndEphot_4000__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9,
            ISRF_integrate__GeVcmm3( dndEphot_7500__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9,
            ISRF_integrate__GeVcmm3( dndEphot_UV__cmm3GeVm1, nphot_params, E_phot_lims__GeV ) * 1.e9 );
//            n_H__cmm3[i]*m_H__g*pow(sig_gas__kmsm1[i],2)*1.e10/GeV__erg * 1.e9,
//            1./3. * spec_integrate( n_Esteps, E_GeV, CR_spec[i] )/(2. * A_Re__pc2[i] * 2. * h__pc[i] * pow(pc__cm, 3))  * 1.e9,
//            M_PI * Sig_gas__Msolpcm2[i]/2. * ( Sig_star__Msolpcm2[i]/2. + Sig_gas__Msolpcm2[i]/2. ) * 4.30091e-3 * 1.e10/GeV__erg * 1.e9 * Msol__g/pow(pc__cm,3),
//            (Sig_SFR__Msolyrm1pcm2[i]/4.9e-10)*( Sig_gas__Msolpcm2[i]/(Sig_gas__Msolpcm2[i]+Sig_star__Msolpcm2[i]) ) / ( pow(Sig_gas__Msolpcm2[i]/10.,2) * sig_gas__kmsm1[i]/10. * ( 1./(pow(2.*chi,0.5) * M_A ) ) ) );
    }
    fclose(fileout);

//    return 0;


    write_radspecs( n_gal, z, T_dust__K, M_star__Msol, SFR__Msolyrm1, Re__kpc, h__pc, outfp );










    //loss times SY, BS, IC, IO, DI
    double **tau_loss_z1_SY = malloc(sizeof *tau_loss_z1_SY * n_gal);
    if (tau_loss_z1_SY){for (i = 0; i < n_gal; i++){tau_loss_z1_SY[i] = malloc(sizeof *tau_loss_z1_SY[i] * n_T_CR);}}
    double **tau_loss_z1_BS = malloc(sizeof *tau_loss_z1_BS * n_gal);
    if (tau_loss_z1_BS){for (i = 0; i < n_gal; i++){tau_loss_z1_BS[i] = malloc(sizeof *tau_loss_z1_BS[i] * n_T_CR);}}
    double **tau_loss_z1_IC = malloc(sizeof *tau_loss_z1_IC * n_gal);
    if (tau_loss_z1_IC){for (i = 0; i < n_gal; i++){tau_loss_z1_IC[i] = malloc(sizeof *tau_loss_z1_IC[i] * n_T_CR);}}
    double **tau_loss_z1_IO = malloc(sizeof *tau_loss_z1_IO * n_gal);
    if (tau_loss_z1_IO){for (i = 0; i < n_gal; i++){tau_loss_z1_IO[i] = malloc(sizeof *tau_loss_z1_IO[i] * n_T_CR);}}
    double **tau_loss_z1_DI = malloc(sizeof *tau_loss_z1_DI * n_gal);
    if (tau_loss_z1_DI){for (i = 0; i < n_gal; i++){tau_loss_z1_DI[i] = malloc(sizeof *tau_loss_z1_DI[i] * n_T_CR);}}

    double **tau_loss_z2_SY = malloc(sizeof *tau_loss_z2_SY * n_gal);
    if (tau_loss_z2_SY){for (i = 0; i < n_gal; i++){tau_loss_z2_SY[i] = malloc(sizeof *tau_loss_z2_SY[i] * n_T_CR);}}
    double **tau_loss_z2_BS = malloc(sizeof *tau_loss_z2_BS * n_gal);
    if (tau_loss_z2_BS){for (i = 0; i < n_gal; i++){tau_loss_z2_BS[i] = malloc(sizeof *tau_loss_z2_BS[i] * n_T_CR);}}
    double **tau_loss_z2_IC = malloc(sizeof *tau_loss_z2_IC * n_gal);
    if (tau_loss_z2_IC){for (i = 0; i < n_gal; i++){tau_loss_z2_IC[i] = malloc(sizeof *tau_loss_z2_IC[i] * n_T_CR);}}
    double **tau_loss_z2_IO = malloc(sizeof *tau_loss_z2_IO * n_gal);
    if (tau_loss_z2_IO){for (i = 0; i < n_gal; i++){tau_loss_z2_IO[i] = malloc(sizeof *tau_loss_z2_IO[i] * n_T_CR);}}
    double **tau_loss_z2_DI = malloc(sizeof *tau_loss_z2_DI * n_gal);
    if (tau_loss_z2_DI){for (i = 0; i < n_gal; i++){tau_loss_z2_DI[i] = malloc(sizeof *tau_loss_z2_DI[i] * n_T_CR);}}

    double **tau_loss_protons_PP = malloc(sizeof *tau_loss_protons_PP * n_gal);
    if (tau_loss_protons_PP){for (i = 0; i < n_gal; i++){tau_loss_protons_PP[i] = malloc(sizeof *tau_loss_protons_PP[i] * n_T_CR);}}
    double **tau_loss_protons_DI = malloc(sizeof *tau_loss_protons_DI * n_gal);
    if (tau_loss_protons_DI){for (i = 0; i < n_gal; i++){tau_loss_protons_DI[i] = malloc(sizeof *tau_loss_protons_DI[i] * n_T_CR);}}



    double **spec_pi = malloc(sizeof *spec_pi * n_gal);
    if (spec_pi){for (i = 0; i < n_gal; i++){spec_pi[i] = malloc(sizeof *spec_pi[i] * n_E_gam);}}

    double **spec_pi_fcal1 = malloc(sizeof *spec_pi_fcal1 * n_gal);
    if (spec_pi_fcal1){for (i = 0; i < n_gal; i++){spec_pi_fcal1[i] = malloc(sizeof *spec_pi_fcal1[i] * n_E_gam);}}

    double **spec_nu = malloc(sizeof *spec_nu * n_gal);
    if (spec_nu){for (i = 0; i < n_gal; i++){spec_nu[i] = malloc(sizeof *spec_nu[i] * n_E_gam);}}



    double **spec_IC_1_z1 = malloc(sizeof *spec_IC_1_z1 * n_gal);
    if (spec_IC_1_z1){for (i = 0; i < n_gal; i++){spec_IC_1_z1[i] = malloc(sizeof *spec_IC_1_z1[i] * n_E_gam);}}
    double **spec_IC_2_z1 = malloc(sizeof *spec_IC_2_z1 * n_gal);
    if (spec_IC_2_z1){for (i = 0; i < n_gal; i++){spec_IC_2_z1[i] = malloc(sizeof *spec_IC_2_z1[i] * n_E_gam);}}
    double **spec_IC_1_z2 = malloc(sizeof *spec_IC_1_z2 * n_gal);
    if (spec_IC_1_z2){for (i = 0; i < n_gal; i++){spec_IC_1_z2[i] = malloc(sizeof *spec_IC_1_z2[i] * n_E_gam);}}
    double **spec_IC_2_z2 = malloc(sizeof *spec_IC_2_z2 * n_gal);
    if (spec_IC_2_z2){for (i = 0; i < n_gal; i++){spec_IC_2_z2[i] = malloc(sizeof *spec_IC_2_z2[i] * n_E_gam);}}

    double **spec_BS_1_z1 = malloc(sizeof *spec_BS_1_z1 * n_gal);
    if (spec_BS_1_z1){for (i = 0; i < n_gal; i++){spec_BS_1_z1[i] = malloc(sizeof *spec_BS_1_z1[i] * n_E_gam);}}
    double **spec_BS_2_z1 = malloc(sizeof *spec_BS_2_z1 * n_gal);
    if (spec_BS_2_z1){for (i = 0; i < n_gal; i++){spec_BS_2_z1[i] = malloc(sizeof *spec_BS_2_z1[i] * n_E_gam);}}
//    double **spec_BS_1_z2 = malloc(sizeof *spec_BS_1_z2 * n_gal);
//    if (spec_BS_1_z2){for (i = 0; i < n_gal; i++){spec_BS_1_z2[i] = malloc(sizeof *spec_BS_1_z2[i] * n_E_gam);}}
//    double **spec_BS_2_z2 = malloc(sizeof *spec_BS_2_z2 * n_gal);
//    if (spec_BS_2_z2){for (i = 0; i < n_gal; i++){spec_BS_2_z2[i] = malloc(sizeof *spec_BS_2_z2[i] * n_E_gam);}}

    double **spec_SY_1_z1 = malloc(sizeof *spec_SY_1_z1 * n_gal);
    if (spec_SY_1_z1){for (i = 0; i < n_gal; i++){spec_SY_1_z1[i] = malloc(sizeof *spec_SY_1_z1[i] * n_E_gam);}}
    double **spec_SY_2_z1 = malloc(sizeof *spec_SY_2_z1 * n_gal);
    if (spec_SY_2_z1){for (i = 0; i < n_gal; i++){spec_SY_2_z1[i] = malloc(sizeof *spec_SY_2_z1[i] * n_E_gam);}}
    double **spec_SY_1_z2 = malloc(sizeof *spec_SY_1_z2 * n_gal);
    if (spec_SY_1_z2){for (i = 0; i < n_gal; i++){spec_SY_1_z2[i] = malloc(sizeof *spec_SY_1_z2[i] * n_E_gam);}}
    double **spec_SY_2_z2 = malloc(sizeof *spec_SY_2_z2 * n_gal);
    if (spec_SY_2_z2){for (i = 0; i < n_gal; i++){spec_SY_2_z2[i] = malloc(sizeof *spec_SY_2_z2[i] * n_E_gam);}}




    double **Q_e_1_z1 = malloc(sizeof *Q_e_1_z1 * n_gal);
    if (Q_e_1_z1){for (i = 0; i < n_gal; i++){Q_e_1_z1[i] = malloc(sizeof *Q_e_1_z1[i] * n_T_CR);}}
    double **Q_e_1_z2 = malloc(sizeof *Q_e_1_z2 * n_gal);
    if (Q_e_1_z2){for (i = 0; i < n_gal; i++){Q_e_1_z2[i] = malloc(sizeof *Q_e_1_z2[i] * n_T_CR);}}
    double **Q_e_2_z1 = malloc(sizeof *Q_e_2_z1 * n_gal);
    if (Q_e_2_z1){for (i = 0; i < n_gal; i++){Q_e_2_z1[i] = malloc(sizeof *Q_e_2_z1[i] * n_T_CR);}}
    double **Q_e_2_z2 = malloc(sizeof *Q_e_2_z2 * n_gal);
    if (Q_e_2_z2){for (i = 0; i < n_gal; i++){Q_e_2_z2[i] = malloc(sizeof *Q_e_2_z2[i] * n_T_CR);}}

    double **q_p_SS_z1 = malloc(sizeof *q_p_SS_z1 * n_gal);
    if (q_p_SS_z1){for (i = 0; i < n_gal; i++){q_p_SS_z1[i] = malloc(sizeof *q_p_SS_z1[i] * n_T_CR);}}
    double **q_e_SS_1_z1 = malloc(sizeof *q_e_SS_1_z1 * n_gal);
    if (q_e_SS_1_z1){for (i = 0; i < n_gal; i++){q_e_SS_1_z1[i] = malloc(sizeof *q_e_SS_1_z1[i] * n_T_CR);}}
    double **q_e_SS_2_z1 = malloc(sizeof *q_e_SS_2_z1 * n_gal);
    if (q_e_SS_2_z1){for (i = 0; i < n_gal; i++){q_e_SS_2_z1[i] = malloc(sizeof *q_e_SS_2_z1[i] * n_T_CR);}}
    double **q_e_SS_1_z2 = malloc(sizeof *q_e_SS_1_z2 * n_gal);
    if (q_e_SS_1_z2){for (i = 0; i < n_gal; i++){q_e_SS_1_z2[i] = malloc(sizeof *q_e_SS_1_z2[i] * n_T_CR);}}
    double **q_e_SS_2_z2 = malloc(sizeof *q_e_SS_2_z2 * n_gal);
    if (q_e_SS_2_z2){for (i = 0; i < n_gal; i++){q_e_SS_2_z2[i] = malloc(sizeof *q_e_SS_2_z2[i] * n_T_CR);}}


    double **tau_FF = malloc(sizeof *tau_FF * n_gal);
    if (tau_FF){for (i = 0; i < n_gal; i++){tau_FF[i] = malloc(sizeof *tau_FF[i] * n_E_gam);}}


    double **spec_FF = malloc(sizeof *spec_FF * n_gal);
    if (spec_FF){for (i = 0; i < n_gal; i++){spec_FF[i] = malloc(sizeof *spec_FF[i] * n_E_gam);}}



    double E_crit__GeV;

    double **E_loss_leptons = malloc(sizeof *E_loss_leptons * n_gal);
    if (E_loss_leptons){for (i = 0; i < n_gal; i++){E_loss_leptons[i] = malloc(sizeof *E_loss_leptons[i] * 16);}}

    double **E_loss_nucrit = malloc(sizeof *E_loss_nucrit * n_gal);
    if (E_loss_nucrit){for (i = 0; i < n_gal; i++){E_loss_nucrit[i] = malloc(sizeof *E_loss_nucrit[i] * 10);}}

    double **Lradio = malloc(sizeof *Lradio * n_gal);
    if (Lradio){for (i = 0; i < n_gal; i++){Lradio[i] = malloc(sizeof *Lradio[i] * 5);}}

    gsl_spline_object_1D gso1D_fcal;



    gsl_spline_object_1D De_gso1D_z1, De_gso1D_z2, gso_1D_Q_inject_1_z1, gso_1D_Q_inject_2_z1, gso_1D_Q_inject_1_z2, gso_1D_Q_inject_2_z2, qe_1_z1_so, qe_2_z1_so, qe_1_z2_so, qe_2_z2_so;

    gsl_spline_object_2D gso2D_IC;
    gsl_spline_object_2D gso2D_IC_Gamma;

    printf("%s \n", "Calculating spectra:");
    fflush(stdout);

    #pragma omp parallel for schedule(guided) private(j, k, qe_1_z1_so, qe_2_z1_so, qe_1_z2_so, qe_2_z2_so, De_gso1D_z1, De_gso1D_z2, nphot_params, gso_1D_Q_inject_1_z1, gso_1D_Q_inject_2_z1, gso_1D_Q_inject_1_z2, gso_1D_Q_inject_2_z2, E_crit__GeV, seconds, gso2D_IC, gso2D_IC_Gamma, gso1D_fcal, t_loss_s ) firstprivate( gso2D_BS, gso1D_SY )

    for (i = 0; i < n_gal; i++)
    {
        seconds = time(NULL);
        printf("Working on galaxy no. %lu\n", i);
        fflush(stdout);

        nphot_params[0] = T_0_CMB__K * (1.+z[i]);
        nphot_params[1] = T_dust__K[i];
        nphot_params[2] = M_star__Msol[i];
        nphot_params[3] = SFR__Msolyrm1[i];
        nphot_params[4] = Re__kpc[i];
        nphot_params[5] = h__pc[i];
        nphot_params[6] = R_vir__kpc( R_half_mass__kpc( Re__kpc[i], z[i] ) );
    
        gso2D_IC = construct_IC_gso2D( 0, ICo, nphot_params );
        gso2D_IC_Gamma = construct_IC_gso2D( 1, ICo, nphot_params );


        //C here is used for the calcultion with dsig/dE, so need t_loss = t_col
        t_loss_s = 1./( n__cmm3[i] * sigma_pp_cm2 * eta_pp * c__cmsm1 );
        C[i] = SFR__Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg__GeV/yr__s * t_loss_s/CnormE[i];


        gso1D_fcal = gsl_so1D( n_T_CR, T_CR__GeV, f_cal[i] );



        //Zone 1 diffusion
        De_gso1D_z1 = gsl_so1D( n_T_CR, E_CRe__GeV, D_e__cm2sm1[i] );


        //primary injection spectrum and steady state spline object per galaxy


        for (j = 0; j < n_T_CR; j++)
        {
            Q_e_1_z1[i][j] = J( T_CR__GeV[j], Ce_Esm1[i], q_e_inject, m_e__GeV, T_e_cutoff__GeV );
if (i == 10){Q_e_1_z1[i][j] = J( T_CR__GeV[j], Ce_Esm1[i], q_e_inject, m_e__GeV, T_e_cutoff__GeV ) * 0.1;}
        }
        //spline object on total energy - not kinetic
        gso_1D_Q_inject_1_z1 = gsl_so1D( n_T_CR, E_CRe__GeV, Q_e_1_z1[i] );

        //secondary injection spectrum and steady state spline object
        for (j = 0; j < n_T_CR; j++)
        {
            //Compute the spectra for the secondary electrons and then interpolate on them
            Q_e_2_z1[i][j] = q_e( T_CR__GeV[j], n_H__cmm3[i], C[i], T_p_cutoff__GeV, gso1D_fcal );
        }
        //spline object on total energy - not kinetic
        gso_1D_Q_inject_2_z1 = gsl_so1D( n_T_CR, E_CRe__GeV, Q_e_2_z1[i] );

//        printf("Before SS %lu\n", i);
//        fflush(stdout);

        CRe_steadystate_solve( 1, E_CRe_lims__GeV, 500, n_H__cmm3[i], B__G[i], h__pc[i], 1, &gso2D_IC_Gamma, gso2D_BS, De_gso1D_z1, 
                               gso_1D_Q_inject_1_z1, gso_1D_Q_inject_2_z1, &qe_1_z1_so, &qe_2_z1_so );
            
 
        //Zone 2 diffusion
        De_gso1D_z2 = gsl_so1D( n_T_CR, E_CRe__GeV, D_e_z2__cm2sm1[i] );

        //primary injection spectrum and steady state spline object per galaxy
        for (j = 0; j < n_T_CR; j++)
        {
            Q_e_1_z2[i][j] = gsl_so1D_eval( qe_1_z1_so, E_CRe__GeV[j] )/tau_diff__s( E_CRe__GeV[j], h__pc[i], De_gso1D_z1 );
        }
        gso_1D_Q_inject_1_z2 = gsl_so1D( n_T_CR, E_CRe__GeV, Q_e_1_z2[i] );


        //secondary injection spectrum and steady state spline object
        for (j = 0; j < n_T_CR; j++)
        {
            Q_e_2_z2[i][j] = gsl_so1D_eval( qe_2_z1_so, E_CRe__GeV[j] )/tau_diff__s( E_CRe__GeV[j], h__pc[i], De_gso1D_z1 );
        }

        gso_1D_Q_inject_2_z2 = gsl_so1D( n_T_CR, E_CRe__GeV, Q_e_2_z2[i] );

        CRe_steadystate_solve( 2, E_CRe_lims__GeV, 500, n_H__cmm3[i]/1000., B_halo__G[i], 50.*h__pc[i], 1, &gso2D_IC_Gamma, gso2D_BS,
                               De_gso1D_z2, gso_1D_Q_inject_1_z2, gso_1D_Q_inject_2_z2, &qe_1_z2_so, &qe_2_z2_so );


//Fix cascade
//        C_gam = norm_casc_C( spec_emit, E_GeV, n_Esteps, z[i], fdata_in );

        for (j = 0; j < n_E_gam; j++)
        {

            tau_gg[i][j] = tau_gg_gal_BW( E_gam__GeV[j], dndEphot_total__cmm3GeVm1, nphot_params, E_phot_lims__GeV, h__pc[i] );

            tau_FF[i][j] = tau_FF_MK( E_gam__GeV[j], Sig_SFR__Msolyrm1pcm2[i], 1.e4 );

            spec_pi[i][j] = eps_pi( E_gam__GeV[j], n_H__cmm3[i], C[i], T_p_cutoff__GeV, gso1D_fcal );


            spec_IC_1_z1[i][j] = eps_IC_3( E_gam__GeV[j], gso2D_IC, qe_1_z1_so ) * exp(-tau_FF[i][j]);
            spec_IC_2_z1[i][j] = eps_IC_3( E_gam__GeV[j], gso2D_IC, qe_2_z1_so ) * exp(-tau_FF[i][j]);

            spec_BS_1_z1[i][j] = eps_BS_3( E_gam__GeV[j], n_H__cmm3[i], gso2D_BS, qe_1_z1_so ) * exp(-tau_FF[i][j]);
            spec_BS_2_z1[i][j] = eps_BS_3( E_gam__GeV[j], n_H__cmm3[i], gso2D_BS, qe_2_z1_so ) * exp(-tau_FF[i][j]);

            spec_SY_1_z1[i][j] = eps_SY_4( E_gam__GeV[j], B__G[i], gso1D_SY, qe_1_z1_so ) * exp(-tau_FF[i][j]);
            spec_SY_2_z1[i][j] = eps_SY_4( E_gam__GeV[j], B__G[i], gso1D_SY, qe_2_z1_so ) * exp(-tau_FF[i][j]);


            spec_IC_1_z2[i][j] = eps_IC_3( E_gam__GeV[j], gso2D_IC, qe_1_z2_so );
            spec_IC_2_z2[i][j] = eps_IC_3( E_gam__GeV[j], gso2D_IC, qe_2_z2_so );

            spec_SY_1_z2[i][j] = eps_SY_4( E_gam__GeV[j], B_halo__G[i], gso1D_SY, qe_1_z2_so );
            spec_SY_2_z2[i][j] = eps_SY_4( E_gam__GeV[j], B_halo__G[i], gso1D_SY, qe_2_z2_so );




            spec_pi_fcal1[i][j] = eps_pi_fcal1( E_gam__GeV[j], n_H__cmm3[i], C[i], T_p_cutoff__GeV, gso1D_fcal );

            spec_nu[i][j] = q_nu( E_gam__GeV[j], n_H__cmm3[i], C[i], T_p_cutoff__GeV, gso1D_fcal );




            spec_FF[i][j] = eps_FF( E_gam__GeV[j], Re__kpc[i], 1.e4, tau_FF[i][j] );


//            specs_obs[i][j] = (Phi_out.Phi + specs_brems_1st[i][j] + specs_brems_2nd[i][j] + specs_IC_1st[i][j] + specs_IC_2nd[i][j] + 
//                               specs_sync_1st[i][j] + specs_sync_2nd[i][j])* mod * exp(-tau_gg[i][j]) * exp(-tau_EBL[i][j]);

            //call the function again this time without redshift
/*
            Phi_out = Phi( E_GeV[j], n_H__cmm3[i], C[i], T_p_cutoff__GeV, gso1D_fcal );


            sIC1 = eps_IC_3( E_GeV[j], gso_2D_total, qe_1_z1_so );
            sIC2 = eps_IC_3( E_GeV[j], gso_2D_total, qe_2_z1_so );
            sb1 = eps_BS_3( E_GeV[j], n_H__cmm3[i], qe_1_z1_so );
            sb2 = eps_BS_3( E_GeV[j], n_H__cmm3[i], qe_2_z1_so );

            ss1 = eps_SY_4( E_GeV[j], B__G[i], gso1D_SY, qe_1_z1_so );
            ss2 = eps_SY_4( E_GeV[j], B__G[i], gso1D_SY, qe_2_z1_so );


            //total energy emitted
            specs_L_emit[i][j] = ( Phi_out.Phi * vol + sb1 + sb2 +sIC1 + sIC2 + ss1 + ss2 ) * exp(-tau_gg[i][j]);
            //for cascade
            spec_emit[j] = ( Phi_out.Phi + sb1 + sb2 +sIC1 + sIC2 + ss1 + ss2 )/vol * exp(-tau_gg[i][j]);

            specs_casc_obs[i][j] = dndE_gam_casc( (1.+z[i]) * E_GeV[j], z[i], C_gam, fdata_in ) * mod;
            specs_casc[i][j] = dndE_gam_casc( E_GeV[j], z[i], C_gam, fdata_in ) * vol;
*/


        }

/*
    if (i == 0)
    {
        gsl_spline_object_2D gso2D_3000, gso2D_4000, gso2D_7500, gso2D_UV, gso2D_FIR, gso2D_CMB;

        return_IC_gso2D( ICo, nphot_params, &gso2D_3000, &gso2D_4000, &gso2D_7500, &gso2D_UV, &gso2D_FIR, &gso2D_CMB );

        write_ICspectra( n_E_gam, E_gam__GeV, qe_1_z1_so, qe_2_z1_so, qe_1_z2_so, qe_2_z2_so, gso2D_3000, gso2D_4000, gso2D_7500, 
                         gso2D_UV, gso2D_FIR, gso2D_CMB, string_cat(outfp, "/IC_spectra.txt") );

        gsl_so2D_free( gso2D_3000 );
        gsl_so2D_free( gso2D_4000 );
        gsl_so2D_free( gso2D_7500 );
        gsl_so2D_free( gso2D_UV );
        gsl_so2D_free( gso2D_FIR );
        gsl_so2D_free( gso2D_CMB );        
    }
*/





        E_crit__GeV = sqrt( (2. * 1.49e9 * m_e__g*c__cmsm1)/(3. * B__G[i]*e__esu) ) * M_PI * m_e__GeV;
        E_loss_nucrit[i][0] = E_crit__GeV/tau_BS_fulltest__s( E_crit__GeV, E_CRe_lims__GeV, n_H__cmm3[i], gso2D_BS );
        E_loss_nucrit[i][1] = E_crit__GeV/tau_sync__s( E_crit__GeV, B__G[i] );
        E_loss_nucrit[i][2] = E_crit__GeV/tau_IC_fulltest__s( E_crit__GeV, E_CRe_lims__GeV, gso2D_IC_Gamma );
        E_loss_nucrit[i][3] = E_crit__GeV/tau_ion__s( E_crit__GeV, n_H__cmm3[i] );
        E_loss_nucrit[i][4] = E_crit__GeV/tau_diff__s( E_crit__GeV, h__pc[i], De_gso1D_z1 );
        E_crit__GeV = sqrt( (2. * 1.49e9 * m_e__g*c__cmsm1)/(3. * B_halo__G[i]*e__esu) ) * M_PI * m_e__GeV;
        E_loss_nucrit[i][5] = E_crit__GeV/tau_BS_fulltest__s( E_crit__GeV, E_CRe_lims__GeV, n_H__cmm3[i]/1000., gso2D_BS );
        E_loss_nucrit[i][6] = E_crit__GeV/tau_sync__s( E_crit__GeV, B_halo__G[i] );
        E_loss_nucrit[i][7] = E_crit__GeV/tau_IC_fulltest__s( E_crit__GeV, E_CRe_lims__GeV, gso2D_IC_Gamma );
        E_loss_nucrit[i][8] = E_crit__GeV/tau_plasma__s( E_crit__GeV, n_H__cmm3[i]/1000. );
        E_loss_nucrit[i][9] = E_crit__GeV/tau_diff__s( E_crit__GeV, 50.*h__pc[i], De_gso1D_z2 );


        Lradio[i][0] = eps_SY_4( 1.49e9 * h__GeVs, B__G[i], gso1D_SY, qe_1_z1_so ) * 1.49e9 * h__Js * h__GeVs * exp(-tau_FF_MK( 1.49e9 * h__GeVs, Sig_SFR__Msolyrm1pcm2[i], 1.e4 ));
        Lradio[i][1] = eps_SY_4( 1.49e9 * h__GeVs, B__G[i], gso1D_SY, qe_2_z1_so ) * 1.49e9 * h__Js * h__GeVs * exp(-tau_FF_MK( 1.49e9 * h__GeVs, Sig_SFR__Msolyrm1pcm2[i], 1.e4 ));
        Lradio[i][2] = eps_SY_4( 1.49e9 * h__GeVs, B_halo__G[i], gso1D_SY, qe_1_z2_so ) * 1.49e9 * h__Js * h__GeVs;
        Lradio[i][3] = eps_SY_4( 1.49e9 * h__GeVs, B_halo__G[i], gso1D_SY, qe_2_z2_so ) * 1.49e9 * h__Js * h__GeVs;
        Lradio[i][4] = eps_FF( 1.49e9 * h__GeVs, Re__kpc[i], 1.e4, tau_FF_MK( 1.49e9 * h__GeVs, Sig_SFR__Msolyrm1pcm2[i], 1.e4 ) ) * 1.49e9 * h__Js * h__GeVs;


        for (j = 0; j < n_T_CR; j++)
        {
            q_p_SS_z1[i][j] = J( T_CR__GeV[j], C[i], q_p_inject, m_p__GeV, T_p_cutoff__GeV ) * f_cal[i][j];
            q_e_SS_1_z1[i][j] = gsl_so1D_eval( qe_1_z1_so, E_CRe__GeV[j] );
            q_e_SS_2_z1[i][j] = gsl_so1D_eval( qe_2_z1_so, E_CRe__GeV[j] );
            q_e_SS_1_z2[i][j] = gsl_so1D_eval( qe_1_z2_so, E_CRe__GeV[j] );
            q_e_SS_2_z2[i][j] = gsl_so1D_eval( qe_2_z2_so, E_CRe__GeV[j] );
        }


        E_loss_leptons[i][0] = spec_integrate( n_T_CR, T_CR__GeV, Q_e_1_z1[i] ); //f_CRe_CRp * SFR__Msolyrm1[i] * n_SN_Msolm1 * f_EtoCR * E_SN_erg * erg__GeV/yr__s;
        E_loss_leptons[i][1] = spec_integrate( n_T_CR, T_CR__GeV, Q_e_2_z1[i] );
        E_loss_leptons[i][2] = spec_integrate( n_E_gam, E_gam__GeV, spec_SY_1_z1[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][3] = spec_integrate( n_E_gam, E_gam__GeV, spec_IC_1_z1[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][4] = spec_integrate( n_E_gam, E_gam__GeV, spec_BS_1_z1[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][5] = spec_integrate( n_E_gam, E_gam__GeV, spec_SY_1_z2[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][6] = spec_integrate( n_E_gam, E_gam__GeV, spec_IC_1_z2[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][7] = 0.; //spec_integrate( n_E_gam, E_gam__GeV, spec_BS_1_z2[i] )/E_loss_leptons[i][0];
        E_loss_leptons[i][8] = spec_integrate( n_T_CR, T_CR__GeV, Q_e_1_z2[i] )/E_loss_leptons[i][0];

        E_loss_leptons[i][9] = spec_integrate( n_E_gam, E_gam__GeV, spec_SY_2_z1[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][10] = spec_integrate( n_E_gam, E_gam__GeV, spec_IC_2_z1[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][11] = spec_integrate( n_E_gam, E_gam__GeV, spec_BS_2_z1[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][12] = spec_integrate( n_E_gam, E_gam__GeV, spec_SY_2_z2[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][13] = spec_integrate( n_E_gam, E_gam__GeV, spec_IC_2_z2[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][14] = 0.; //spec_integrate( n_E_gam, E_gam__GeV, spec_BS_2_z2[i] )/E_loss_leptons[i][1];
        E_loss_leptons[i][15] = spec_integrate( n_T_CR, T_CR__GeV, Q_e_2_z2[i] )/E_loss_leptons[i][1];



        //Compute some loss times SY, BS, IC, IO, DI
        for (j = 0; j < n_T_CR; j++)
        {
            tau_loss_z1_SY[i][j] = tau_sync__s( E_CRe__GeV[j], B__G[i] );
            tau_loss_z1_BS[i][j] = tau_BS_fulltest__s( E_CRe__GeV[j], E_CRe_lims__GeV, n_H__cmm3[i], gso2D_BS );
            tau_loss_z1_IC[i][j] = tau_IC_fulltest__s( E_CRe__GeV[j], E_CRe_lims__GeV, gso2D_IC_Gamma );
            tau_loss_z1_IO[i][j] = tau_ion__s( E_CRe__GeV[j], n_H__cmm3[i] );
            tau_loss_z1_DI[i][j] = tau_diff__s( E_CRe__GeV[j], h__pc[i], De_gso1D_z1 );

            tau_loss_z2_SY[i][j] = tau_sync__s( E_CRe__GeV[j], B_halo__G[i] );
            tau_loss_z2_BS[i][j] = tau_BS_fulltest__s( E_CRe__GeV[j], E_CRe_lims__GeV, n_H__cmm3[i]/1000., gso2D_BS );
            tau_loss_z2_IC[i][j] = tau_IC_fulltest__s( E_CRe__GeV[j], E_CRe_lims__GeV, gso2D_IC_Gamma );
            tau_loss_z2_IO[i][j] = tau_ion__s( E_CRe__GeV[j], n_H__cmm3[i]/1000. );
            tau_loss_z2_DI[i][j] = tau_diff__s( E_CRe__GeV[j], h__pc[i]*50., De_gso1D_z2 );

            tau_loss_protons_PP[i][j] = 1./( n__cmm3[i] * sigma_pp_cm2 * eta_pp * c__cmsm1 );
            tau_loss_protons_DI[i][j] = pow(h__pc[i] * pc__cm,2)/D__cm2sm1[i][j];

        }


        gsl_so2D_free( gso2D_IC );
        gsl_so2D_free( gso2D_IC_Gamma );


        gsl_so1D_free( gso1D_fcal );
        gsl_so1D_free( qe_1_z1_so );
        gsl_so1D_free( qe_2_z1_so );
        gsl_so1D_free( qe_1_z2_so );
        gsl_so1D_free( qe_2_z2_so );

        gsl_so1D_free( gso_1D_Q_inject_1_z1 );
        gsl_so1D_free( gso_1D_Q_inject_2_z1 );
        gsl_so1D_free( gso_1D_Q_inject_1_z2 );
        gsl_so1D_free( gso_1D_Q_inject_2_z2 );
        gsl_so1D_free( De_gso1D_z1 );
        gsl_so1D_free( De_gso1D_z2 );

        seconds = time(NULL) - seconds;
        printf("Galaxy %lu log10(time/sec): %le \n", i, log10(seconds));
        fflush(stdout);



//test_IC_interp_obj( E_CRe_lims__GeV, gso_2D_total_low, gso_2D_total, gso_2D_IC_Gamma, qe_1_z1_so );

    }


    double **array2Dlist[4] = { f_cal, D__cm2sm1, D_e__cm2sm1, D_e_z2__cm2sm1 };
    for (i = 0; i < 4; i++)
    {
        free2D( n_gal, array2Dlist[i] );
    }

    gsl_so2D_free( gso2D_BS );
    gsl_so1D_free( gso1D_SY );
    IC_object_free( ICo );


/*################################################################################################################################*/
    printf("%s \n", "Writing output files 1...");
/*################################################################################################################################*/


    write_2D_file( n_gal, n_E_gam, tau_FF, "free-free optical depth at wavelength of emission", string_cat(outfp, "/tau_ff.txt") );


    char taufile[14][34] = { "/tau_loss/tau_loss_z1_SY.txt", "/tau_loss/tau_loss_z2_SY.txt", 
                              "/tau_loss/tau_loss_z1_BS.txt", "/tau_loss/tau_loss_z2_BS.txt",
                              "/tau_loss/tau_loss_z1_IC.txt", "/tau_loss/tau_loss_z2_IC.txt",
                              "/tau_loss/tau_loss_z1_DI.txt", "/tau_loss/tau_loss_z2_DI.txt",
                              "/tau_loss/tau_loss_z1_IO.txt", "/tau_loss/tau_loss_z2_IO.txt",
                              "/tau_loss/tau_loss_protons_PP.txt", "/tau_loss/tau_loss_protons_DI.txt" };

    double **taulist[12] = { tau_loss_z1_SY, tau_loss_z2_SY, tau_loss_z1_BS, tau_loss_z2_BS, 
                              tau_loss_z1_IC, tau_loss_z2_IC, tau_loss_z1_DI, tau_loss_z2_DI, 
                              tau_loss_z1_IO, tau_loss_z2_IO, tau_loss_protons_PP, tau_loss_protons_DI };

    for (k = 0; k < 12; k++)
    {
        write_2D_file( n_gal, n_T_CR, taulist[k], "s", string_cat(outfp, taufile[k]) );
    }

    for (k = 0; k < 12; k++)
    {
        free2D( n_gal, taulist[k] );
    }
    

    write_2D_file( n_gal, 16, E_loss_leptons, "E_inj_1 E_inj_2 fSY1z1 fIC1z1 fBS1z1 fSY1z2 fIC1z2 fBS1z2 fE_inj_1_z2 fSY2z1 fIC2z1 fBS2z1 fSY2z2 fIC2z2 fBS2z2 fE_inj_2_z2", string_cat(outfp, "/E_loss_leptons.txt") );
    free2D( n_gal, E_loss_leptons );
    

    write_2D_file( n_gal, 10, E_loss_nucrit, "BS_z1 SY_z1 IC_z1 IO_z1 DI_z1 BS_z2 SY_z2 IC_z2 IO_z2 DI_z2", string_cat(outfp, "/E_loss_nucrit.txt") );
    free2D( n_gal, E_loss_nucrit );

    FILE *CR_specs_out = fopen( string_cat(outfp, "/CR_specs.txt"), "w+" );
    for (i = 0; i < n_gal; i++)
    {
        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_out, "%e ", q_p_SS_z1[i][j] );
        }
        fprintf( CR_specs_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_out, "%e ", q_e_SS_1_z1[i][j] );
        }
        fprintf( CR_specs_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_out, "%e ", q_e_SS_2_z1[i][j] );
        }
        fprintf( CR_specs_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_out, "%e ", q_e_SS_1_z2[i][j] );
        }
        fprintf( CR_specs_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_out, "%e ", q_e_SS_2_z2[i][j] );
        }
        fprintf( CR_specs_out, "\n" );
    }
    fclose(CR_specs_out);

    free2D( n_gal, q_p_SS_z1 );
    free2D( n_gal, q_e_SS_1_z1 );
    free2D( n_gal, q_e_SS_2_z1 );
    free2D( n_gal, q_e_SS_1_z2 );
    free2D( n_gal, q_e_SS_2_z2 );


    FILE *CR_specs_inj_out = fopen( string_cat(outfp, "/CR_specs_inj.txt"), "w+" );
    for (i = 0; i < n_gal; i++)
    {
//        for (j = 0; j < n_T_CR; j++)
//        {
//            fprintf( CR_specs_out, "%e ", q_p_SS_z1[i][j] );
//        }
//        fprintf( CR_specs_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_inj_out, "%e ", Q_e_1_z1[i][j] );
        }
        fprintf( CR_specs_inj_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_inj_out, "%e ", Q_e_2_z1[i][j] );
        }
        fprintf( CR_specs_inj_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_inj_out, "%e ", Q_e_1_z2[i][j] );
        }
        fprintf( CR_specs_inj_out, "\n" );

        for (j = 0; j < n_T_CR; j++)
        {
            fprintf( CR_specs_inj_out, "%e ", Q_e_2_z2[i][j] );
        }
        fprintf( CR_specs_inj_out, "\n" );
    }
    fclose(CR_specs_inj_out);

//    free2D( n_gal, q_p_SS_z1 );
    free2D( n_gal, Q_e_1_z1 );
    free2D( n_gal, Q_e_2_z1 );
    free2D( n_gal, Q_e_1_z2 );
    free2D( n_gal, Q_e_2_z2 );


    write_2D_file( n_gal, 5, Lradio, "L_radio_1_z1__WHzm1 L_radio_2_z1__WHzm1 L_radio_1_z2__WHzm1 L_radio_2_z2__WHzm1 L_FF__WHzm1", string_cat(outfp, "/L_radio.txt") );
    free2D( n_gal, Lradio );


    double **data = malloc(sizeof *data * n_gal);
    if (data){for (i = 0; i < n_gal; i++){data[i] = malloc(sizeof *data[i] * n_E_gam);}}

    double distmod[n_gal];

    double **E_gam_z0__GeV = malloc(sizeof *E_gam_z0__GeV * n_gal);
    if (E_gam_z0__GeV){for (i = 0; i < n_gal; i++){E_gam_z0__GeV[i] = malloc(sizeof *E_gam_z0__GeV[i] * n_E_gam);}}
    


    for (i = 0; i < n_gal; i++)
    {
        distmod[i] = pow( (1.+z[i]), 2 )/( 4. * M_PI * pow( d_l_MPc( z[i] ) * Mpc__cm, 2 ) );
        for (j = 0; j < n_E_gam; j++)
        {
            E_gam_z0__GeV[i][j] = E_gam__GeV[j]/(1.+z[i]);
            tau_EBL[i][j] = fmax(0., gsl_interp2d_eval_extrap( fdata_in.interp, fdata_in.xa, fdata_in.ya, fdata_in.za, 
                            fmin(E_gam_z0__GeV[i][j] * 1e9, 1e15), z[i], fdata_in.xacc, fdata_in.yacc ));
        }
    }

    write_1D_file( n_gal, distmod, "dist_modulus__cmm2", string_cat(outfp, "/distmod.txt") );
    write_2D_file( n_gal, n_E_gam, tau_FF, "free-free optical depth at wavelength of emission", string_cat(outfp, "/tau_ff.txt") );

    gsl_spline_object_1D gso1D_data;

    char filelist[14][19] = { "/spec_pi.txt", "/spec_IC_1_z1.txt", "/spec_IC_2_z1.txt", "/spec_BS_1_z1.txt", "/spec_BS_2_z1.txt", 
                        "/spec_SY_1_z1.txt", "/spec_SY_2_z1.txt", "/spec_IC_1_z2.txt", "/spec_IC_2_z2.txt", "/spec_SY_1_z2.txt", 
                        "/spec_SY_2_z2.txt", "/spec_pi_fcal1.txt", "/spec_nu.txt", "/spec_FF.txt" };

    double **speclist[14] = { spec_pi, spec_IC_1_z1, spec_IC_2_z1, spec_BS_1_z1, spec_BS_2_z1, spec_SY_1_z1, spec_SY_2_z1,
                           spec_IC_1_z2, spec_IC_2_z2, spec_SY_1_z2, spec_SY_2_z2, spec_pi_fcal1, spec_nu, spec_FF };

    for (k = 0; k < 14; k++)
    {
        for (i = 0; i < n_gal; i++)
        {
            gso1D_data = gsl_so1D( n_E_gam, E_gam_z0__GeV[i], speclist[k][i] );
            for (j = 0; j < n_E_gam; j++)
            {
                data[i][j] = gsl_so1D_eval( gso1D_data, E_gam__GeV[j] );
            }
            gsl_so1D_free( gso1D_data );
        }
        write_2D_spec_file( n_gal, n_E_gam, data, E_gam__GeV, tau_gg, tau_EBL, distmod, "E^2 dN/dE [GeV cmm2 sm1]", string_cat(outfp, filelist[k]) );
    }



    double L_gamma[n_gal];
    double spectot[n_E_gam];
    for (i = 0; i < n_gal; i++)
    {
        for (j = 0; j < n_E_gam; j++)
        {
            spectot[j] = (spec_pi[i][j] + spec_IC_1_z1[i][j] + spec_IC_2_z1[i][j] + spec_BS_1_z1[i][j] + spec_BS_2_z1[i][j] + spec_SY_1_z1[i][j] + spec_SY_2_z1[i][j] + 
                           spec_IC_1_z2[i][j] + spec_IC_2_z2[i][j] + spec_SY_1_z2[i][j] + spec_SY_2_z2[i][j]) * exp(-tau_gg[i][j]);
        }
        gso1D_data = gsl_so1D( n_E_gam, E_gam__GeV, spectot );
        L_gamma[i] = spec_integrate_gso1D_lim( gso1D_data, 0.1, 100. );
        gsl_so1D_free( gso1D_data );
    }

    write_1D_file( n_gal, L_gamma, "L_gamma__GeVsm1", string_cat(outfp, "/L_gamma.txt") );


    for (k = 0; k < 14; k++)
    {
        free2D( n_gal, speclist[k] );
    }

    free2D( n_gal, data );
    free2D( n_gal, tau_gg );
    free2D( n_gal, tau_EBL );
    free2D( n_gal, tau_FF );
    free2D( n_gal, E_gam_z0__GeV );

    free( fdata_in.xa );
    free( fdata_in.ya );
    free( fdata_in.za );

    gsl_interp_accel_free( fdata_in.xacc );
    gsl_interp_accel_free( fdata_in.yacc );

/*################################################################################################################################*/
// Calculate total emission observed
/*
    gsl_spline_object_1D gso1D_spec;

    double N_gam_emit_obs[n_gal];

    printf("%s \n", "Calculating photon luminosity...");

    #pragma omp parallel for schedule(dynamic) private( j, gso1D_spec, spec_emit )
    for (i = 0; i < n_gal; i++)
    {
        for (j = 0; j < n_T_CR; j++)
        {
            spec_emit[j] = specs_obs[i][j] + specs_casc_obs[i][j];
        }
        gso1D_spec = gsl_so1D( n_T_CR, E_CR__GeV, spec_emit );

        N_gam_emit_obs[i] = N_gam_tot_obs( gso1D_spec, 1., 100. );
    }
*/
/*################################################################################################################################*/
/*
    double L_gam_emit[n_gal];
    double N_gam_emit[n_gal];

    printf("%s \n", "Calculating total luminosity...");



    #pragma omp parallel for schedule(dynamic) private( j, gso1D_spec, spec_emit )
    for (i = 0; i < n_gal; i++)
    {
        for (j = 0; j < n_T_CR; j++)
        {
            spec_emit[j] = specs_L_emit[i][j] + specs_casc[i][j];
        }
        gso1D_spec = gsl_so1D( n_T_CR, E_CR__GeV, spec_emit );
        L_gam_emit[i] = E_gam_tot( gso1D_spec, 0.1, 100. );
        N_gam_emit[i] = N_gam_tot( gso1D_spec, 1., 100. );
    }

*/
/*################################################################################################################################*/

    printf("%s \n", "Done");

/*################################################################################################################################*/

    return 0;
}
