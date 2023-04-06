#ifndef data_calc_h
#define data_calc_h

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "gen_funcs.h"
#include "data_objects.h"

#include "inverse_Compton.h"
#include "synchrotron.h"
#include "bremsstrahlung.h"

#include "gal_rad.h"
#include "file_io.h"

int CMB_num( double T_CMB_min__K_dummy, double T_CMB_max__K, double Delta_T__K, size_t * num )
{
    double Delta = (T_CMB_max__K - T_0_CMB__K)/Delta_T__K;
    *(num) = (size_t) ceil(Delta) + 1;
    return 0;
}

int FIR_num( double T_FIR_min__K, double T_FIR_max__K, double Delta_T__K, size_t * num )
{
    double Delta = (T_FIR_max__K - T_FIR_min__K)/Delta_T__K;
    *(num) = (size_t) ceil(Delta) + 1;
    return 0;
}

int CMB_lim( double T_CMB_min__K_dummy, double T_CMB_max__K, double * lim )
{
    lim[0] = floor(T_0_CMB__K*10.)/10.;
    lim[1] = ceil(T_CMB_max__K*10.)/10.;
    return 0;
}

int FIR_lim( double T_FIR_min__K, double T_FIR_max__K, double * lim )
{
    lim[0] = floor(T_FIR_min__K);
    lim[1] = ceil(T_FIR_max__K);
    return 0;
}

data_object_3D do3D_FIELD_array( int (*FIELD_num)( double, double, double, size_t * ), int (*FIELD_lim)( double, double, double * ), double T_FIELD_min__K, double T_FIELD_max__K, double Delta_T__K )
{
    unsigned int i;

    data_object_3D do3D;
    double Delta;
    double lim[2];

    FIELD_num( T_FIELD_min__K, T_FIELD_max__K, Delta_T__K, &(do3D.nf) );
    FIELD_lim( T_FIELD_min__K, T_FIELD_max__K, lim );
    do3D.f_data = malloc(sizeof(double) * do3D.nf);
    do3D.f_data[0] = lim[0];
    do3D.f_data[do3D.nf-1] = lim[1];
//printf("%zu %le %le\n", do3D.nf, lim[0], lim[1] );
    do3D.f_lim[0] = do3D.f_data[0];
    do3D.f_lim[1] = do3D.f_data[do3D.nf-1];
//printf("%i %le %le %le\n", nT, Delta, T_CMB_array__K[0], T_CMB_array__K[nT-1] );
    Delta = (do3D.f_data[do3D.nf-1]-do3D.f_data[0])/(do3D.nf-1);

    for (i = 0; i < do3D.nf-2; i++)
    {
        do3D.f_data[i+1] = do3D.f_data[0] + (i+1) * Delta;
    }
    return do3D;
}



data_object_2D do2D_IC( 
    data_object_2D (*init_do_2D_IC_version)( double (*)(double *, double), double *, double *, double *, double *, size_t *), 
    double (*n_phot)(double *, double), 
    double T_BB__K, 
    double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], size_t n_pts[2] )
{
    return init_do_2D_IC_version( n_phot, &(T_BB__K), E_gam__GeV_lims, E_e__GeV_lims, E_phot__GeV_lims, n_pts);
}



data_object_3D do3D_IC( 
    data_object_2D (*init_do_2D_IC_version)( double (*)(double *, double), double *, double *, double *, double *, size_t *),
    double (*n_phot)(double *, double),
    int (*FIELD_num)( double, double, double, size_t * ), int (*FIELD_lim)( double, double, double * ), 
    double T_FIELD_min__K, double T_FIELD_max__K, double Delta_T__K, 
    double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], size_t n_pts[2] )
{
    unsigned int i, j, k;
    data_object_3D do3D = do3D_FIELD_array( FIELD_num, FIELD_lim, T_FIELD_min__K, T_FIELD_max__K, Delta_T__K );

    do3D.nx = n_pts[0];
    do3D.ny = n_pts[1];
    do3D.x_data = malloc(sizeof(double) * do3D.nx);
    do3D.y_data = malloc(sizeof(double) * do3D.ny);
    do3D.z_data = malloc(sizeof *do3D.z_data * do3D.nf);
    if (do3D.z_data){for (i = 0; i < do3D.nf; i++){do3D.z_data[i] = malloc(sizeof *do3D.z_data[i] * n_pts[0] * n_pts[1]);}}

    data_object_2D dummy_do2D;

//printf("%zu %le %le %le %le %le\n", do3D.nf, Delta_T__K, do3D.f_lim[0], do3D.f_lim[1], do3D.f_data[0], do3D.f_data[do3D.nf-1] );

    for (i = 0; i < do3D.nf; i++)
    {

        dummy_do2D = init_do_2D_IC_version( n_phot, &(do3D.f_data[i]), E_gam__GeV_lims, E_e__GeV_lims, E_phot__GeV_lims, n_pts);

        for (j = 0; j < do3D.nx; j++)
        {
            for (k = 0; k < do3D.ny; k++)
            {
                do3D.z_data[i][k * do3D.nx + j] = dummy_do2D.z_data[k * do3D.nx + j];
            }
        }
    }

    for (i = 0; i < do3D.nx; i++)
    {
        do3D.x_data[i] = dummy_do2D.x_data[i];
    }
    for (i = 0; i < do3D.ny; i++)
    {
        do3D.y_data[i] = dummy_do2D.y_data[i];
    }

    do3D.x_lim[0] = do3D.x_data[0];
    do3D.x_lim[1] = do3D.x_data[do3D.nx-1];
    do3D.y_lim[0] = do3D.y_data[0];
    do3D.y_lim[1] = do3D.y_data[do3D.ny-1];

    return do3D;
}




IC_object load_IC_do_files( size_t n_pts[2], double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2],
                            double T_CMB_max__K, double Delta_T_CMB__K, double T_FIR_min__K, double T_FIR_max__K, double Delta_T_FIR__K,
                            char * datadir )
{
    IC_object ICo;
    unsigned short int i, j;
//    char *datadir = "data/"

    char files_do2D_strings[2][4][24] = { { "/IC_3000_do2D.txt", "/IC_4000_do2D.txt", "/IC_7500_do2D.txt", "/IC_UV_do2D.txt" },
                                       { "/IC_3000_Gamma_do2D.txt", "/IC_4000_Gamma_do2D.txt", "/IC_7500_Gamma_do2D.txt", 
                                         "/IC_UV_Gamma_do2D.txt" } };

    char files_do3D_strings[2][2][23] = { { "/IC_CMB_do3D.txt", "/IC_FIR_do3D.txt" }, 
                                          { "/IC_CMB_Gamma_do3D.txt", "/IC_FIR_Gamma_do3D.txt" } };

    char files_IC_do2D[2][4][strlen(datadir)+23+1];
    char files_IC_do3D[2][2][strlen(datadir)+22+1];

    for (j = 0; j < 2; j++)
    {
        for (i = 0; i < 4; i++)
        {
            strcpy( files_IC_do2D[j][i], string_cat( datadir, files_do2D_strings[j][i] ) );
        }
    }

    for (j = 0; j < 2; j++)
    {
        for (i = 0; i < 2; i++)
        {
            strcpy( files_IC_do3D[j][i], string_cat( datadir, files_do3D_strings[j][i] ) );
        }
    }

    typedef double (*PhotFuncArray)(double *, double);
    PhotFuncArray dndEphot2D[] = { dndEphot_BB__cmm3GeVm1, dndEphot_BB__cmm3GeVm1, dndEphot_BB__cmm3GeVm1, dndEphot_UVMattis__cmm3GeVm1 };
    PhotFuncArray dndEphot3D[] = { dndEphot_BB__cmm3GeVm1, dndEphot_modBB__cmm3GeVm1 };
    double T_comp__K[] = { 3000., 4000., 7500., 0. };

    typedef data_object_2D (*ICFuncArray)( double (*)(double *, double), double *, double *, double *, double *, size_t *);
    ICFuncArray init_do_2D_IC_version[] = { init_do_2D_IC, init_do_2D_IC_Gamma };

    typedef int (*FIELD_numArray)( double, double, double, size_t * );
    FIELD_numArray FIELD_num[] = { CMB_num, FIR_num };
    typedef int (*FIELD_limArray)( double, double, double * );
    FIELD_limArray FIELD_lim[] = { CMB_lim, FIR_lim };

    double T_FIELD_min__K[] = { 0., T_FIR_min__K };
    double T_FIELD_max__K[] = { T_CMB_max__K, T_FIR_max__K };
    double Delta_T__K[] = { Delta_T_CMB__K, Delta_T_FIR__K };

    unsigned short int cobjint;
    char * filename;

    //IC do2D
    for (j = 0; j < 2; j++)
    {
        for (i = 0; i < 4; i++)
        {
            filename = files_IC_do2D[j][i];
            ICo.do_2D_IC[j][i] = read_do2D( filename );

            cobjint = check_do2D( n_pts[0], n_pts[1], E_gam__GeV_lims, E_e__GeV_lims, ICo.do_2D_IC[j][i] );

            if (cobjint == 1)
            {
                printf("Trying to calc/write file %s\n", filename);
                fflush(stdout);
                ICo.do_2D_IC[j][i] = init_do_2D_IC_version[j]( dndEphot2D[i], &(T_comp__K[i]), E_gam__GeV_lims, E_e__GeV_lims, E_phot__GeV_lims, n_pts );
                write_do2D( ICo.do_2D_IC[j][i], filename );
            }
            else
            {
                printf("Successfully read and imported file %s\n", filename);
                fflush(stdout);
            }
        }

        //IC do3D
        for (i = 0; i < 2; i++)
        {
            filename = files_IC_do3D[j][i];
            ICo.do_3D_IC[j][i] = read_do3D( filename );

            cobjint = check_do3D_FIELD( n_pts[0], n_pts[1], E_gam__GeV_lims, E_e__GeV_lims,
                              FIELD_num[i], FIELD_lim[i], T_FIELD_min__K[i], T_FIELD_max__K[i], Delta_T__K[i], ICo.do_3D_IC[j][i] );

            if (cobjint == 1)
            {
                printf("Trying to calc/write file %s\n", filename);
                fflush(stdout);
                ICo.do_3D_IC[j][i] = do3D_IC( init_do_2D_IC_version[i], dndEphot3D[i], FIELD_num[i], FIELD_lim[i], 
                    T_FIELD_min__K[i], T_FIELD_max__K[i], Delta_T__K[i], E_gam__GeV_lims, E_e__GeV_lims, E_phot__GeV_lims, n_pts );
                write_do3D( ICo.do_3D_IC[j][i], filename );
            }
            else
            {
                printf("Successfully read and imported file %s\n", filename);
                fflush(stdout);
            }
        }
    }

    return ICo;
}

data_object_2D load_BS_do_files( size_t n_pts[2], double E_gam__GeV_lims[2], double E_e__GeV_lims[2], char *datadir )
{
    data_object_2D do_2D_BS;

    char file_BS_do2D[1][strlen(datadir)+12+1];
    strcpy( file_BS_do2D[0], string_cat(datadir, "/BS_do2D.txt") );

    unsigned short int cobjint;
    char * filename;

    //BS
    filename = file_BS_do2D[0];
    do_2D_BS = read_do2D( filename );

    cobjint = check_do2D( n_pts[0], n_pts[1], E_gam__GeV_lims, E_e__GeV_lims, do_2D_BS );

    if (cobjint == 1)
    {
        printf("Trying to calc/write file %s\n", filename);
        fflush(stdout);
        do_2D_BS = init_do_2D_BS( E_gam__GeV_lims, E_e__GeV_lims, n_pts );;
        write_do2D( do_2D_BS, filename );
    }
    else
    {
        printf("Successfully read and imported file %s\n", filename);
        fflush(stdout);
    }

    return do_2D_BS;
}

data_object_1D load_SY_do_files( size_t n_pts[1], double x_sync_lims[2], char *datadir )
{
    data_object_1D do_1D_SY;

    char file_SY_do1D[1][strlen(datadir)+12+1];
    strcpy( file_SY_do1D[0], string_cat(datadir, "/SY_do1D.txt") );

    unsigned short int cobjint;
    char * filename;

    //SY
    filename = file_SY_do1D[0];
    do_1D_SY = read_do1D( filename );

    cobjint = check_do1D( n_pts[0], x_sync_lims, do_1D_SY );

    if (cobjint == 1)
    {
        printf("Trying to calc/write file %s\n", filename);
        fflush(stdout);
        do_1D_SY = init_do_1D_sync( x_sync_lims, n_pts[0] );
        write_do1D( do_1D_SY, filename );
    }
    else
    {
        printf("Successfully read and imported file %s\n", filename);
        fflush(stdout);
    }

    return do_1D_SY;
}

//intj is 0 for IC and 1 for IC_Gamma
gsl_spline_object_2D construct_IC_gso2D( unsigned short int intj, IC_object ICo, double * nphot_params )
{

    double C_dil_gal[5];

    //Draine 3000K
    C_dil_gal[0] = C_dil( u_rad_BB__GeVcmm3( 3000. ), L3000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
    //Draine 4000K
    C_dil_gal[1] = C_dil( u_rad_BB__GeVcmm3( 4000. ),  L4000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
    //Draine 7500K
    C_dil_gal[2] = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
    //Draine UV Mattis field
    C_dil_gal[3] = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
    //FIR field from dust
    C_dil_gal[4] = C_dil( u_rad_modBB__GeVcmm3( nphot_params[1] ),  LFIR__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );


    unsigned int i, j;

    //CMB j = 0, FIR j = 1

    double **z_data_CF = malloc(sizeof( *z_data_CF ) * 2);
    if ( z_data_CF ){ for (i = 0; i < 2; i++){ z_data_CF[i] = malloc(sizeof( * z_data_CF[i] ) * ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny); } }


    double *doub_int_array;
    gsl_spline_object_1D gso1D_FIELD_array;
    double posTFIELD;
    int intTFIELD;
    double fracTFIELD;

    for (j = 0; j < 2; j++)
    {
        doub_int_array = malloc( sizeof(double) * ICo.do_3D_IC[intj][j].nf);
        double_integer_array( 0, ICo.do_3D_IC[intj][j].nf-1, doub_int_array );

        gso1D_FIELD_array = gsl_so1D( ICo.do_3D_IC[intj][j].nf, ICo.do_3D_IC[intj][j].f_data, doub_int_array );

        posTFIELD = gsl_so1D_eval( gso1D_FIELD_array, nphot_params[j] );

        intTFIELD = (int) floor(posTFIELD);
        fracTFIELD = posTFIELD - (double) intTFIELD;

    

        for (i = 0; i < ICo.do_3D_IC[intj][j].nx * ICo.do_3D_IC[intj][j].ny; i++)
        {
            z_data_CF[j][i] = ICo.do_3D_IC[intj][j].z_data[intTFIELD][i] * fracTFIELD + 
                           ICo.do_3D_IC[intj][j].z_data[intTFIELD+1][i] * (1.-fracTFIELD);
        }
        gsl_so1D_free( gso1D_FIELD_array );
        free( doub_int_array );
    }

    double * z_data = malloc(sizeof(double) * ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny);


    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = z_data_CF[0][i] + 
                    C_dil_gal[4] * z_data_CF[1][i] + 
                    C_dil_gal[0] * ICo.do_2D_IC[intj][0].z_data[i] + 
                    C_dil_gal[1] * ICo.do_2D_IC[intj][1].z_data[i] + 
                    C_dil_gal[2] * ICo.do_2D_IC[intj][2].z_data[i] + 
                    C_dil_gal[3] * ICo.do_2D_IC[intj][3].z_data[i];
    }

    free2D( 2, z_data_CF );
    gsl_spline_object_2D gso2D_out = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );
    free( z_data );
    return gso2D_out;
}


void return_IC_gso2D( IC_object ICo, double * nphot_params, 
                      gsl_spline_object_2D* gso2D_3000, gsl_spline_object_2D* gso2D_4000, gsl_spline_object_2D* gso2D_7500, 
                      gsl_spline_object_2D* gso2D_UV, gsl_spline_object_2D* gso2D_FIR, gsl_spline_object_2D* gso2D_CMB )
{

    double C_dil_gal[5];

    //Draine 3000K
    C_dil_gal[0] = C_dil( u_rad_BB__GeVcmm3( 3000. ), L3000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
    //Draine 4000K
    C_dil_gal[1] = C_dil( u_rad_BB__GeVcmm3( 4000. ),  L4000K__Lsol( nphot_params[2] ), nphot_params[4], nphot_params[5] );
    //Draine 7500K
    C_dil_gal[2] = C_dil( u_rad_BB__GeVcmm3( 7500. ), L7500K__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
    //Draine UV Mattis field
    C_dil_gal[3] = C_dil( u_rad_UVMattis__GeVcmm3(), LUV__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );
    //FIR field from dust
    C_dil_gal[4] = C_dil( u_rad_modBB__GeVcmm3( nphot_params[1] ),  LFIR__Lsol( nphot_params[3] ), nphot_params[4], nphot_params[5] );


    unsigned int i, j;
    unsigned short int intj = 0;
    //CMB j = 0, FIR j = 1

    double **z_data_CF = malloc(sizeof( *z_data_CF ) * 2);
    if ( z_data_CF ){ for (i = 0; i < 2; i++){ z_data_CF[i] = malloc(sizeof( * z_data_CF[i] ) * ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny); } }


    double *doub_int_array;
    gsl_spline_object_1D gso1D_FIELD_array;
    double posTFIELD;
    int intTFIELD;
    double fracTFIELD;

    for (j = 0; j < 2; j++)
    {
        doub_int_array = malloc( sizeof(double) * ICo.do_3D_IC[intj][j].nf);
        double_integer_array( 0, ICo.do_3D_IC[intj][j].nf-1, doub_int_array );

        gso1D_FIELD_array = gsl_so1D( ICo.do_3D_IC[intj][j].nf, ICo.do_3D_IC[intj][j].f_data, doub_int_array );

        posTFIELD = gsl_so1D_eval( gso1D_FIELD_array, nphot_params[j] );

        intTFIELD = (int) floor(posTFIELD);
        fracTFIELD = posTFIELD - (double) intTFIELD;

    

        for (i = 0; i < ICo.do_3D_IC[intj][j].nx * ICo.do_3D_IC[intj][j].ny; i++)
        {
            z_data_CF[j][i] = ICo.do_3D_IC[intj][j].z_data[intTFIELD][i] * fracTFIELD + 
                           ICo.do_3D_IC[intj][j].z_data[intTFIELD+1][i] * (1.-fracTFIELD);
        }
        gsl_so1D_free( gso1D_FIELD_array );
        free( doub_int_array );
    }

    double * z_data = malloc(sizeof(double) * ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny);

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = z_data_CF[0][i];
    }
    *gso2D_CMB = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = C_dil_gal[4] * z_data_CF[1][i];
    }
    *gso2D_FIR = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = C_dil_gal[0] * ICo.do_2D_IC[intj][0].z_data[i];
    }
    *gso2D_3000 = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = C_dil_gal[1] * ICo.do_2D_IC[intj][1].z_data[i];
    }
    *gso2D_4000 = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = C_dil_gal[2] * ICo.do_2D_IC[intj][2].z_data[i];
    }
    *gso2D_7500 = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );

    for (i = 0; i < ICo.do_3D_IC[intj][0].nx * ICo.do_3D_IC[intj][0].ny; i++)
    {
        z_data[i] = C_dil_gal[3] * ICo.do_2D_IC[intj][3].z_data[i];
    }
    *gso2D_UV = gsl_so2D( ICo.do_3D_IC[intj][0].nx, ICo.do_3D_IC[intj][0].ny, ICo.do_3D_IC[intj][0].x_data, ICo.do_3D_IC[intj][0].y_data, z_data );


    free2D( 2, z_data_CF );
    free( z_data );
    return;
}


void write_ICspectra( unsigned int n_E_gam, double E_gam__GeV[n_E_gam], gsl_spline_object_1D qe_1_z1_so, gsl_spline_object_1D qe_2_z1_so, 
                                          gsl_spline_object_1D qe_1_z2_so, gsl_spline_object_1D qe_2_z2_so, 
                      gsl_spline_object_2D gso2D_3000, gsl_spline_object_2D gso2D_4000, gsl_spline_object_2D gso2D_7500, 
                      gsl_spline_object_2D gso2D_UV, gsl_spline_object_2D gso2D_FIR, gsl_spline_object_2D gso2D_CMB, char *filepath )
{
    gsl_spline_object_1D qe_so[4] = { qe_1_z1_so, qe_2_z1_so, qe_1_z2_so, qe_2_z2_so };
    gsl_spline_object_2D gso2D_IC[6] = { gso2D_3000, gso2D_4000, gso2D_7500, gso2D_UV, gso2D_FIR, gso2D_CMB };

    unsigned short int i, j, k;


    double **data = malloc(sizeof *data * 24);
    if (data){for (i = 0; i < 24; i++){data[i] = malloc(sizeof *data[i] * n_E_gam);}}

    for (i = 0; i < 6; i++)
    {
        for (k = 0; k < 4; k++)
        {
            for (j = 0; j < n_E_gam; j++)
            {
                data[i*4+k][j] = eps_IC_3( E_gam__GeV[j], gso2D_IC[i], qe_so[k] ) * pow(E_gam__GeV[j],2);
            }
        }
    }
    write_2D_file( 24, n_E_gam, data, "IC - prim. disc; sec. disc; prim. halo; sec. halo - 3000K; 4000K; 7500K; UV; FIR; CMB", filepath );

    free2D( 24, data );
    free( filepath );
    return;
}





#endif
