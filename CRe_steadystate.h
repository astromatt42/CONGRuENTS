#ifndef CRe_steadystate
#define CRe_steadystate

#include <cubature.h>
#include <gsl_linalg.h>

#include "gen_funcs.h"
#include "gsl_decs.h"

#include "CR_funcs.h"


struct F_int_data
{
    double E_e__GeV;
    double E_f__GeV;
    double n_H__cmm3;
    double B__G;
    double h__pc;
    int n_gso2D;
    gsl_spline_object_2D * gso_2D_radfield;
    gsl_spline_object_2D gso2D_BS;
    gsl_spline_object_1D gso_1D_Q;
    gsl_spline_object_1D gso_1D_D__cm2sm1;
    //Pointer to photon field function
    double (*n_phot)( double *, double );
    double (*E_func)( double, double, double, double );
    //Pointer to array containing the arguments for the above
    double *n_phot_params;
};



int F_Gamma_2D_2_log( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}

int F_Gamma_i0_2D_2_log( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}


int F_Gamma_ii_2D_2_log( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            fval[j * fdim + 0] = (exp(x[j*ndim+0])-exp(x[j*ndim+1])) * dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), 
                                 exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
        }
    }
    return 0;
}


int F_Q_i_log( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.gso_1D_Q, exp(x[j*ndim+0]) );
    }
    return 0;
}


int F_D_2_log( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = gsl_so1D_eval( fdata_in.gso_1D_D__cm2sm1, exp(x[j*ndim+0]) );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}



//Additional functions for the energy conserving scheme



int F_Gamma_i0_2D_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            fval[j * fdim + 0] = exp(x[j*ndim+0]) * dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), 
                                 exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
        }
    }
    return 0;
}

int F_Gamma_ii_2D_4_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            fval[j * fdim + 0] = (exp(x[j*ndim+0])-exp(x[j*ndim+1])) * dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), 
                                 exp(x[j*ndim+1]), fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
        }
    }
    return 0;
}

//ij is 0,1 and j,i is 2,3
int F_Gamma_2D_4_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    double dummy1;
    for (j = 0; j < npts; ++j)
    {
        if (exp(x[j*ndim+0]) > exp(x[j*ndim+1]))
        {
            dummy1 = dGammadlogEf_total__logGeVm1sm1( exp(x[j*ndim+0]), exp(x[j*ndim+1]), 
                                                      fdata_in.n_gso2D, fdata_in.gso_2D_radfield, fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
            fval[j * fdim + 0] = exp(x[j*ndim+0]) * dummy1;
            fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
            fval[j * fdim + 2] = exp(x[j*ndim+1]) * dummy1;
            fval[j * fdim + 3] = x[j*ndim+0] * fval[j * fdim + 2];
        }
        else
        {
            fval[j * fdim + 0] = 0.;
            fval[j * fdim + 1] = 0.;
            fval[j * fdim + 2] = 0.;
            fval[j * fdim + 3] = 0.;
        }
    }
    return 0;
}

int F_EdotDE_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.gso_1D_D__cm2sm1, exp(x[j*ndim+0]) )/pow(fdata_in.h__pc*pc__cm,2) - 
                             fdata_in.E_func( exp(x[j*ndim+0]), fdata_in.B__G, fdata_in.n_H__cmm3, fdata_in.h__pc );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}
/*
int F_DE_2_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j * fdim + 0] = exp(x[j*ndim+0]) * gsl_so1D_eval( fdata_in.gso_1D_D__cm2sm1, exp(x[j*ndim+0]) );
        fval[j * fdim + 1] = x[j*ndim+0] * fval[j * fdim + 0];
    }
    return 0;
}
*/
int F_QE2_i_log_E( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
{
    unsigned j;
    struct F_int_data fdata_in = *((struct F_int_data *)fdata);
    for (j = 0; j < npts; ++j)
    {
        fval[j] = pow(exp(x[j*ndim+0]),2) * gsl_so1D_eval( fdata_in.gso_1D_Q, exp(x[j*ndim+0]) );
    }
    return 0;
}


int CRe_steadystate_solve_3( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc, 
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1, 
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2, gsl_spline_object_1D gso_1D_Q_inject_3, 
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D, gsl_spline_object_1D * qe_3_so_1D )
{

    int solve_system( int n_E, double A[n_E*n_E], double b[n_E], double x_out[n_E])
    {
        gsl_matrix_view m_gsl = gsl_matrix_view_array( A, n_E, n_E );
        gsl_vector_view b_gsl = gsl_vector_view_array( b, n_E );
        gsl_vector *x = gsl_vector_alloc( n_E );
        int s;
        gsl_permutation * p = gsl_permutation_alloc( n_E );
        gsl_linalg_LU_decomp( &m_gsl.matrix, p, &s );
        gsl_linalg_LU_solve( &m_gsl.matrix, p, &b_gsl.vector, x );

        int i;
        for (i = 0; i < n_E; ++i)
        {   
            x_out[i] = fmax(0.,gsl_vector_get( x, i ));
        }
        gsl_permutation_free(p);
        gsl_vector_free(x);
        return 0;
    }

    double E__GeV[n_E+1];
    logspace_array( n_E+1, E_e_lim__GeV[0], E_e_lim__GeV[1], E__GeV );

    double DeltalogE = log(E_e_lim__GeV[1]/E_e_lim__GeV[0])/n_E; //log(E__GeV[1]/E__GeV[0]);

    int i,j;

    double lnE_i__GeV[n_E];
    for (i = 0; i < n_E; ++i)
    {
        lnE_i__GeV[i] = (log(E__GeV[i])+log(E__GeV[i+1]))/2.;
    }

    double E_out__GeV[n_E+2];
    for (i = 0; i < n_E; ++i)
    {
        E_out__GeV[i+1] = exp((log(E__GeV[i])+log(E__GeV[i+1]))/2.);
    }
    E_out__GeV[0] = E__GeV[0];
    E_out__GeV[n_E+1] = E__GeV[n_E];


    double xmin[1], xmax[1];
    double xmin2D[2], xmax2D[2];
    double res, abserr;
    double res2[2], abserr2[2];
    double res4[4], abserr4[4];


    struct F_int_data fdata;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.B__G = B__G;
    fdata.h__pc = h__pc;
    fdata.n_gso2D = n_gso2D;
    fdata.gso_2D_radfield = gso_2D_radfields;
    fdata.gso2D_BS = gso2D_BS;
    fdata.gso_1D_D__cm2sm1 = gso_1D_D__cm2sm1;


    //calculate losses to anything below min energy down to m_e
    double Gamma_i0[n_E];
    double Gamma_i0_prime[n_E];
    xmin2D[1] = log(m_e__GeV);
    xmax2D[1] = log(E__GeV[0]);
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_Gamma_i0_2D_2_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_i0[i] = res2[0]/DeltalogE;
        Gamma_i0_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_i0[i];
    }



    //Energy correction for intra-bin losses
    double Gamma_ii[n_E];
    double Gamma_ii_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        xmin2D[1] = log(E__GeV[i]);
        xmax2D[1] = log(E__GeV[i+1]);
        hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );
        Gamma_ii[i] = (res4[0]/DeltalogE) - (res4[2]/DeltalogE) ;
        Gamma_ii_prime[i] = (res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[0]/DeltalogE)) - 
                            (res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[2]/DeltalogE));
    }

    //transitions from bin i to bin j
//    double Gamma_ij[n_E][n_E];
//    double Gamma_ij_prime[n_E][n_E];
//    double Gamma_ji_prime[n_E][n_E];
//    double Gamma_ji[n_E][n_E];

    double **Gamma_ij = malloc(sizeof *Gamma_ij * n_E);
    if (Gamma_ij){for (i = 0; i < n_E; i++){Gamma_ij[i] = malloc(sizeof *Gamma_ij[i] * n_E);}}
    double **Gamma_ij_prime = malloc(sizeof *Gamma_ij_prime * n_E);
    if (Gamma_ij_prime){for (i = 0; i < n_E; i++){Gamma_ij_prime[i] = malloc(sizeof *Gamma_ij_prime[i] * n_E);}}
    double **Gamma_ji_prime = malloc(sizeof *Gamma_ji_prime * n_E);
    if (Gamma_ji_prime){for (i = 0; i < n_E; i++){Gamma_ji_prime[i] = malloc(sizeof *Gamma_ji_prime[i] * n_E);}}
    double **Gamma_ji = malloc(sizeof *Gamma_ji * n_E);
    if (Gamma_ji){for (i = 0; i < n_E; i++){Gamma_ji[i] = malloc(sizeof *Gamma_ji[i] * n_E);}}

    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        for (j = 0; j < i; ++j)
        {
            xmin2D[1] = log(E__GeV[j]);
            xmax2D[1] = log(E__GeV[j+1]);

            hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );

            Gamma_ij[i][j] = res4[0]/DeltalogE;
            Gamma_ij_prime[i][j] = res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ij[i][j];
            Gamma_ji[i][j] = res4[2]/DeltalogE;
            Gamma_ji_prime[i][j] = res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ji[i][j];
        }
        for (j = i; j < n_E; ++j)
        {
            Gamma_ij[i][j] = 0.;
            Gamma_ij_prime[i][j] = 0.;
            Gamma_ji[i][j] = 0.;
            Gamma_ji_prime[i][j] = 0.;
        }
    }

    //add up the transitions out of bin i
    double Gamma_i[n_E];
    double Gamma_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        Gamma_i[i] = Gamma_i0[i] + Gamma_ii[i];
        Gamma_i_prime[i] = Gamma_i0_prime[i] + Gamma_ii_prime[i];
        for (j = 0; j < i; ++j)
        {
            Gamma_i[i] += Gamma_ij[i][j];
            Gamma_i_prime[i] += Gamma_ij_prime[i][j];
        }
    }

    double Edot_i[n_E+1];
    if (structure == 1)
    {
        fdata.E_func = dEdtm1_total_disc__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_disc__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );//E__GeV[i];
        }
    }
    else if (structure == 2)
    {
        fdata.E_func = dEdtm1_total_halo__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_halo__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );//E__GeV[i];
        }
    }
    else
    {
        printf("No valid structure specified in CRe steady state solver!");
    }

    //set diffusion D_i
    double D_i[n_E];
    double D_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_EdotDE_2_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        D_i[i] = res2[0]/DeltalogE;
        D_i_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * D_i[i];
    }



    //for readability and debugging, we populate the matrix separately before changing to a 1D array
    double **M = malloc(sizeof *M * n_E);
    if (M){for (i = 0; i < n_E; i++){M[i] = malloc(sizeof *M[i] * n_E);}}

    //Second order scheme
    for (i = 0; i < n_E-2; ++i)
    {

        for (j = 0; j < i-1; ++j)
        {
            M[i][j] = 0.;
        }

        if (i > 0)
        {
            M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;
        }
        
        if (i == 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] 
                      - Gamma_ji_prime[i+1][i]/2. - Edot_i[i]/(2.*DeltalogE) + D_i_prime[i] + Gamma_i_prime[i];
        }
        else if (i > 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;
        }


        M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                    + Gamma_ji[i+1][i] - Gamma_i_prime[i]/2. - Gamma_ji_prime[i+2][i]/2.;

        if (i == n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2.;
        }
        else if (i < n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2. - Gamma_ji_prime[i+3][i]/2.;
        }

        for (j = i+3; j < n_E; ++j)
        {
            if (j == n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2.;
            }
            else if (j < n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2. - Gamma_ji_prime[j+1][i]/2.;
            }

        }
    }

    i = n_E-2;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;

    M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                + Gamma_ji[i+1][i] - Gamma_i_prime[i];

    i = n_E-1;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i];



    //set injection Q_i
    double Q_i_1[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_1;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_1[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }

    double Q_i_2[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_2;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_2[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }

    double Q_i_3[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_3;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_3[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }



    double q_e_1[n_E+2];

    if (Q_i_1[n_E-1] == 0.)
    {
        int n_Estar1 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_1[i] == 0.)
            {
                n_Estar1 = n_Estar1 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_1star[n_Estar1];
        for (i = 0; i < n_Estar1; ++i)
        {
            Q_i_1star[i] = Q_i_1[i];
        }
        double *A1 = malloc(sizeof A1 * n_Estar1*n_Estar1);
        for (i = 0; i < n_Estar1; ++i)
        {
            for (j = 0; j < n_Estar1; ++j)
            {
                A1[i*n_Estar1+j] = M[i][j];
            }
        }
        double x_out_1[n_Estar1];
        solve_system( n_Estar1, A1, Q_i_1star, x_out_1);
        for (i = 0; i < n_Estar1; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar1; i < n_E; ++i)
        {   
            q_e_1[i+1] = 0.;
        }
        free( A1 );

    }
    else
    {
        double *A1 = malloc(sizeof A1 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A1[i*n_E+j] = M[i][j];
            }
        }
        double x_out_1[n_E];
        solve_system( n_E, A1, Q_i_1, x_out_1);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        free( A1 );
    }


    double q_e_2[n_E+2];

    if (Q_i_2[n_E-1] == 0.)
    {
        int n_Estar2 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_2[i] == 0.)
            {
                n_Estar2 = n_Estar2 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_2star[n_Estar2];
        for (i = 0; i < n_Estar2; ++i)
        {
            Q_i_2star[i] = Q_i_2[i];
        }
        double *A2 = malloc(sizeof A2 * n_Estar2*n_Estar2);
        for (i = 0; i < n_Estar2; ++i)
        {
            for (j = 0; j < n_Estar2; ++j)
            {
                A2[i*n_Estar2+j] = M[i][j];
            }
        }
        double x_out_2[n_Estar2];
        solve_system( n_Estar2, A2, Q_i_2star, x_out_2);
        for (i = 0; i < n_Estar2; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar2; i < n_E; ++i)
        {   
            q_e_2[i+1] = 0.;
        }
        free( A2 );
    }
    else
    {
        double *A2 = malloc(sizeof A2 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A2[i*n_E+j] = M[i][j];
            }
        }
        double x_out_2[n_E];
        solve_system( n_E, A2, Q_i_2, x_out_2);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        free( A2 );
    }

    double q_e_3[n_E+2];

    if (Q_i_3[n_E-1] == 0.)
    {
        int n_Estar3 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_3[i] == 0.)
            {
                n_Estar3 = n_Estar3 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_3star[n_Estar3];
        for (i = 0; i < n_Estar3; ++i)
        {
            Q_i_3star[i] = Q_i_3[i];
        }
        double *A3 = malloc(sizeof A3 * n_Estar3*n_Estar3);
        for (i = 0; i < n_Estar3; ++i)
        {
            for (j = 0; j < n_Estar3; ++j)
            {
                A3[i*n_Estar3+j] = M[i][j];
            }
        }
        double x_out_3[n_Estar3];
        solve_system( n_Estar3, A3, Q_i_3star, x_out_3);
        for (i = 0; i < n_Estar3; ++i)
        {   
            q_e_3[i+1] = x_out_3[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar3; i < n_E; ++i)
        {
            q_e_3[i+1] = 0.;
        }
        free( A3 );
    }
    else
    {
        double *A3 = malloc(sizeof A3 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A3[i*n_E+j] = M[i][j];
            }
        }
        double x_out_3[n_E];
        solve_system( n_E, A3, Q_i_3, x_out_3);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_3[i+1] = x_out_3[i]/E_out__GeV[i+1];
        }
        free( A3 );
    }
 
    q_e_1[0] = fmax(0.,exp( ((log(q_e_1[2])-log(q_e_1[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_1[1]) ));
    q_e_1[n_E+1] = fmax(0.,exp( ((log(q_e_1[n_E])-log(q_e_1[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_1[n_E]) ));
    q_e_2[0] = fmax(0.,exp( ((log(q_e_2[2])-log(q_e_2[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_2[1]) ));
    q_e_2[n_E+1] = fmax(0.,exp( ((log(q_e_2[n_E])-log(q_e_2[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_2[n_E]) ));
    q_e_3[0] = fmax(0.,exp( ((log(q_e_3[2])-log(q_e_3[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_3[1]) ));
    q_e_3[n_E+1] = fmax(0.,exp( ((log(q_e_3[n_E])-log(q_e_3[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_3[n_E]) ));


    *qe_1_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_1 );
    *qe_2_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_2 );
    *qe_3_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_3 );

    free2D( n_E, Gamma_ij );
    free2D( n_E, Gamma_ij_prime );
    free2D( n_E, Gamma_ji );
    free2D( n_E, Gamma_ji_prime );

    free2D( n_E, M );

    return 0;

}


int CRe_steadystate_solve( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc, 
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1, 
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2, 
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D )
{

    int solve_system( int n_E, double A[n_E*n_E], double b[n_E], double x_out[n_E])
    {
        gsl_matrix_view m_gsl = gsl_matrix_view_array( A, n_E, n_E );
        gsl_vector_view b_gsl = gsl_vector_view_array( b, n_E );
        gsl_vector *x = gsl_vector_alloc( n_E );
        int s;
        gsl_permutation * p = gsl_permutation_alloc( n_E );
        gsl_linalg_LU_decomp( &m_gsl.matrix, p, &s );
        gsl_linalg_LU_solve( &m_gsl.matrix, p, &b_gsl.vector, x );

        int i;
        for (i = 0; i < n_E; ++i)
        {   
            x_out[i] = fmax(0.,gsl_vector_get( x, i ));
        }
        gsl_permutation_free(p);
        gsl_vector_free(x);
        return 0;
    }

    double E__GeV[n_E+1];
    logspace_array( n_E+1, E_e_lim__GeV[0], E_e_lim__GeV[1], E__GeV );

    double DeltalogE = log(E_e_lim__GeV[1]/E_e_lim__GeV[0])/n_E; //log(E__GeV[1]/E__GeV[0]);

    int i,j;

    double lnE_i__GeV[n_E];
    for (i = 0; i < n_E; ++i)
    {
        lnE_i__GeV[i] = (log(E__GeV[i])+log(E__GeV[i+1]))/2.;
    }

    double E_out__GeV[n_E+2];
    for (i = 0; i < n_E; ++i)
    {
        E_out__GeV[i+1] = exp((log(E__GeV[i])+log(E__GeV[i+1]))/2.);
    }
    E_out__GeV[0] = E__GeV[0];
    E_out__GeV[n_E+1] = E__GeV[n_E];


    double xmin[1], xmax[1];
    double xmin2D[2], xmax2D[2];
    double res, abserr;
    double res2[2], abserr2[2];
    double res4[4], abserr4[4];


    struct F_int_data fdata;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.B__G = B__G;
    fdata.h__pc = h__pc;
    fdata.n_gso2D = n_gso2D;
    fdata.gso_2D_radfield = gso_2D_radfields;
    fdata.gso2D_BS = gso2D_BS;
    fdata.gso_1D_D__cm2sm1 = gso_1D_D__cm2sm1;


    //calculate losses to anything below min energy down to m_e
    double Gamma_i0[n_E];
    double Gamma_i0_prime[n_E];
    xmin2D[1] = log(m_e__GeV);
    xmax2D[1] = log(E__GeV[0]);
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_Gamma_i0_2D_2_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_i0[i] = res2[0]/DeltalogE;
        Gamma_i0_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_i0[i];
    }



    //Energy correction for intra-bin losses
    double Gamma_ii[n_E];
    double Gamma_ii_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        xmin2D[1] = log(E__GeV[i]);
        xmax2D[1] = log(E__GeV[i+1]);
        hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );
        Gamma_ii[i] = (res4[0]/DeltalogE) - (res4[2]/DeltalogE) ;
        Gamma_ii_prime[i] = (res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[0]/DeltalogE)) - 
                            (res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * (res4[2]/DeltalogE));
    }

    //transitions from bin i to bin j
//    double Gamma_ij[n_E][n_E];
//    double Gamma_ij_prime[n_E][n_E];
//    double Gamma_ji_prime[n_E][n_E];
//    double Gamma_ji[n_E][n_E];

    double **Gamma_ij = malloc(sizeof *Gamma_ij * n_E);
    if (Gamma_ij){for (i = 0; i < n_E; i++){Gamma_ij[i] = malloc(sizeof *Gamma_ij[i] * n_E);}}
    double **Gamma_ij_prime = malloc(sizeof *Gamma_ij_prime * n_E);
    if (Gamma_ij_prime){for (i = 0; i < n_E; i++){Gamma_ij_prime[i] = malloc(sizeof *Gamma_ij_prime[i] * n_E);}}
    double **Gamma_ji_prime = malloc(sizeof *Gamma_ji_prime * n_E);
    if (Gamma_ji_prime){for (i = 0; i < n_E; i++){Gamma_ji_prime[i] = malloc(sizeof *Gamma_ji_prime[i] * n_E);}}
    double **Gamma_ji = malloc(sizeof *Gamma_ji * n_E);
    if (Gamma_ji){for (i = 0; i < n_E; i++){Gamma_ji[i] = malloc(sizeof *Gamma_ji[i] * n_E);}}

    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        for (j = 0; j < i; ++j)
        {
            xmin2D[1] = log(E__GeV[j]);
            xmax2D[1] = log(E__GeV[j+1]);

            hcubature_v( 4, F_Gamma_2D_4_log_E, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res4, abserr4 );

            Gamma_ij[i][j] = res4[0]/DeltalogE;
            Gamma_ij_prime[i][j] = res4[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ij[i][j];
            Gamma_ji[i][j] = res4[2]/DeltalogE;
            Gamma_ji_prime[i][j] = res4[3]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_ji[i][j];
        }
        for (j = i; j < n_E; ++j)
        {
            Gamma_ij[i][j] = 0.;
            Gamma_ij_prime[i][j] = 0.;
            Gamma_ji[i][j] = 0.;
            Gamma_ji_prime[i][j] = 0.;
        }
    }

    //add up the transitions out of bin i
    double Gamma_i[n_E];
    double Gamma_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        Gamma_i[i] = Gamma_i0[i] + Gamma_ii[i];
        Gamma_i_prime[i] = Gamma_i0_prime[i] + Gamma_ii_prime[i];
        for (j = 0; j < i; ++j)
        {
            Gamma_i[i] += Gamma_ij[i][j];
            Gamma_i_prime[i] += Gamma_ij_prime[i][j];
        }
    }

    double Edot_i[n_E+1];
    if (structure == 1)
    {
        fdata.E_func = dEdtm1_total_disc__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_disc__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );//E__GeV[i];
        }
    }
    else if (structure == 2)
    {
        fdata.E_func = dEdtm1_total_halo__GeVsm1;
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_halo__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc );//E__GeV[i];
        }
    }
    else
    {
        printf("No valid structure specified in CRe steady state solver!");
    }

    //set diffusion D_i
    double D_i[n_E];
    double D_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_EdotDE_2_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        D_i[i] = res2[0]/DeltalogE;
        D_i_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * D_i[i];
    }



    //for readability and debugging, we populate the matrix separately before changing to a 1D array
    double **M = malloc(sizeof *M * n_E);
    if (M){for (i = 0; i < n_E; i++){M[i] = malloc(sizeof *M[i] * n_E);}}

    //Second order scheme
    for (i = 0; i < n_E-2; ++i)
    {

        for (j = 0; j < i-1; ++j)
        {
            M[i][j] = 0.;
        }

        if (i > 0)
        {
            M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;
        }
        
        if (i == 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] 
                      - Gamma_ji_prime[i+1][i]/2. - Edot_i[i]/(2.*DeltalogE) + D_i_prime[i] + Gamma_i_prime[i];
        }
        else if (i > 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;
        }


        M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                    + Gamma_ji[i+1][i] - Gamma_i_prime[i]/2. - Gamma_ji_prime[i+2][i]/2.;

        if (i == n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2.;
        }
        else if (i < n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2. - Gamma_ji_prime[i+3][i]/2.;
        }

        for (j = i+3; j < n_E; ++j)
        {
            if (j == n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2.;
            }
            else if (j < n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2. - Gamma_ji_prime[j+1][i]/2.;
            }

        }
    }

    i = n_E-2;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;

    M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                + Gamma_ji[i+1][i] - Gamma_i_prime[i];

    i = n_E-1;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i];



    //set injection Q_i
    double Q_i_1[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_1;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_1[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }

    double Q_i_2[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_2;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_QE2_i_log_E, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_2[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }



    double q_e_1[n_E+2];

    if (Q_i_1[n_E-1] == 0.)
    {
        int n_Estar1 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_1[i] == 0.)
            {
                n_Estar1 = n_Estar1 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_1star[n_Estar1];
        for (i = 0; i < n_Estar1; ++i)
        {
            Q_i_1star[i] = Q_i_1[i];
        }
        double *A1 = malloc(sizeof A1 * n_Estar1*n_Estar1);
        for (i = 0; i < n_Estar1; ++i)
        {
            for (j = 0; j < n_Estar1; ++j)
            {
                A1[i*n_Estar1+j] = M[i][j];
            }
        }
        double x_out_1[n_Estar1];
        solve_system( n_Estar1, A1, Q_i_1star, x_out_1);
        for (i = 0; i < n_Estar1; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar1; i < n_E; ++i)
        {   
            q_e_1[i+1] = 0.;
        }
        free( A1 );

    }
    else
    {
        double *A1 = malloc(sizeof A1 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A1[i*n_E+j] = M[i][j];
            }
        }
        double x_out_1[n_E];
        solve_system( n_E, A1, Q_i_1, x_out_1);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        free( A1 );
    }


    double q_e_2[n_E+2];

    if (Q_i_2[n_E-1] == 0.)
    {
        int n_Estar2 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_2[i] == 0.)
            {
                n_Estar2 = n_Estar2 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_2star[n_Estar2];
        for (i = 0; i < n_Estar2; ++i)
        {
            Q_i_2star[i] = Q_i_2[i];
        }
        double *A2 = malloc(sizeof A2 * n_Estar2*n_Estar2);
        for (i = 0; i < n_Estar2; ++i)
        {
            for (j = 0; j < n_Estar2; ++j)
            {
                A2[i*n_Estar2+j] = M[i][j];
            }
        }
        double x_out_2[n_Estar2];
        solve_system( n_Estar2, A2, Q_i_2star, x_out_2);
        for (i = 0; i < n_Estar2; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar2; i < n_E; ++i)
        {   
            q_e_2[i+1] = 0.;
        }
        free( A2 );
    }
    else
    {
        double *A2 = malloc(sizeof A2 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A2[i*n_E+j] = M[i][j];
            }
        }
        double x_out_2[n_E];
        solve_system( n_E, A2, Q_i_2, x_out_2);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        free( A2 );
    }

 
    q_e_1[0] = fmax(0.,exp( ((log(q_e_1[2])-log(q_e_1[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_1[1]) ));
    q_e_1[n_E+1] = fmax(0.,exp( ((log(q_e_1[n_E])-log(q_e_1[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_1[n_E]) ));
    q_e_2[0] = fmax(0.,exp( ((log(q_e_2[2])-log(q_e_2[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_2[1]) ));
    q_e_2[n_E+1] = fmax(0.,exp( ((log(q_e_2[n_E])-log(q_e_2[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_2[n_E]) ));

    *qe_1_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_1 );
    *qe_2_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_2 );

    free2D( n_E, Gamma_ij );
    free2D( n_E, Gamma_ij_prime );
    free2D( n_E, Gamma_ji );
    free2D( n_E, Gamma_ji_prime );

    free2D( n_E, M );

    return 0;

}


int CRe_steadystate_solve_number( int structure, double E_e_lim__GeV[2], int n_E, double n_H__cmm3, double B__G, double h__pc, 
    unsigned int n_gso2D, gsl_spline_object_2D * gso_2D_radfields, gsl_spline_object_2D gso2D_BS, gsl_spline_object_1D gso_1D_D__cm2sm1, 
    gsl_spline_object_1D gso_1D_Q_inject_1, gsl_spline_object_1D gso_1D_Q_inject_2, 
    gsl_spline_object_1D * qe_1_so_1D, gsl_spline_object_1D * qe_2_so_1D )
{

    int solve_system( int n_E, double A[n_E*n_E], double b[n_E], double x_out[n_E])
    {
        gsl_matrix_view m_gsl = gsl_matrix_view_array( A, n_E, n_E );
        gsl_vector_view b_gsl = gsl_vector_view_array( b, n_E );
        gsl_vector *x = gsl_vector_alloc( n_E );
        int s;
        gsl_permutation * p = gsl_permutation_alloc( n_E );
        gsl_linalg_LU_decomp( &m_gsl.matrix, p, &s );
        gsl_linalg_LU_solve( &m_gsl.matrix, p, &b_gsl.vector, x );

        int i;
        for (i = 0; i < n_E; ++i)
        {   
            x_out[i] = fmax(0.,gsl_vector_get( x, i ));
        }
        gsl_permutation_free(p);
        gsl_vector_free(x);
        return 0;
    }

    double E__GeV[n_E+1];
    logspace_array( n_E+1, E_e_lim__GeV[0], E_e_lim__GeV[1], E__GeV );

    double DeltalogE = log(E_e_lim__GeV[1]/E_e_lim__GeV[0])/n_E;

    int i,j;

    double lnE_i__GeV[n_E];
    for (i = 0; i < n_E; ++i)
    {
        lnE_i__GeV[i] = (log(E__GeV[i])+log(E__GeV[i+1]))/2.;
    }

    double E_out__GeV[n_E+2];
    for (i = 0; i < n_E; ++i)
    {
        E_out__GeV[i+1] = exp((log(E__GeV[i])+log(E__GeV[i+1]))/2.);
    }
    E_out__GeV[0] = E__GeV[0];
    E_out__GeV[n_E+1] = E__GeV[n_E];


    double xmin[1], xmax[1];
    double xmin2D[2], xmax2D[2];
    double res, abserr;
    double res2[2], abserr2[2];


    struct F_int_data fdata;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.n_gso2D = n_gso2D;
    fdata.gso_2D_radfield = gso_2D_radfields;
    fdata.gso2D_BS = gso2D_BS;
    fdata.gso_1D_D__cm2sm1 = gso_1D_D__cm2sm1;

    //calculate losses to anything below min energy down to m_e
    double Gamma_i0[n_E];
    double Gamma_i0_prime[n_E];
    xmin2D[1] = log(m_e__GeV);
    xmax2D[1] = log(E__GeV[0]);
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);

        hcubature_v( 2, F_Gamma_i0_2D_2_log, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_i0[i] = res2[0]/DeltalogE;
        Gamma_i0_prime[i] = res2[1]/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * Gamma_i0[i];
    }

/*
    //Energy correction for intra-bin losses
    double Gamma_ii[n_E];
    double Gamma_ii_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        xmin2D[1] = log(E__GeV[i]);
        xmax2D[1] = log(E__GeV[i+1]);
        hcubature_v( 2, F_Gamma_ii_2D_2_log, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        Gamma_ii[i] = res2[0]/DeltalogE / E_out__GeV[1+i];
        Gamma_ii_prime[i] = res2[1]/pow(DeltalogE,2) / E_out__GeV[1+i] - lnE_i__GeV[i]/DeltalogE * Gamma_ii[i];
    }
*/
    //transitions from bin i to bin j
//    double Gamma_ji[n_E][n_E]; // = Gamma_ij
//    double Gamma_ji_prime[n_E][n_E];

    double **Gamma_ji_prime = malloc(sizeof *Gamma_ji_prime * n_E);
    if (Gamma_ji_prime){for (i = 0; i < n_E; i++){Gamma_ji_prime[i] = malloc(sizeof *Gamma_ji_prime[i] * n_E);}}
    double **Gamma_ji = malloc(sizeof *Gamma_ji * n_E);
    if (Gamma_ji){for (i = 0; i < n_E; i++){Gamma_ji[i] = malloc(sizeof *Gamma_ji[i] * n_E);}}

    for (i = 0; i < n_E; ++i)
    {
        xmin2D[0] = log(E__GeV[i]);
        xmax2D[0] = log(E__GeV[i+1]);
        for (j = 0; j < i; ++j)
        {
            xmin2D[1] = log(E__GeV[j]);
            xmax2D[1] = log(E__GeV[j+1]);
            hcubature_v( 2, F_Gamma_2D_2_log, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
            Gamma_ji[i][j] = res2[0]/log(E__GeV[j+1]/E__GeV[j]);
            Gamma_ji_prime[i][j] = res2[1]/pow(log(E__GeV[j+1]/E__GeV[j]),2) - lnE_i__GeV[j]/log(E__GeV[j+1]/E__GeV[j]) * Gamma_ji[i][j];
        }
        for (j = i; j < n_E; ++j)
        {
            Gamma_ji[i][j] = 0.;
            Gamma_ji_prime[i][j] = 0.;
        }
    }

    //add up the transitions out of bin i
    double Gamma_i[n_E];
    double Gamma_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        Gamma_i[i] = Gamma_i0[i];// + Gamma_ii[i];
        Gamma_i_prime[i] = Gamma_i0_prime[i];// + Gamma_ii_prime[i];
        for (j = 0; j < i; ++j)
        {
            Gamma_i[i] += Gamma_ji[i][j];
            Gamma_i_prime[i] += Gamma_ji_prime[i][j];
        }
    }


    //set diffusion D_i
    double D_i[n_E];
    double D_i_prime[n_E];
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 2, F_D_2_log, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, res2, abserr2 );
        D_i[i] = (res2[0]/pow(h__pc*pc__cm,2))/DeltalogE;
        D_i_prime[i] = (res2[1]/pow(h__pc*pc__cm,2))/pow(DeltalogE,2) - lnE_i__GeV[i]/DeltalogE * D_i[i];
    }


    double Edot_i[n_E+1];
    if (structure == 1)
    {
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_disc__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc )/E__GeV[i];
        }
    }
    else if (structure == 2)
    {
        for (i = 0; i < n_E+1; ++i)
        {
            Edot_i[i] = dEdtm1_total_halo__GeVsm1( E__GeV[i], B__G, n_H__cmm3, h__pc )/E__GeV[i];
        }
    }
    else
    {
        printf("No valid structure specified in CRe steady state solver!");
    }


    //for readability and debugging, we populate the matrix separately before changing to a 1D array
    double **M = malloc(sizeof *M * n_E);
    if (M){for (i = 0; i < n_E; i++){M[i] = malloc(sizeof *M[i] * n_E);}}

    //Second order scheme
    for (i = 0; i < n_E-2; ++i)
    {

        for (j = 0; j < i-1; ++j)
        {
            M[i][j] = 0.;
        }

        if (i > 0)
        {
            M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;
        }
        
        if (i == 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] 
                      - Gamma_ji_prime[i+1][i]/2. - Edot_i[i]/(2.*DeltalogE) + D_i_prime[i] + Gamma_i_prime[i];
        }
        else if (i > 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;
        }

        M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                    + Gamma_ji[i+1][i] - Gamma_i_prime[i]/2. - Gamma_ji_prime[i+2][i]/2.;

        if (i == n_E-3)
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2.;
        }
        else
        {
            M[i][i+2] = - Edot_i[i+1]/(4.*DeltalogE) + Gamma_ji[i+2][i] + Gamma_ji_prime[i+1][i]/2. - Gamma_ji_prime[i+3][i]/2.;
        }
        

        for (j = i+3; j < n_E; ++j)
        {
            if (j == n_E-1)
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2.;
            }
            else
            {
                M[i][j] = Gamma_ji[j][i] + Gamma_ji_prime[j-1][i]/2. - Gamma_ji_prime[j+1][i]/2.;
            }
        }
    }

    i = n_E-2;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i] - Gamma_ji_prime[i+1][i]/2.;

    M[i][i+1] = - Edot_i[i+1]/DeltalogE + Edot_i[i]/(4.*DeltalogE) - D_i_prime[i]/2. 
                + Gamma_ji[i+1][i] - Gamma_i_prime[i];

    i = n_E-1;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = - Edot_i[i]/(4.*DeltalogE) + D_i_prime[i]/2. + Gamma_i_prime[i]/2.;

    M[i][i] = Edot_i[i]/DeltalogE + Edot_i[i+1]/(4.*DeltalogE) - D_i[i] - Gamma_i[i];


    //set injection Q_i
    double Q_i_1[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_1;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_Q_i_log, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_1[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }

    double Q_i_2[n_E];
    fdata.gso_1D_Q = gso_1D_Q_inject_2;
    for (i = 0; i < n_E; ++i)
    {
        xmin[0] = log(E__GeV[i]);
        xmax[0] = log(E__GeV[i+1]);
        hcubature_v( 1, F_Q_i_log, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        Q_i_2[i] = -1.*res/DeltalogE; //-ve as on RHS of Eqn in linalg system
    }



    double q_e_1[n_E+2];

    if (Q_i_1[n_E-1] == 0.)
    {
        int n_Estar1 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_1[i] == 0.)
            {
                n_Estar1 = n_Estar1 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_1star[n_Estar1];
        for (i = 0; i < n_Estar1; ++i)
        {
            Q_i_1star[i] = Q_i_1[i];
        }
        double *A1 = malloc(sizeof A1 * n_Estar1*n_Estar1);
        for (i = 0; i < n_Estar1; ++i)
        {
            for (j = 0; j < n_Estar1; ++j)
            {
                A1[i*n_Estar1+j] = M[i][j];
            }
        }
        double x_out_1[n_Estar1];
        solve_system( n_Estar1, A1, Q_i_1star, x_out_1);
        for (i = 0; i < n_Estar1; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar1; i < n_E; ++i)
        {   
            q_e_1[i+1] = 0.;
        }

        free( A1 );

    }
    else
    {
        double *A1 = malloc(sizeof A1 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A1[i*n_E+j] = M[i][j];
            }
        }
        double x_out_1[n_E];
        solve_system( n_E, A1, Q_i_1, x_out_1);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_1[i+1] = x_out_1[i]/E_out__GeV[i+1];
        }
        free( A1 );
    }


    double q_e_2[n_E+2];

    if (Q_i_2[n_E-1] == 0.)
    {
        int n_Estar2 = n_E - 1;
        for (i = 0; i < n_E-2; ++i)
        {
            if (Q_i_2[i] == 0.)
            {
                n_Estar2 = n_Estar2 - 1;
            }
            else
            {
                break;
            }
        }

        double Q_i_2star[n_Estar2];
        for (i = 0; i < n_Estar2; ++i)
        {
            Q_i_2star[i] = Q_i_2[i];
        }
        double *A2 = malloc(sizeof A2 * n_Estar2*n_Estar2);
        for (i = 0; i < n_Estar2; ++i)
        {
            for (j = 0; j < n_Estar2; ++j)
            {
                A2[i*n_Estar2+j] = M[i][j];
            }
        }
        double x_out_2[n_Estar2];
        solve_system( n_Estar2, A2, Q_i_2star, x_out_2);
        for (i = 0; i < n_Estar2; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        for (i = n_Estar2; i < n_E; ++i)
        {   
            q_e_2[i+1] = 0.;
        }
        free( A2 );
    }
    else
    {
        double *A2 = malloc(sizeof A2 * n_E*n_E);
        for (i = 0; i < n_E; ++i)
        {
            for (j = 0; j < n_E; ++j)
            {
                A2[i*n_E+j] = M[i][j];
            }
        }
        double x_out_2[n_E];
        solve_system( n_E, A2, Q_i_2, x_out_2);
        for (i = 0; i < n_E; ++i)
        {   
            q_e_2[i+1] = x_out_2[i]/E_out__GeV[i+1];
        }
        free( A2 );
    }
 

    q_e_1[0] = fmax(0.,exp( ((log(q_e_1[2])-log(q_e_1[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_1[1]) ));
    q_e_1[n_E+1] = fmax(0.,exp( ((log(q_e_1[n_E])-log(q_e_1[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_1[n_E]) ));
    q_e_2[0] = fmax(0.,exp( ((log(q_e_2[2])-log(q_e_2[1]))/(log(E_out__GeV[2])-log(E_out__GeV[1]))) * (log(E_out__GeV[0]) - log(E_out__GeV[1])) + log(q_e_2[1]) ));
    q_e_2[n_E+1] = fmax(0.,exp( ((log(q_e_2[n_E])-log(q_e_2[n_E-1]))/(log(E_out__GeV[n_E])-log(E_out__GeV[n_E-1]))) * (log(E_out__GeV[n_E+1]) - log(E_out__GeV[n_E])) + log(q_e_2[n_E]) ));

    *qe_1_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_1 );
    *qe_2_so_1D = gsl_so1D( n_E+2, E_out__GeV, q_e_2 );

    free2D( n_E, Gamma_ji );
    free2D( n_E, Gamma_ji_prime );

    free2D( n_E, M );

    return 0;

}


/*
    //First order scheme below, comment out above scheme
    for (i = 0; i < n_E-2; ++i)
    {

        for (j = 0; j < i-1; ++j)
        {
            M[i][j] = 0.;
        }

        if (i > 0)
        {
            M[i][i-1] = 0.;
        }
        
        if (i == 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE - D_i[i] - Gamma_i[i];
        }
        else if (i > 0)
        {
            M[i][i] = Edot_i[i]/DeltalogE - D_i[i] - Gamma_i[i];
        }

        M[i][i+1] = - Edot_i[i+1]/DeltalogE + Gamma_ji[i+1][i];

        M[i][i+2] = Gamma_ji[i+2][i];
        

        for (j = i+3; j < n_E; ++j)
        {
            M[i][j] = Gamma_ji[j][i];

        }
    }

    i = n_E-2;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = 0.;

    M[i][i] = Edot_i[i]/DeltalogE - D_i[i] - Gamma_i[i];

    M[i][i+1] = - Edot_i[i+1]/DeltalogE + Gamma_ji[i+1][i];

    i = n_E-1;

    for (j = 0; j < i-1; ++j)
    {
        M[i][j] = 0.;
    }

    M[i][i-1] = 0.;

    M[i][i] = Edot_i[i]/DeltalogE - D_i[i] - Gamma_i[i];

*/










#endif
