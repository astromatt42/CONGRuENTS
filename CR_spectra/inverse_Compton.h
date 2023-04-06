/*
 * This file contains the full set of routines that allow the calculation of inverse Compton emission of an electron of energy E_e__GeV
 * on a photon of energy E_phot__GeV, producing a photon of energy E_gam__GeV.
 * Permitted range of inputs:

1e-3 GeV <= E_e__GeV <= 100 PeV

 * Author: Matt Roth
 * Purpose: To show a comment that spans multiple lines.
 * Language:  C

 * This code was tested with GNU Scientific Library version 2.6
 * https://www.gnu.org/software/gsl/doc/html/index.html

 */


#ifndef inverse_Compton_h
#define inverse_Compton_h

//IO for testing, don't need this after
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cubature.h>
#include <gsl_spline2d.h>

#include "math_funcs.h"
#include "physical_constants.h"
#include "gsl_decs.h"
#include "data_objects.h"


struct hcub_data_IC 
{
    double E_e__GeV;
    double E_gam__GeV;
    double E_f__GeV;
    double Delta_E__GeV;
    //Pointer to photon field function
    double (*n_phot)( double *, double );
    //Pointer to array containing the arguments for the above
    double *n_phot_params;
    gsl_spline_object_2D gso_2D_so;
};


//This tests the IC interpolation object energies
void test_IC_interp_obj( double E_e__GeV_lims[2], 
                           gsl_spline_object_2D gso_2D_so_IC_spec_low, gsl_spline_object_2D gso_2D_so_IC_spec_high, 
                           gsl_spline_object_2D gso_2D_so_IC_Gamma, gsl_spline_object_1D qe_so_1D )
{

    struct test_IC 
    {
        double E_e__GeV;
        double E_gam__GeV;
        double E_f__GeV;
        gsl_spline_object_2D gso_1D_so;
        gsl_spline_object_2D gso_2D_so;
    };


    int int_F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct test_IC fdata_in = *((struct test_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+1]) * exp(x[j*ndim+0]) * exp(x[j*ndim+0]) * 
                      gsl_so2D_eval( fdata_in.gso_2D_so, exp(x[j*ndim+0]), exp(x[j*ndim+1]) );
        }
        return 0;
    }


    struct test_IC fdata;


    double res_Gamma, res_spec_low, res_spec_high;
    double abserr;
    double xmin2D[2],xmax2D[2];

    fdata.gso_2D_so = gso_2D_so_IC_Gamma;
    xmin2D[0] = log(fdata.gso_2D_so.x_lim[0]);
    xmax2D[0] = log(fdata.gso_2D_so.x_lim[1]);
    xmin2D[1] = log(fdata.gso_2D_so.y_lim[0]);
    xmax2D[1] = log(fdata.gso_2D_so.y_lim[1]);
    hcubature_v( 1, int_F_IC, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_Gamma, &abserr );

    fdata.gso_2D_so = gso_2D_so_IC_spec_low;
    xmin2D[0] = log(fdata.gso_2D_so.x_lim[0]);
    xmax2D[0] = log(fdata.gso_2D_so.x_lim[1]);
    xmin2D[1] = log(fdata.gso_2D_so.y_lim[0]);
    xmax2D[1] = log(fdata.gso_2D_so.y_lim[1]);
    hcubature_v( 1, int_F_IC, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_spec_low, &abserr );

    fdata.gso_2D_so = gso_2D_so_IC_spec_high;
    xmin2D[0] = log(fdata.gso_2D_so.x_lim[0]);
    xmax2D[0] = log(fdata.gso_2D_so.x_lim[1]);
    xmin2D[1] = log(fdata.gso_2D_so.y_lim[0]);
    xmax2D[1] = log(fdata.gso_2D_so.y_lim[1]);
    hcubature_v( 1, int_F_IC, &fdata, 2, xmin2D, xmax2D, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_spec_high, &abserr );

    printf("Gamma: %e Spec_sum: %e (low: %e high: %e)\n", res_Gamma, res_spec_low + res_spec_high, res_spec_low, res_spec_high);

}







/*
 * The following function gives the scattered photon spectrum according to Blumenthal and Gould 1970 (Jones 1968) Eq. 2.48
 * E1 = E_gam/E_e
 * E2 = E_phot/m_e c^2
 * gamma_e = E_e/ m_e c^2
 *
 */

double d3NdtdEgamdEphot_sm1GeVm2( double E_e__GeV, double E_phot__GeV, double E_gam__GeV, double n_phot__Ephotm1cmm3 )
{
    double G( double q, double Gam_e )
    {
        return 2.*q*log(q) + (1.+2.*q)*(1.-q) + pow(Gam_e*q,2)*(1.-q)/(2.*(1.+Gam_e*q));
    }
    double Gam_e = 4. * E_phot__GeV/m_e__GeV * E_e__GeV/m_e__GeV;
    double q = E_gam__GeV/(Gam_e*(E_e__GeV - E_gam__GeV));

    if (q > pow(m_e__GeV/(2.*E_e__GeV),2) && q < 1. )
    {
//        return 3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * 1./(E_phot__GeV * pow(E_e__GeV/m_e__GeV, 2)) * n_phot__Ephotm1cmm3 * G( q, Gam_e );
        return n_phot__Ephotm1cmm3 * G( q, Gam_e );
    }
    else
    {
        return 0.;
    }
}


//Sigl p.545 Eq. 8.38
double d3NdtdEfdEphot_sm1GeVm2( double E_e__GeV, double E_phot__GeV, double E_f__GeV, double n_phot__Ephotm1cmm3 )
{
    double beta_e = sqrt(1. - 1./pow(E_e__GeV/m_e__GeV,2));
    double s = 2.*E_e__GeV*E_phot__GeV * (1.- beta_e*2./M_PI) + pow(m_e__GeV,2);
    double beta = (s-pow(m_e__GeV,2))/(s+pow(m_e__GeV,2));
    double EfEem1 = E_f__GeV/E_e__GeV;

    if ( EfEem1 >= (1.-beta)/(1.+beta) && EfEem1 <= 1. )
    {
        return 3./8. * sigma_T__mb * mb__cm2 * c__cmsm1 * n_phot__Ephotm1cmm3 * pow(m_e__GeV,2)/(s*E_e__GeV) * (1.+beta)/beta *
        ( EfEem1 + 1./EfEem1 + 2.*(1.-beta)/beta * (1.-1./EfEem1) + pow((1.-beta),2)/pow(beta,2) * pow((1.-1./EfEem1),2) );
    }
    else
    {
        return 0.;
    }
}







/*
 * Pass in function for photon spectrum to integrator as well as integration limits and limits for interpolation
 * Limits: 
 *         E_e__GeV >= m_e c^2
 * The output spline object is d2Ndtm1dEgamm1(E_gam__GeV,E_e__GeV)
 */

/*
gsl_spline_object_2D init_gso_2D_IC_new( double (*n_phot)(double *, double), double *n_phot_params, 
                                     double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                     size_t n_pts[2], char type[3] )
{
    int i,j;
//    double zdata[ n_pts[0] * n_pts[1] ];
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEfdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), 
                       (fdata_in.E_e__GeV - fdata_in.E_gam__GeV) + exp(x[j*ndim+0]), 
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) ) * (fdata_in.E_e__GeV-fdata_in.E_gam__GeV+exp(x[j*ndim+0]))/fdata_in.E_gam__GeV;
        }
        return 0;
    }

    //construct the arrays
    double E_gam__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
        io = logspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
    }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            zdata[ j * n_pts[0] + i ] = 0.;

            //set the integration limits for photon field

//            xmin[0] = fmax( log( E_gam__GeV_array[i]*pow(m_e__GeV,2)/(4.* E_e__GeV_array[j] * (E_e__GeV_array[j] - E_gam__GeV_array[i]))), log( E_phot__GeV_lims[0]) );
//            xmax[0] = fmin( log( E_gam__GeV_array[i]*E_e__GeV_array[j]/(E_e__GeV_array[j] - E_gam__GeV_array[i]) ), log( E_phot__GeV_lims[1] ) );

            xmin[0] = log( E_phot__GeV_lims[0]);
            xmax[0] = log( E_phot__GeV_lims[1] );


//            if ( (E_gam__GeV_array[i] < E_e__GeV_array[j]) && (xmin[0] < xmax[0])  )
            if ( (xmin[0] < xmax[0])  )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.E_gam__GeV = E_gam__GeV_array[i];
                fdata.n_phot = n_phot;
                fdata.n_phot_params = n_phot_params;

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
                    zdata[j * n_pts[0] + i] =  res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res); 
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 
        }
    }

    return gsl_so2D( n_pts[0], n_pts[1], E_gam__GeV_array, E_e__GeV_array, zdata );
}
*/

gsl_spline_object_2D init_gso_2D_IC( double (*n_phot)(double *, double), double *n_phot_params, 
                                     double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                     size_t n_pts[2], char type[3] )
{
    int i,j;
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.E_gam__GeV,
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.E_gam__GeV,
//                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );
        }
        return 0;
    }

    //construct the arrays
    double E_gam__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
        io = logspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
    }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;

    double gamma_e;

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            //set the integration limits for photon field
            xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( E_gam__GeV_array[i]*pow(m_e__GeV,2)/(4.* E_e__GeV_array[j] * (E_e__GeV_array[j] - E_gam__GeV_array[i]))) );
            xmax[0] = fmin( log( E_phot__GeV_lims[1] ), log( E_gam__GeV_array[i]*E_e__GeV_array[j]/(E_e__GeV_array[j] - E_gam__GeV_array[i]) ) );

//xmin[0] = log( E_phot__GeV_lims[0]);
//xmax[0] = log( E_phot__GeV_lims[1] );

            if ( (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.E_gam__GeV = E_gam__GeV_array[i];
                gamma_e = fdata.E_e__GeV/m_e__GeV;

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
//                    zdata[j * n_pts[0] + i] = res;
                    zdata[j * n_pts[0] + i] = 3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * 1./(pow(gamma_e, 2)) * res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res); 
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            }
        }
    }
    gsl_spline_object_2D gso2D_out = gsl_so2D( n_pts[0], n_pts[1], E_gam__GeV_array, E_e__GeV_array, zdata );
    free( zdata );
    return gso2D_out;
}

gsl_spline_object_2D init_gso_2D_IC_Gamma( double (*n_phot)(double *, double), double *n_phot_params, 
                                  double E_f__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                  size_t n_pts[2], char type[3] )
{
    int i,j;
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), (fdata_in.Delta_E__GeV + exp(x[j*ndim+0])),
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), 
//                       (fdata_in.Delta_E__GeV + exp(x[j*ndim+0])), fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEfdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), (fdata_in.E_e__GeV - fdata_in.Delta_E__GeV), 
//                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

        }
        return 0;
    }

    //construct the arrays
    double Delta_E__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
        io = logspace_array( n_pts[0], E_phot__GeV_lims[0], E_e__GeV_lims[1], Delta_E__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_f__GeV_lims[0], E_f__GeV_lims[1], Delta_E__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
    }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;
    double gamma_e;

    //set the integration limits for photon field
//            xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( (E_e__GeV_array[j]-E_f__GeV_array[i])*pow(m_e__GeV,2)/
//                      (4.*E_e__GeV_array[j]*E_f__GeV_array[i]-pow(m_e__GeV,2)) ) );
    xmin[0] = log( E_phot__GeV_lims[0]);
//            xmax[0] = fmin( log( 0.5 * ( ( E_f__GeV_array[i] - (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ) + 
//                      sqrt( pow(( E_f__GeV_array[i] + (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ),2) - 
//                      pow(m_e__GeV,2) ) ) ), log( E_phot__GeV_lims[1]) );
    xmax[0] = log( E_phot__GeV_lims[1] );

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            if ( (Delta_E__GeV_array[i] < E_e__GeV_array[j]) && (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.Delta_E__GeV = Delta_E__GeV_array[i];
                gamma_e = fdata.E_e__GeV/m_e__GeV;

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
//                    zdata[j * n_pts[0] + i] = res;
                    zdata[j * n_pts[0] + i] =  3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * 1./pow(gamma_e, 2) * res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res);
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 
        }
    }

    gsl_spline_object_2D gso2D_out = gsl_so2D( n_pts[0], n_pts[1], Delta_E__GeV_array, E_e__GeV_array, zdata );
    free( zdata );
    return gso2D_out;
}

data_object_2D init_do_2D_IC( double (*n_phot)(double *, double), double *n_phot_params, 
                                     double E_gam__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                     size_t n_pts[2] )
{
    int i,j;
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.E_gam__GeV,
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.E_gam__GeV,
//                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );
        }
        return 0;
    }

    //construct the arrays
    double E_gam__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;


    io = logspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
    if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
    if (io == 1){ printf("Error assigning log array in inverse Compton"); }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;

    double gamma_e;

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            //set the integration limits for photon field
            xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( E_gam__GeV_array[i]*pow(m_e__GeV,2)/(4.* E_e__GeV_array[j] * (E_e__GeV_array[j] - E_gam__GeV_array[i]))) );
            xmax[0] = fmin( log( E_phot__GeV_lims[1] ), log( E_gam__GeV_array[i]*E_e__GeV_array[j]/(E_e__GeV_array[j] - E_gam__GeV_array[i]) ) );

//xmin[0] = log( E_phot__GeV_lims[0]);
//xmax[0] = log( E_phot__GeV_lims[1] );

            if ( (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.E_gam__GeV = E_gam__GeV_array[i];
                gamma_e = fdata.E_e__GeV/m_e__GeV;
//printf("%i %i\n", i,j);
                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
//                    zdata[j * n_pts[0] + i] = res;
                    zdata[j * n_pts[0] + i] = 3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * 1./(pow(gamma_e, 2)) * res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res); 
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            }
        }
    }

    data_object_2D do_2D_IC = init_do2D( n_pts[0], n_pts[1], E_gam__GeV_array, E_e__GeV_array, zdata );

    free( zdata );

    return do_2D_IC;
}

data_object_2D init_do_2D_IC_Gamma( double (*n_phot)(double *, double), double *n_phot_params, 
                                  double E_f__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                  size_t n_pts[2] )
{
    int i,j;
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), (fdata_in.Delta_E__GeV + exp(x[j*ndim+0])),
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), 
//                       (fdata_in.Delta_E__GeV + exp(x[j*ndim+0])), fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

//            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEfdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), (fdata_in.E_e__GeV - fdata_in.Delta_E__GeV), 
//                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );

        }
        return 0;
    }

    //construct the arrays
    double Delta_E__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    io = logspace_array( n_pts[0], E_phot__GeV_lims[0], E_e__GeV_lims[1], Delta_E__GeV_array);
    if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
    if (io == 1){ printf("Error assigning log array in inverse Compton"); }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;
    double gamma_e;

    //set the integration limits for photon field
//            xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( (E_e__GeV_array[j]-E_f__GeV_array[i])*pow(m_e__GeV,2)/
//                      (4.*E_e__GeV_array[j]*E_f__GeV_array[i]-pow(m_e__GeV,2)) ) );
    xmin[0] = log( E_phot__GeV_lims[0]);
//            xmax[0] = fmin( log( 0.5 * ( ( E_f__GeV_array[i] - (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ) + 
//                      sqrt( pow(( E_f__GeV_array[i] + (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ),2) - 
//                      pow(m_e__GeV,2) ) ) ), log( E_phot__GeV_lims[1]) );
    xmax[0] = log( E_phot__GeV_lims[1] );

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            if ( (Delta_E__GeV_array[i] < E_e__GeV_array[j]) && (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.Delta_E__GeV = Delta_E__GeV_array[i];
                gamma_e = fdata.E_e__GeV/m_e__GeV;

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
//                    zdata[j * n_pts[0] + i] = res;
                    zdata[j * n_pts[0] + i] =  3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * 1./pow(gamma_e, 2) * res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res);
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 
        }
    }

    data_object_2D do_2D_IC_Gamma = init_do2D( n_pts[0], n_pts[1], Delta_E__GeV_array, E_e__GeV_array, zdata );

    free( zdata );

    return do_2D_IC_Gamma;
}

/*
gsl_spline_object_2D init_gso_2D_IC_Gamma_new( double (*n_phot)(double *, double), double *n_phot_params, 
                                  double E_f__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                  size_t n_pts[2], char type[3] )
{
    int i,j;
//    double zdata[ n_pts[0] * n_pts[1] ];
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
//            if (fdata_in.E_e__GeV + exp(x[j*ndim+0]) - fdata_in.E_f__GeV < 1.e-15){printf("Here\n");}
            fval[j] =  exp(x[j*ndim+0]) * d3NdtdEfdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.E_f__GeV, 
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );
        }
        return 0;
    }

    //construct the arrays
    double E_f__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
//        io = logspace_array( n_pts[0], E_f__GeV_lims[0], E_f__GeV_lims[1], E_f__GeV_array);
        io = logspace_array( n_pts[0], m_e__GeV, E_f__GeV_lims[1], E_f__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
//        io = logspace_array( n_pts[1], m_e__GeV, E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_f__GeV_lims[0], E_f__GeV_lims[1], E_f__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
    }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;


    //set the integration limits for photon field
//    xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( (E_e__GeV_array[j]-E_f__GeV_array[i])*pow(m_e__GeV,2)/
//                      (4.*E_e__GeV_array[j]*E_f__GeV_array[i]-pow(m_e__GeV,2)) ) );
    xmin[0] = log( E_phot__GeV_lims[0]);
//            xmax[0] = fmin( log( 0.5 * ( ( E_f__GeV_array[i] - (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ) + 
//                      sqrt( pow(( E_f__GeV_array[i] + (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ),2) - 
//                      pow(m_e__GeV,2) ) ) ), log( E_phot__GeV_lims[1]) );
    xmax[0] = log( E_phot__GeV_lims[1] );

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    for (i = 0; i < n_pts[0]; ++i)
    {

        for (j = 0; j < n_pts[1]; ++j)
        {
            zdata[ j * n_pts[0] + i ] = 0.;


//            if ( (xmin[0] < xmax[0]) )
            if ( (E_f__GeV_array[i] < E_e__GeV_array[j]) && (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.E_f__GeV = E_f__GeV_array[i];

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

                if (res >= 0.)
                {
                    zdata[j * n_pts[0] + i] =  res;
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res);
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 

//printf("%e ", zdata[j * n_pts[0] + i]);
        }
//printf("\n");
    }

    return gsl_so2D( n_pts[0], n_pts[1], E_f__GeV_array, E_e__GeV_array, zdata );

//    return gsl_so2D_temp( n_pts[0], n_pts[1], E_f__GeV_array, E_e__GeV_array, E_f__GeV_lims, E_e__GeV_lims, zdata );
}
*/

/*
gsl_spline_object_2D init_gso_2D_IC_Gamma_old( double (*n_phot)(double *, double), double *n_phot_params, 
                                  double E_f__GeV_lims[2], double E_e__GeV_lims[2], double E_phot__GeV_lims[2], 
                                  size_t n_pts[2], char type[3] )
{
    int i,j;
//    double zdata[ n_pts[0] * n_pts[1] ];
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    int F_IC( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct hcub_data_IC fdata_in = *((struct hcub_data_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            
            fval[j] =  d3NdtdEgamdEphot_sm1GeVm2( fdata_in.E_e__GeV, exp(x[j*ndim+0]),
                       ((fdata_in.E_e__GeV - fdata_in.E_f__GeV) + exp(x[j*ndim+0])),
                       fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) ) );// * ((fdata_in.E_e__GeV - fdata_in.E_f__GeV) + exp(x[j*ndim+0]))/fdata_in.E_f__GeV;// *
        }
        return 0;
    }

    //construct the arrays
    double E_f__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
//        io = logspace_array( n_pts[0], E_f__GeV_lims[0], E_f__GeV_lims[1], E_f__GeV_array);
        io = logspace_array( n_pts[0], m_e__GeV, E_f__GeV_lims[1], E_f__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
//        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        io = logspace_array( n_pts[1], m_e__GeV, E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in inverse Compton"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_f__GeV_lims[0], E_f__GeV_lims[1], E_f__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in inverse Compton"); }
    }

    struct hcub_data_IC fdata;
    double xmin[1], xmax[1];

    double res;
    double abserr;

    double gamma_e, E_phot;

    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            zdata[ j * n_pts[0] + i ] = 0.;
            gamma_e = E_e__GeV_array[j]/m_e__GeV;

            //set the integration limits for photon field
//            xmin[0] = fmax( log( E_phot__GeV_lims[0] ), log( (E_e__GeV_array[j]-E_f__GeV_array[i])*pow(m_e__GeV,2)/
//                      (4.*E_e__GeV_array[j]*E_f__GeV_array[i]-pow(m_e__GeV,2)) ) );
            xmin[0] = log( E_phot__GeV_lims[0]);
//            xmax[0] = fmin( log( 0.5 * ( ( E_f__GeV_array[i] - (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ) + 
//                      sqrt( pow(( E_f__GeV_array[i] + (pow(m_e__GeV,2)/(4.*E_e__GeV_array[j])) ),2) - 
//                      pow(m_e__GeV,2) ) ) ), log( E_phot__GeV_lims[1]) );
            xmax[0] = log( E_phot__GeV_lims[1] );

            if ( (E_f__GeV_array[i] < E_e__GeV_array[j]) && (xmin[0] < xmax[0]) )
            {
                fdata.E_e__GeV = E_e__GeV_array[j];
                fdata.E_f__GeV = E_f__GeV_array[i];
                fdata.n_phot = n_phot;
                fdata.n_phot_params = n_phot_params;

                //call integrator
                hcubature_v( 1, F_IC, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
//printf("%e\n", res);
                if (res >= 0.)
                {
                    zdata[j * n_pts[0] + i] =  3./4. * sigma_T__mb * mb__cm2 * c__cmsm1 * res/pow( gamma_e, 2 );
                }
                else 
                {
                    zdata[j * n_pts[0] + i] = 0.;
                    printf("%e - Integrator in IC routine returned value < 0 - The code monkey fucked up! \n", res);
                }
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 
        }
    }

//    return gsl_so2D( n_pts[0], n_pts[1], E_f__GeV_array, E_e__GeV_array, zdata );

    return gsl_so2D_temp( n_pts[0], n_pts[1], E_f__GeV_array, E_e__GeV_array, E_f__GeV_lims, E_e__GeV_lims, zdata );
}
*/

/* Gamma is d2N/(dE dt) */
double P_IC__GeVm1sm1( double E_e__GeV, double E_f__GeV, gsl_spline_object_2D gso2D_IC )
{
    if (E_f__GeV < E_e__GeV)
    {
        double Delta_E__GeV = E_e__GeV - E_f__GeV;
        return gsl_so2D_eval( gso2D_IC, Delta_E__GeV, E_e__GeV );
    }
    else
    {
        return 0.;
    }
}



/*
 * The following is an approximation taken from Fang+2020 to do this properly need to integrate the above
 */

double tau_IC__s( double E_e__GeV, double urad__GeVcmm3, double E_peak__GeV )
{
    double gamma_e = E_e__GeV/m_e__GeV;
    double beta = sqrt(1. - 1./pow(gamma_e,2));
//    double gam_K = 0.27*m_e__GeV/E_peak__GeV;
  //KN turned off
//  return (3. * pow(m_e__GeV,2))/( 4. * sigma_T_mb * mb_cm2 * c_cmsm1 * urad__GeVcmm3 * pow(beta,2) * E_e__GeV );
  //Klein Nishina correction Schlickeiser 2010
//  return 3. * m_e__GeV/( 4. * sigma_T__mb * mb__cm2 * c__cmsm1 * urad__GeVcmm3 * pow(beta,2) * (E_e__GeV/m_e__GeV) ) * (pow(E_e__GeV/m_e__GeV,2) + pow(gam_K,2))/pow(gam_K,2);

  //KN according to Fang et al. 2020
    double x = 4.*gamma_e*E_peak__GeV/m_e__GeV;
    double Y( double x )
    {
        if (x <= 1.5e-3)
        {
            return pow(M_PI,4)/15.;
        }
        else if (x < 150.)
        {
            return exp( -3.996e-2*pow(log(x),0) - 9.1e-1*pow(log(x),1) - 1.197e-1*pow(log(x),2) + 3.305e-3*pow(log(x),3) +
                   1.044e-3*pow(log(x),4) - 7.013e-5*pow(log(x),5) - 9.618e-6*pow(log(x),6));
        }
        else
        {
            return 3./4. * pow(M_PI/x,2) * (log(x) - 1.9805);
        }
    }

    return pow(M_PI,4) * m_e__GeV/( 20. * sigma_T__mb * mb__cm2 * c__cmsm1 * urad__GeVcmm3 * pow(beta,2) * gamma_e * Y(x) );
}


double tau_IC_fulltest__s( double E_e__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so )
{

    struct fdata_IC
    {
        double E_e__GeV;
        gsl_spline_object_2D gso_2D_so;
    };

    int F_IC_out( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct fdata_IC fdata_in = *((struct fdata_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * (fdata_in.E_e__GeV - exp(x[j*ndim+0])) *
                      P_IC__GeVm1sm1( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.gso_2D_so );
        }
        return 0;
    }

    struct fdata_IC fdata;
    double xmin[1], xmax[1];
    double res;
    double abserr;

    fdata.E_e__GeV = E_e__GeV;
    fdata.gso_2D_so = gso_2D_so;

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV);
    hcubature_v( 1, F_IC_out, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );

    return E_e__GeV/res;

}




double tau_IC_gamspec__s( double E_e__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so_high, gsl_spline_object_2D gso_2D_so_low )
{

    struct fdata_IC
    {
        double E_e__GeV;
        gsl_spline_object_2D gso_2D_so;
    };

    int F_IC_out( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct fdata_IC fdata_in = *((struct fdata_IC *)fdata);

        for (j = 0; j < npts; ++j)
        {

            fval[j] = exp(x[j*ndim+0]) * exp(x[j*ndim+0]) *
                      P_IC__GeVm1sm1( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.gso_2D_so );
        }
        return 0;
    }

    struct fdata_IC fdata;
    double xmin[1], xmax[1];
    double res_high, res_low;
    double abserr;

    fdata.E_e__GeV = E_e__GeV;

    fdata.gso_2D_so = gso_2D_so_high;
    xmin[0] = log(gso_2D_so_high.x_lim[0]);
    xmax[0] = log(gso_2D_so_high.x_lim[1]);
    hcubature_v( 1, F_IC_out, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_high, &abserr );

    fdata.gso_2D_so = gso_2D_so_low;
    xmin[0] = log(gso_2D_so_low.x_lim[0]);
    xmax[0] = log(gso_2D_so_low.x_lim[1]);
    hcubature_v( 1, F_IC_out, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_low, &abserr );

    return E_e__GeV/(res_high + res_low);

}

double dEdtm1_IC__GeVsm1( double E_e__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so_high, gsl_spline_object_2D gso_2D_so_low )
{
    return -1.* E_e__GeV/tau_IC_gamspec__s( E_e__GeV, E_e__GeV_lims, gso_2D_so_high, gso_2D_so_low );
}

int test_IC_integration( double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_ICGamma, gsl_spline_object_2D gso_2D_so_low, gsl_spline_object_2D gso_2D_so_high )
{

    struct F_data 
    {
        gsl_spline_object_2D gso_2D_so;
        gsl_spline_object_1D gso_1D_so;
        double E_e__GeV, E_gam__GeV;
    };
    
    int F_IC_Spec( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * gsl_spline_eval( fdata_in.gso_1D_so.spline, exp(x[j*ndim+0]), fdata_in.gso_1D_so.acc ) * exp(x[j*ndim+0]);
        }
        return 0;
    }

    int F_IC_ai( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * gsl_spline2d_eval( fdata_in.gso_2D_so.spline, fdata_in.E_gam__GeV, exp(x[j*ndim+0]), fdata_in.gso_2D_so.xacc, fdata_in.gso_2D_so.yacc );
        }
        return 0;
    }

    int F_IC_P( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * gsl_spline_eval( fdata_in.gso_1D_so.spline, exp(x[j*ndim+0]), fdata_in.gso_1D_so.acc );
        }
        return 0;
    }

    int F_IC_Gamma( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * P_IC__GeVm1sm1(fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.gso_2D_so) * (fdata_in.E_e__GeV - exp(x[j*ndim+0]));
        }
        return 0;
    }

    struct F_data fdata;
    double xmin[1], xmax[1];
    int i;

    double res, res_speclow, res_spechigh, res_Gamma_out;
    double abserr;

    int n_Esteps = 1001;


    double E__GeV_low[n_Esteps], dNdEm1_low[n_Esteps];
    logspace_array( n_Esteps, gso_2D_so_low.x_lim[0], gso_2D_so_low.x_lim[1], E__GeV_low);

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV_lims[1]);
    fdata.gso_2D_so = gso_2D_so_low;
    for (i = 0; i < n_Esteps; ++i)
    {
        res = 0.;
        fdata.E_gam__GeV = E__GeV_low[i];
        hcubature_v( 1, F_IC_ai, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        if ( res > 0. )
        {
            dNdEm1_low[i] = res;
        }
        else
        {
            dNdEm1_low[i] = 0.;
        }
    }

    gsl_spline_object_1D gso_1D_low = gsl_so1D( n_Esteps, E__GeV_low, dNdEm1_low );
    xmin[0] = log(gso_2D_so_low.x_lim[0]);
    xmax[0] = log(gso_2D_so_low.x_lim[1]);
    fdata.gso_1D_so = gso_1D_low;
    hcubature_v( 1, F_IC_Spec, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_speclow, &abserr );


    double E__GeV_high[n_Esteps], dNdEm1_high[n_Esteps];
    logspace_array( n_Esteps, gso_2D_so_high.x_lim[0], gso_2D_so_high.x_lim[1], E__GeV_high);

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV_lims[1]);
    fdata.gso_2D_so = gso_2D_so_high;
    for (i = 0; i < n_Esteps; ++i)
    {
        res = 0.;
        fdata.E_gam__GeV = E__GeV_high[i];
        hcubature_v( 1, F_IC_ai, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        if ( res > 0. )
        {
            dNdEm1_high[i] = res;
        }
        else
        {
            dNdEm1_high[i] = 0.;
        }
    }

    gsl_spline_object_1D gso_1D_high = gsl_so1D( n_Esteps, E__GeV_high, dNdEm1_high );
    xmin[0] = log(gso_2D_so_high.x_lim[0]);
    xmax[0] = log(gso_2D_so_high.x_lim[1]);
    fdata.gso_1D_so = gso_1D_high;
    hcubature_v( 1, F_IC_Spec, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_spechigh, &abserr );





    double E__GeV[n_Esteps], dNdEm1[n_Esteps];
    logspace_array( n_Esteps, E_e__GeV_lims[0], E_e__GeV_lims[1], E__GeV);

    xmin[0] = log(m_e__GeV);
    fdata.gso_2D_so = gso_2D_ICGamma;

    for (i = 0; i < n_Esteps; ++i)
    {
        res = 0.;
        xmax[0] = log(E__GeV[i]);
        fdata.E_e__GeV = E__GeV[i];
        hcubature_v( 1, F_IC_Gamma, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        if ( res > 0. )
        {
            dNdEm1[i] = res;
        }
        else
        {
            dNdEm1[i] = 0.;
        }

    }

    gsl_spline_object_1D gso_1D_P_out = gsl_so1D( n_Esteps, E__GeV, dNdEm1 );

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV_lims[1]);
    fdata.gso_1D_so = gso_1D_P_out;
    hcubature_v( 1, F_IC_P, &fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res_Gamma_out, &abserr );


    printf("Inverse Compton - total spec: %e total Gamma: %e ratio: %e\n", res_speclow+res_spechigh, res_Gamma_out, 
        (res_speclow+res_spechigh)/(res_Gamma_out));

    return 0.;

}


#endif
