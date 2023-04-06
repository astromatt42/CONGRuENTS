#ifndef synchrotron_h
#define synchrotron_h

#include <math.h>
#include <string.h>
#include <gsl_sf_bessel.h>
#include <cubature.h>
#include "physical_constants.h"
#include "gsl_decs.h"
#include "math_funcs.h"

#include "data_objects.h"

/* Ghisellini, Blumenthal or Schlickeiser
 * B__G in cgs so [B^2/8pi] = erg cm^{-3}
 */


//This populates the spline object for the function F(x)
gsl_spline_object_1D init_gso_1D_sync( double x_lims[2], size_t n_pts, char type[3] )
{

    int i;
    double zdata[n_pts];

    int F_sync( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * gsl_sf_bessel_Knu( 5./3., exp(x[j*ndim+0]));
        }
        return 0;
    }

    //construct the array
    double x_array[n_pts];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
        io = logspace_array( n_pts, x_lims[0], x_lims[1], x_array);
        if (io == 1){ printf("Error assigning log array in synchrotron"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts, x_lims[0], x_lims[1], x_array);
        if (io == 1){ printf("Error assigning lin array in synchrotron"); }
    }

    double xmin[1], xmax[1];
    double res;
    double abserr;

    double * fdata = NULL;

    for (i = 0; i < n_pts; ++i)
    {
        xmin[0] = log(x_array[i]);
        xmax[0] = fmin(log(x_lims[1]), log(150.));
        //call integrator
        if (xmin[0] < xmax[0])
        {
            hcubature_v( 1, F_sync, fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        }
        else
        {
            res = 0.;
        }

        if (res >= 0.)
        {
            zdata[i] =  x_array[i] * res;
        }
        else 
        {
            zdata[i] = 0.;
            printf("%e - Integrator in sync routine returned value < 0 - The code monkey fucked up! \n", res); 
        }
    }

    return gsl_so1D( n_pts, x_array, zdata );
}

data_object_1D init_do_1D_sync( double x_lims[2], size_t n_pts )
{

    int i;
    double zdata[n_pts];

    int F_sync( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        for (j = 0; j < npts; ++j)
        {
            fval[j] =  exp(x[j*ndim+0]) * gsl_sf_bessel_Knu( 5./3., exp(x[j*ndim+0]));
        }
        return 0;
    }

    //construct the array
    double x_array[n_pts];
    int io;

    io = logspace_array( n_pts, x_lims[0], x_lims[1], x_array);
    if (io == 1){ printf("Error assigning log array in synchrotron"); }


    double xmin[1], xmax[1];
    double res;
    double abserr;

    double * fdata = NULL;

    for (i = 0; i < n_pts; ++i)
    {
        xmin[0] = log(x_array[i]);
        xmax[0] = fmin(log(x_lims[1]), log(150.));
        //call integrator
        if (xmin[0] < xmax[0])
        {
            hcubature_v( 1, F_sync, fdata, 1, xmin, xmax, 100000, 0., 1e-8, ERROR_INDIVIDUAL, &res, &abserr );
        }
        else
        {
            res = 0.;
        }

        if (res >= 0.)
        {
            zdata[i] =  x_array[i] * res;
        }
        else 
        {
            zdata[i] = 0.;
            printf("%e - Integrator in sync routine returned value < 0 - The code monkey fucked up! \n", res); 
        }
    }

    data_object_1D do_1D_SY = init_do1D( n_pts, x_array, zdata );
    return do_1D_SY;
}



double dEdtm1_sync__GeVsm1( double E_e__GeV, double B__G )
{
    double gamma = E_e__GeV/m_e__GeV;
    double beta = sqrt(1. - 1./pow(gamma,2));
//    return -4./3. * sigma_T__mb * mb__cm2 * c__cmsm1 * pow(gamma,2) * pow(beta,2) * pow(B__G, 2)/(8.*M_PI) * erg__GeV;
    return -1. * sigma_T__mb * mb__cm2 * c__cmsm1 * pow(gamma,2) * pow(beta,2) * pow(B__G, 2)/(8.*M_PI) * erg__GeV;
}


double deldelEm1dEdtm1_sync__sm1( double E_e__GeV, double B__G )
{
    double gamma = E_e__GeV/m_e__GeV;
    double beta = sqrt(1. - 1./pow(gamma,2));
    return dEdtm1_sync__GeVsm1( E_e__GeV, B__G ) * 2./(gamma * m_e__GeV * pow(beta,2));
}


double tau_sync__s( double E_e__GeV, double B__G )
{
    return -1.* E_e__GeV/dEdtm1_sync__GeVsm1( E_e__GeV, B__G );
}

#endif
