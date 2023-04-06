#ifndef spec_integrate_h
#define spec_integrate_h

//IO for testing, don't need this after
//#include <stdio.h>
#include <math.h>
#include <gsl_spline.h>
#include <cubature.h>
#include "gsl_decs.h"

double spec_integrate_gso1D( gsl_spline_object_1D qess_so )
{
    struct F_data 
    {
        gsl_spline_object_1D spec_so;
    };

    int F( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  pow( exp(x[j*ndim+0]), 2 ) * gsl_spline_eval( fdata_in.spec_so.spline, exp(x[j*ndim+0]), fdata_in.spec_so.acc );
        }
        return 0;
    }

    double res, abserr;
    double xmin[1], xmax[1];
    struct F_data fdata;

    xmin[0] = log(qess_so.x_lim[0]);
    xmax[0] = log(qess_so.x_lim[1]);

    fdata.spec_so = qess_so;

    hcubature_v( 1, F, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    gsl_so1D_free( fdata.spec_so );

    return res;
}

double spec_integrate_gso1D_lim( gsl_spline_object_1D qess_so, double xlow, double xhigh )
{
    struct F_data 
    {
        gsl_spline_object_1D spec_so;
    };

    int F( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  pow( exp(x[j*ndim+0]), 2 ) * gsl_spline_eval( fdata_in.spec_so.spline, exp(x[j*ndim+0]), fdata_in.spec_so.acc );
        }
        return 0;
    }

    double res, abserr;
    double xmin[1], xmax[1];
    struct F_data fdata;

    xmin[0] = log(xlow);
    xmax[0] = log(xhigh);

    fdata.spec_so = qess_so;
    hcubature_v( 1, F, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );
    return res;
}

double spec_integrate( int n_E, double * E__GeV, double * dNdE__GeVm1 )
{
    struct F_data 
    {
        gsl_spline_object_1D spec_so;
    };

    int F( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data fdata_in = *((struct F_data *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  pow( exp(x[j*ndim+0]), 2 ) * gsl_spline_eval( fdata_in.spec_so.spline, exp(x[j*ndim+0]), fdata_in.spec_so.acc );
        }
        return 0;
    }

    double res, abserr;
    double xmin[1], xmax[1];
    struct F_data fdata;

    xmin[0] = log(E__GeV[0]);
    xmax[0] = log(E__GeV[n_E-1]);

    fdata.spec_so = gsl_so1D( n_E, E__GeV, dNdE__GeVm1 );

    hcubature_v( 1, F, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    gsl_so1D_free( fdata.spec_so );

    return res;
}

double ISRF_integrate__GeVcmm3( double (*n_phot)(double *, double), double *n_phot_params, double E_phot__GeV_lims[2] )
{
    struct F_data_ISRF
    {
        double *n_phot_params;
        double (*n_phot)( double *, double );
    };

    int F( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct F_data_ISRF fdata_in = *((struct F_data_ISRF *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] =  pow( exp(x[j*ndim+0]), 2 ) * fdata_in.n_phot( fdata_in.n_phot_params, exp(x[j*ndim+0]) );
        }
        return 0;
    }

    double res, abserr;
    double xmin[1], xmax[1];
    struct F_data_ISRF fdata;

    xmin[0] = log(E_phot__GeV_lims[0]);
    xmax[0] = log(E_phot__GeV_lims[1]);

    fdata.n_phot = n_phot;
    fdata.n_phot_params = n_phot_params;

    hcubature_v( 1, F, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    return res;
}

#endif
