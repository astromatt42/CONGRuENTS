#ifndef bremsstrahlung_h
#define bremsstrahlung_h

#include <math.h>
#include "physical_constants.h"
#include "gsl_decs.h"
#include <gsl_spline2d.h>

#include "data_objects.h"







double dsigma_BSdEgamm1__mbGeVm1( double E_gam__GeV, double E_e__GeV, gsl_spline_object_1D Phi_1H_gso1D, gsl_spline_object_1D Phi_2H_gso1D )
{


    double phi_1, phi_2;
    double Delta = E_gam__GeV * m_e__GeV/(4.*alpha*E_e__GeV*(E_e__GeV - E_gam__GeV));
    if (Delta > 2)
    {
        double Z = 1.; //adjust this for a mix of H and He?
        double phi_u = 4.*( log( 2.*E_e__GeV/m_e__GeV * ((E_e__GeV - E_gam__GeV)/E_gam__GeV) ) - 0.5 );   
        phi_1 = phi_u * pow(Z,2);
        phi_2 = phi_u * pow(Z,2);
    }
    else
    {
        //for future reference, include He table as well and mix for a more realistic ISM and ionised species
        phi_1 = gsl_so1D_eval( Phi_1H_gso1D, Delta );
        phi_2 = gsl_so1D_eval( Phi_2H_gso1D, Delta );
    }

    return 3./(8.*M_PI) * sigma_T__mb * alpha * fmax(( (1. + pow(1. - E_gam__GeV/E_e__GeV, 2)) * phi_1 - 2./3. * (1. - E_gam__GeV/E_e__GeV) * phi_2 ),0.) * 1./E_gam__GeV;

}



gsl_spline_object_2D init_gso2D_BS( double E_gam__GeV_lims[2], double E_e__GeV_lims[2], size_t n_pts[2], char type[3] )
{
    double Delta_x[11] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10. };
    double Phi_1H[11] = { 45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65 };
    double Phi_2H[11] = { 44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41 };

    gsl_spline_object_1D Phi_1H_gso1D = gsl_so1D( 11, Delta_x, Phi_1H );
    gsl_spline_object_1D Phi_2H_gso1D = gsl_so1D( 11, Delta_x, Phi_2H );

    int i,j;
//    double zdata[ n_pts[0] * n_pts[1] ];
    double *zdata = malloc( sizeof zdata * n_pts[0] * n_pts[1] );

    //construct the arrays
    double E_gam__GeV_array[n_pts[0]];
    double E_e__GeV_array[n_pts[1]];
    int io;

    if ( strcmp( type, "log" ) == 0 )
    {
        io = logspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning log array in bremsstrahlung"); }
        io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning log array in bremsstrahlung"); }
    }
    else if ( strcmp( type, "lin" ) == 0 )
    {
        io = linspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], E_gam__GeV_array);
        if (io == 1){ printf("Error assigning lin array in bremsstrahlung"); }
        io = linspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], E_e__GeV_array);
        if (io == 1){ printf("Error assigning lin array in bremsstrahlung"); }
    }


    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            if ( E_gam__GeV_array[i] < E_e__GeV_array[j] )
            {
                zdata[ j * n_pts[0] + i ] = dsigma_BSdEgamm1__mbGeVm1( E_gam__GeV_array[i], E_e__GeV_array[j], Phi_1H_gso1D, Phi_2H_gso1D );
            }
            else
            {
                zdata[j * n_pts[0] + i] = 0.;
            } 
        }
    }

    gsl_spline_object_2D gso_2D_BS = gsl_so2D( n_pts[0], n_pts[1], E_gam__GeV_array, E_e__GeV_array, zdata );
    free( zdata );
    return gso_2D_BS;
}

data_object_2D init_do_2D_BS( double E_gam__GeV_lims[2], double E_e__GeV_lims[2], size_t n_pts[2] )
{
    double Delta_x[11] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10. };
    double Phi_1H[11] = { 45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65 };
    double Phi_2H[11] = { 44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41 };

    gsl_spline_object_1D Phi_1H_gso1D = gsl_so1D( 11, Delta_x, Phi_1H );
    gsl_spline_object_1D Phi_2H_gso1D = gsl_so1D( 11, Delta_x, Phi_2H );

    int i,j;

    data_object_2D do_2D_BS;

    do_2D_BS.nx = n_pts[0];
    do_2D_BS.ny = n_pts[1];
    for (i = 0; i < 2; ++i)
    {
        do_2D_BS.x_lim[i] = E_gam__GeV_lims[i];
        do_2D_BS.y_lim[i] = E_e__GeV_lims[i];
    }
    do_2D_BS.x_data = malloc( sizeof(double) * n_pts[0] );
    do_2D_BS.y_data = malloc( sizeof(double) * n_pts[1] );
    do_2D_BS.z_data = malloc( sizeof(double) * n_pts[0] * n_pts[1] );



    //construct the arrays
    int io;

    io = logspace_array( n_pts[0], E_gam__GeV_lims[0], E_gam__GeV_lims[1], do_2D_BS.x_data);
    if (io == 1){ printf("Error assigning log array in bremsstrahlung"); }
    io = logspace_array( n_pts[1], E_e__GeV_lims[0], E_e__GeV_lims[1], do_2D_BS.y_data);
    if (io == 1){ printf("Error assigning log array in bremsstrahlung"); }



    for (i = 0; i < n_pts[0]; ++i)
    {
        for (j = 0; j < n_pts[1]; ++j)
        {
            if ( do_2D_BS.x_data[i] < do_2D_BS.y_data[j] )
            {
                do_2D_BS.z_data[ j * n_pts[0] + i ] = dsigma_BSdEgamm1__mbGeVm1( do_2D_BS.x_data[i], do_2D_BS.y_data[j], Phi_1H_gso1D, Phi_2H_gso1D );
            }
            else
            {
                do_2D_BS.z_data[j * n_pts[0] + i] = 0.;
            } 
        }
    }

    return do_2D_BS;
}









/*
double sigma_BS__mb( double E_gam__GeV, double E_e__GeV )
{
    double phi_1, phi_2;
    double Delta = E_gam__GeV * m_e__GeV/(4.*alpha*E_e__GeV*(E_e__GeV - E_gam__GeV));
    if (Delta > 2)
    {
        double Z = 1.; //adjust this for a mix of H and He?
        double phi_u = 4.*( log( 2.*E_e__GeV/m_e__GeV * ((E_e__GeV - E_gam__GeV)/E_gam__GeV) ) - 0.5 );   
        phi_1 = phi_u * pow(Z,2);
        phi_2 = phi_u * pow(Z,2);
    }
    else
    {
        const double Delta_x[11] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10. };
        const double Phi_1_H[11] = { 45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65 };
        const double Phi_2_H[11] = { 44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41 };
        //for future reference, include He table as well and mix for a more realistic ISM and ionised species

        double Phi1( double Delta )
        {
            int i;
            for (i = 0; i < 10; i++)
            {
                if ( Delta > Delta_x[i] && Delta < Delta_x[i+1] )
                {
                    return (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_1_H[i+1] - Phi_1_H[i]) + Phi_1_H[i];
                }
            }
            return 0.;
        }

        double Phi2( double Delta )
        {
            int i;
            for (i = 0; i < 10; i++)
            {
                if ( Delta > Delta_x[i] && Delta < Delta_x[i+1] )
                {
                    return  (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_2_H[i+1] - Phi_2_H[i]) + Phi_2_H[i];
                }
            }
            return 0.;
        }
        phi_1 = Phi1( Delta );
        phi_2 = Phi2( Delta );
    }

    return 3./(8.*M_PI) * sigma_T__mb * alpha * ( (1. + pow(1. - E_gam__GeV/E_e__GeV, 2)) * phi_1 - 2./3. * (1. - E_gam__GeV/E_e__GeV) * phi_2 );

}
*/




double P_BS__GeVm1sm1( double E_e__GeV, double E_f__GeV, double n_H__cmm3, gsl_spline_object_2D gso2D_BS  )
{
    if (E_f__GeV < E_e__GeV)
    {
        double E_gam__GeV = E_e__GeV - E_f__GeV;
        return c__cmsm1 * n_H__cmm3 * gsl_so2D_eval( gso2D_BS, E_gam__GeV, E_e__GeV ) * mb__cm2;
    }
    else
    {
        return 0.;
    }
}


/*
 * The following is an approximation taken from Blumenthal(?) to do this properly need to integrate the above
 */

double tau_BS__s( double E_e__GeV, double n_H__cmm3 )
{
    double Delta_x[11] = { 0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10. };
    double Phi_1H[11] = { 45.79, 45.43, 45.09, 44.11, 42.64, 40.16, 34.97, 29.97, 24.73, 18.09, 13.65 };
//    double Phi_2H[11] = { 44.46, 44.38, 44.24, 43.65, 42.49, 40.19, 34.93, 29.78, 24.34, 17.28, 12.41 };

    //Schlickeiser 2002
    double fac_Delta( double Delta )
    {
        int i;
        for (i = 0; i < 10; ++i)
        {
            if ( Delta > Delta_x[i] && Delta < Delta_x[i+1] )
            {
                return (Delta - Delta_x[i])/(Delta_x[i+1] - Delta_x[i]) * (Phi_1H[i+1] - Phi_1H[i]) + Phi_1H[i];
            }
        }
        printf("Error in tauBS");
        return 1.;
    }

    double gamma = E_e__GeV/m_e__GeV;
    if (gamma < 15.)
    {
        return 2.*M_PI/(3. * alpha * c__cmsm1*sigma_T__mb * mb__cm2 * n_H__cmm3 * 2. * (log(gamma) + log(2.) - 1./3.));
    }
    else
    {
        double Delta = 1./( 4. * alpha * gamma );
        return 8.*M_PI/(3.9 * alpha * c__cmsm1*sigma_T__mb * mb__cm2 * n_H__cmm3 * fac_Delta( Delta ));
    }
}

double tau_BS_fulltest__s( double E_e__GeV, double E_e__GeV_lims[2], double n_H__cmm3, gsl_spline_object_2D gso2D_BS )
{

    struct fdata_BS
    {
        double E_e__GeV;
        double n_H__cmm3;
        gsl_spline_object_2D gso2D_BS;
    };

    int F_BS_out( unsigned ndim, size_t npts, const double *x, void *fdata, unsigned fdim, double *fval )
    {
        unsigned j;
        struct fdata_BS fdata_in = *((struct fdata_BS *)fdata);

        for (j = 0; j < npts; ++j)
        {
            fval[j] = exp(x[j*ndim+0]) * (fdata_in.E_e__GeV - exp(x[j*ndim+0])) * 
                      P_BS__GeVm1sm1( fdata_in.E_e__GeV, exp(x[j*ndim+0]), fdata_in.n_H__cmm3, fdata_in.gso2D_BS );
        }
        return 0;
    }

    struct fdata_BS fdata;
    double xmin[1], xmax[1];
    double res;
    double abserr;

    fdata.E_e__GeV = E_e__GeV;
    fdata.n_H__cmm3 = n_H__cmm3;
    fdata.gso2D_BS = gso2D_BS;

    xmin[0] = log(E_e__GeV_lims[0]);
    xmax[0] = log(E_e__GeV);
    hcubature_v( 1, F_BS_out, &fdata, 1, xmin, xmax, 100000, 0., 1e-6, ERROR_INDIVIDUAL, &res, &abserr );

    return E_e__GeV/res;

}


#endif

