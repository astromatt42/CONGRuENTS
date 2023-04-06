#ifndef CR_funcs_h
#define CR_funcs_h

#include "ionisation.h"
#include "diffusion.h"
#include "synchrotron.h"
#include "inverse_Compton.h"
#include "bremsstrahlung.h"
#include "gsl_decs.h"

double P_total__GeVm1sm1( double E_e__GeV, double E_f__GeV, double n_gso2D, gsl_spline_object_2D * gso_2D_radfield, double n_H__cmm3, gsl_spline_object_2D gso2D_BS )
{
    double P_total = 0.;
    int i;

    //inverse Compton on all photon fields
    if ( E_e__GeV > E_f__GeV )
    {
        for (i = 0; i < n_gso2D; ++i)
        {
//if (isnan(P_IC__GeVm1sm1( E_e__GeV, E_f__GeV, gso_2D_radfield[i] ))){printf("Nan in IC\n");}
            P_total += P_IC__GeVm1sm1( E_e__GeV, E_f__GeV, gso_2D_radfield[i] );
        }
        //add Bremsstrahlung
//if (isnan(P_BS__GeVm1sm1( E_e__GeV, E_f__GeV, n_H__cmm3, gso2D_BS ))){printf("Nan in BS\n");}
        P_total += P_BS__GeVm1sm1( E_e__GeV, E_f__GeV, n_H__cmm3, gso2D_BS );
        return P_total;
    }
    else
    {
        return 0.;
    }
}

double dGammadlogEf_total__logGeVm1sm1( double E_e__GeV, double E_f__GeV, double n_gso2D, gsl_spline_object_2D * gso_2D_radfield, double n_H__cmm3, gsl_spline_object_2D gso2D_BS )
{
    double P_total = 0.;
    unsigned short int i;

    //inverse Compton on all photon fields
    if ( E_e__GeV > E_f__GeV )
    {
        for (i = 0; i < n_gso2D; ++i)
        {
            P_total += P_IC__GeVm1sm1( E_e__GeV, E_f__GeV, gso_2D_radfield[i] ) * E_f__GeV;
        }
        //add Bremsstrahlung
        P_total += P_BS__GeVm1sm1( E_e__GeV, E_f__GeV, n_H__cmm3, gso2D_BS ) * E_f__GeV;
        return P_total;
    }
    else
    {
        return 0.;
    }
}


double dEdtm1_total_disc__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return dEdtm1_sync__GeVsm1( E_e__GeV, B__G ) + dEdtm1_ion__GeVsm1( E_e__GeV, n_H__cmm3 );
}

double dEdtm1_total_halo__GeVsm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return dEdtm1_sync__GeVsm1( E_e__GeV, B__G ) + dEdtm1_plasma__GeVsm1( E_e__GeV, n_H__cmm3 );
}



double deldelEm1dEdtm1_total__sm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return deldelEm1dEdtm1_sync__sm1( E_e__GeV, B__G ) + deldelEm1dEdtm1_ion__sm1( E_e__GeV, n_H__cmm3 );
}

double deldelEm1dEdtm1_total_z2__sm1( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc )
{
    return deldelEm1dEdtm1_sync__sm1( E_e__GeV, B__G ) + deldelEm1dEdtm1_plasma__sm1( E_e__GeV, n_H__cmm3 );
}


double tau_total__s( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1, int n_radfields, 
                     double * u_rad__GeVcmm3, double * E_peak__GeV, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so )
{
    unsigned int i;
    double tauICm1 = 0.;
    for (i = 0; i < n_radfields; ++i)
    {
        tauICm1 += 1./tau_IC__s( E_e__GeV, u_rad__GeVcmm3[i], E_peak__GeV[i] );
    }

tauICm1 = 1./tau_IC_fulltest__s( E_e__GeV, E_e__GeV_lims, gso_2D_so ) ;

    double tau_tot__s = 1./(tauICm1 + 1./tau_BS__s( E_e__GeV, n_H__cmm3 ) + 1./tau_diff__s( E_e__GeV, h__pc, gso_1D_D__cm2sm1 ) + 
                        1./tau_sync__s( E_e__GeV, B__G ) + 1./tau_ion__s( E_e__GeV, n_H__cmm3 ));

    return tau_tot__s;
}

double tau_total_nodiff__s( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so_high, gsl_spline_object_2D gso_2D_so_low )
{
    return 1./(1./tau_IC_gamspec__s( E_e__GeV, E_e__GeV_lims, gso_2D_so_high, gso_2D_so_low ) + 1./tau_BS__s( E_e__GeV, n_H__cmm3 ) + 
                        1./tau_sync__s( E_e__GeV, B__G ) + 1./tau_ion__s( E_e__GeV, n_H__cmm3 ));
}

double tau_total_all__s( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so )
{
    return 1./(1./tau_IC_fulltest__s( E_e__GeV, E_e__GeV_lims, gso_2D_so ) + 1./tau_BS__s( E_e__GeV, n_H__cmm3 ) + 
                        1./tau_sync__s( E_e__GeV, B__G ) + 1./tau_ion__s( E_e__GeV, n_H__cmm3 ) + 
                        1./tau_diff__s( E_e__GeV, h__pc, gso_1D_D__cm2sm1 ));
}

double tau_total_all_z2__s( double E_e__GeV, double B__G, double n_H__cmm3, double h__pc, gsl_spline_object_1D gso_1D_D__cm2sm1, double E_e__GeV_lims[2], gsl_spline_object_2D gso_2D_so )
{
    return 1./(1./tau_IC_fulltest__s( E_e__GeV, E_e__GeV_lims, gso_2D_so ) + 1./tau_BS__s( E_e__GeV, n_H__cmm3 ) + 
                        1./tau_sync__s( E_e__GeV, B__G ) + 1./tau_plasma__s( E_e__GeV, n_H__cmm3 ) + 
                        1./tau_diff__s( E_e__GeV, h__pc, gso_1D_D__cm2sm1 ));
}

#endif











