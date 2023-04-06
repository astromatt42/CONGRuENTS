#ifndef data_objects_h
#define data_objects_h

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "astro_const.h"
#include "gal_rad.h"
#include "gen_funcs.h"

typedef struct data_obj_1Ds
{
    size_t nx;
    double x_lim[2];
    double* x_data;
    double* z_data;
} data_object_1D;

void data_object_1D_free( data_object_1D do1D )
{
    free( do1D.x_data );
    free( do1D.z_data );
}

typedef struct data_obj_2Ds
{
    size_t nx, ny;
    double x_lim[2], y_lim[2];
    double* x_data;
    double* y_data;
    double* z_data;
} data_object_2D;

void data_object_2D_free( data_object_2D do2D )
{
    free( do2D.x_data );
    free( do2D.y_data );
    free( do2D.z_data );
}

typedef struct data_obj_3Ds
{
    size_t nx, ny, nf;
    double x_lim[2], y_lim[2], f_lim[2];
    double *x_data;
    double *y_data;
    double *f_data;
    double **z_data;
} data_object_3D;


void data_object_3D_free( data_object_3D do3D )
{
    free( do3D.x_data );
    free( do3D.y_data );
    free( do3D.f_data );
    free2D( do3D.nf, do3D.z_data );
}


typedef struct IC_objects
{
    //[0] is IC, [1] is Gamma
    data_object_2D do_2D_IC[2][4];
    data_object_3D do_3D_IC[2][2];
} IC_object;

void IC_object_free( IC_object ICo )
{
    unsigned short int i, j;
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 4; j++)
        {
            data_object_2D_free( ICo.do_2D_IC[i][j] );
        }
        for (j = 0; j < 2; j++)
        {
            data_object_3D_free( ICo.do_3D_IC[i][j] );
        }
    }
}


int write_do1D( data_object_1D do1D, char * filename )
{
    FILE * outfile;
    short int cfint;
//    unsigned short int fw;
    unsigned int i;
     
    // open file for writing
    outfile = fopen ( filename, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filename);
        cfint = 1;
    }
    else
    {
        fprintf( outfile, "%zu \n", do1D.nx );
        fprintf( outfile, "%le %le \n", do1D.x_lim[0], do1D.x_lim[1] );
        for (i = 0; i < do1D.nx; i++)
        {
            fprintf( outfile, "%le ", do1D.x_data[i] );
        }
        fprintf( outfile, "\n" );
        for (i = 0; i < do1D.nx; i++)
        {
            fprintf( outfile, "%le ", do1D.z_data[i] );
        }
        fprintf( outfile, "\n" );

        fclose( outfile );
    }

    return cfint;

}

int write_do2D( data_object_2D do2D, char * filename )
{
    FILE * outfile;
    short int cfint;
//    unsigned short int fw;
    unsigned int i;
     
    // open file for writing
    outfile = fopen ( filename, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filename);
        cfint = 1;
    }
    else
    {
        fprintf( outfile, "%zu \n", do2D.nx );
        fprintf( outfile, "%zu \n", do2D.ny );
        fprintf( outfile, "%le %le \n", do2D.x_lim[0], do2D.x_lim[1] );
        fprintf( outfile, "%le %le \n", do2D.y_lim[0], do2D.y_lim[1] );
        for (i = 0; i < do2D.nx; i++)
        {
            fprintf( outfile, "%le ", do2D.x_data[i] );
        }
        fprintf( outfile, "\n" );
        for (i = 0; i < do2D.ny; i++)
        {
            fprintf( outfile, "%le ", do2D.y_data[i] );
        }
        fprintf( outfile, "\n" );
        for (i = 0; i < do2D.nx * do2D.ny; i++)
        {
            fprintf( outfile, "%le ", do2D.z_data[i] );
        }
        fprintf( outfile, "\n" );

        fclose( outfile );
    }

    return cfint;

}

int write_do3D( data_object_3D do3D, char * filename )
{
    FILE * outfile;
    short int cfint;
//    unsigned short int fw;
    unsigned int i,j;
     
    // open file for writing
    outfile = fopen ( filename, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filename);
        cfint = 1;
    }
    else
    {
        fprintf( outfile, "%zu \n", do3D.nx );
        fprintf( outfile, "%zu \n", do3D.ny );
        fprintf( outfile, "%zu \n", do3D.nf );
        fprintf( outfile, "%le %le \n", do3D.x_lim[0], do3D.x_lim[1] );
        fprintf( outfile, "%le %le \n", do3D.y_lim[0], do3D.y_lim[1] );
        fprintf( outfile, "%le %le \n", do3D.f_lim[0], do3D.f_lim[1] );
        for (i = 0; i < do3D.nx; i++)
        {
            fprintf( outfile, "%le ", do3D.x_data[i] );
        }
        fprintf( outfile, "\n" );
        for (i = 0; i < do3D.ny; i++)
        {
            fprintf( outfile, "%le ", do3D.y_data[i] );
        }
        fprintf( outfile, "\n" );
        for (i = 0; i < do3D.nf; i++)
        {
            fprintf( outfile, "%le ", do3D.f_data[i] );
        }
        fprintf( outfile, "\n" );
        for (j = 0; j < do3D.nf; j++)
        {
            for (i = 0; i < do3D.nx * do3D.ny; i++)
            {
                fprintf( outfile, "%le ", do3D.z_data[j][i] );
            }
        }
        fprintf( outfile, "\n" );

        fclose( outfile );
    }
    return cfint;
}

gsl_spline_object_1D do1D_to_gso1D( data_object_1D do1D )
{
    return gsl_so1D( do1D.nx, do1D.x_data, do1D.z_data );
}

gsl_spline_object_2D do2D_to_gso2D( data_object_2D do2D )
{
    return gsl_so2D( do2D.nx, do2D.ny, do2D.x_data, do2D.y_data, do2D.z_data );
}

data_object_1D read_do1D( char * filename )
{
    data_object_1D do1D;
    FILE *infile;
//    short int cfint;
//    unsigned short int fr;
    unsigned int i;

    if ( access( filename, F_OK ) == 0 )
    {
        infile = fopen( filename, "r" );
        if (infile == NULL)
        {
            printf( "Error opening file %s: creating file\n", filename);
        }
        else
        {
            fscanf( infile, "%zu \n", &(do1D.nx) );
            fscanf( infile, "%le %le \n", &(do1D.x_lim[0]), &(do1D.x_lim[1]) );

            do1D.x_data = malloc( sizeof(double) * do1D.nx );
            do1D.z_data = malloc( sizeof(double) * do1D.nx );

            for (i = 0; i < do1D.nx; i++)
            {
                fscanf( infile, "%le ", &(do1D.x_data[i]) );
            }
            fscanf( infile, "\n" );
            for (i = 0; i < do1D.nx; i++)
            {
                fscanf( infile, "%le ", &(do1D.z_data[i]) );
            }
            fscanf( infile, "\n" );

            fclose (infile);
        }

    }
    else
    {
        printf( "File %s doesn't exist\n", filename);
    }

    return do1D;
}


data_object_2D read_do2D( char * filename )
{
    data_object_2D do2D;
    FILE *infile;
//    short int cfint;
//    unsigned short int fr;
    unsigned int i;

    if ( access( filename, F_OK ) == 0 )
    {
        infile = fopen( filename, "r" );
        if (infile == NULL)
        {
            printf( "Error opening file %s: creating file\n", filename);
        }
        else
        {
            fscanf( infile, "%zu \n", &(do2D.nx) );
            fscanf( infile, "%zu \n", &(do2D.ny) );
            fscanf( infile, "%le %le \n", &(do2D.x_lim[0]), &(do2D.x_lim[1]) );
            fscanf( infile, "%le %le \n", &(do2D.y_lim[0]), &(do2D.y_lim[1]) );

            do2D.x_data = malloc( sizeof(double) * do2D.nx );
            do2D.y_data = malloc( sizeof(double) * do2D.ny );
            do2D.z_data = malloc( sizeof(double) * do2D.nx * do2D.ny );

            for (i = 0; i < do2D.nx; i++)
            {
                fscanf( infile, "%le ", &(do2D.x_data[i]) );
            }
            fscanf( infile, "\n" );
            for (i = 0; i < do2D.ny; i++)
            {
                fscanf( infile, "%le ", &(do2D.y_data[i]) );
            }
            fscanf( infile, "\n" );
            for (i = 0; i < do2D.nx * do2D.ny; i++)
            {
                fscanf( infile, "%le ", &(do2D.z_data[i]) );
            }
            fscanf( infile, "\n" );

            fclose (infile);
        }

    }
    else
    {
        printf( "File %s doesn't exist\n", filename);
    }

    return do2D;
}


data_object_3D read_do3D( char * filename )
{
    data_object_3D do3D;
    FILE * infile;
//    short int cfint;
//    unsigned short int fw;
    unsigned int i,j;

    if ( access( filename, F_OK ) == 0 )
    {
        infile = fopen( filename, "r" );
        if (infile == NULL)
        {
            printf( "Error opening file %s: creating file\n", filename);
        }
        else
        {

            fscanf( infile, "%zu \n", &(do3D.nx) );
            fscanf( infile, "%zu \n", &(do3D.ny) );
            fscanf( infile, "%zu \n", &(do3D.nf) );
            fscanf( infile, "%le %le \n", &(do3D.x_lim[0]), &(do3D.x_lim[1]) );
            fscanf( infile, "%le %le \n", &(do3D.y_lim[0]), &(do3D.y_lim[1]) );
            fscanf( infile, "%le %le \n", &(do3D.f_lim[0]), &(do3D.f_lim[1]) );

            do3D.x_data = malloc(sizeof(double) * do3D.nx);
            do3D.y_data = malloc(sizeof(double) * do3D.ny);
            do3D.f_data = malloc(sizeof(double) * do3D.nf);
            do3D.z_data = malloc(sizeof *do3D.z_data * do3D.nf);
            if (do3D.z_data){for (i = 0; i < do3D.nf; i++){do3D.z_data[i] = malloc(sizeof *do3D.z_data[i] * do3D.nx * do3D.ny);}}

            for (i = 0; i < do3D.nx; i++)
            {
                fscanf( infile, "%le ", &(do3D.x_data[i]) );
            }
            fscanf( infile, "\n" );
            for (i = 0; i < do3D.ny; i++)
            {
                fscanf( infile, "%le ", &(do3D.y_data[i]) );
            }
            fscanf( infile, "\n" );
            for (i = 0; i < do3D.nf; i++)
            {
                fscanf( infile, "%le ", &(do3D.f_data[i]) );
            }
            fscanf( infile, "\n" );
            for (j = 0; j < do3D.nf; j++)
            {
                for (i = 0; i < do3D.nx * do3D.ny; i++)
                {
                    fscanf( infile, "%le ", &(do3D.z_data[j][i]) );
                }
            }
            fscanf( infile, "\n" );
            fclose( infile );
        }

    }
    else
    {
        printf( "File %s doesn't exist\n", filename);
    }
    return do3D;
}

int check_do1D( size_t nx, double x_lim[2], data_object_1D do1D )
{
    if ( do1D.nx == nx && abs(1. - do1D.x_lim[0]/x_lim[0]) < 1.e-6 && abs(1. - do1D.x_lim[1]/x_lim[1]) < 1.e-6 )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int check_do2D( size_t nx, size_t ny, double x_lim[2], double y_lim[2], data_object_2D do2D )
{
    if ( do2D.nx == nx && do2D.ny == ny && 
         abs(1. - do2D.x_lim[0]/x_lim[0]) < 1.e-6 && abs(1. - do2D.x_lim[1]/x_lim[1]) < 1.e-6 && 
         abs(1. - do2D.y_lim[0]/y_lim[0]) < 1.e-6 && abs(1. - do2D.y_lim[1]/y_lim[1]) < 1.e-6 )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int check_do3D( size_t nx, size_t ny, size_t nf, double x_lim[2], double y_lim[2], double f_lim[2], data_object_3D do3D )
{
//printf( "%zu %zu %zu %zu %zu %zu %le %le %le %le %le %le %le %le %le %le %le %le\n", do3D.nx, nx, do3D.ny, ny, do3D.nf, nf, do3D.x_lim[0], x_lim[0], do3D.x_lim[1], x_lim[1], do3D.y_lim[0], y_lim[0], do3D.y_lim[1], y_lim[1], do3D.f_lim[0], f_lim[0], do3D.f_lim[1], f_lim[1]);
    if ( do3D.nx == nx && do3D.ny == ny && do3D.nf >= nf &&
         abs(1. - do3D.x_lim[0]/x_lim[0]) < 1.e-6 && abs(1. - do3D.x_lim[1]/x_lim[1]) < 1.e-6 &&
         abs(1. - do3D.y_lim[0]/y_lim[0]) < 1.e-6 && abs(1. - do3D.y_lim[1]/y_lim[1]) < 1.e-6 &&
//         abs(1. - do3D.f_lim[0]/f_lim[0]) < 1.e-6 && abs(1. - do3D.f_lim[1]/f_lim[1]) < 1.e-6 )
//         do3D.f_lim[0] <= f_lim[0] && do3D.f_lim[1] >= f_lim[1] )
         do3D.f_lim[0] - f_lim[0] < 1.e-6 && f_lim[1] - do3D.f_lim[1] < 1.e6 )
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int check_do3D_FIELD( size_t nx, size_t ny, double x_lim[2], double y_lim[2], 
                int (*FIELD_num)( double, double, double, size_t * ), int (*FIELD_lim)( double, double, double * ),
                double T_FIELD_min__K, double T_FIELD_max__K, double Delta_T__K, data_object_3D do3D )
{
    size_t nf;
    FIELD_num( T_FIELD_min__K, T_FIELD_max__K, Delta_T__K, &nf );
    double f_lim[2];
    FIELD_lim( T_FIELD_min__K, T_FIELD_max__K, f_lim );

    return check_do3D( nx, ny, nf, x_lim, y_lim, f_lim, do3D );
}


data_object_2D init_do2D( size_t nx, size_t ny, double *x_data, double *y_data, double *z_data )
{
    data_object_2D do2D;
    unsigned int i,j;
    do2D.nx = nx;
    do2D.ny = ny; 
    do2D.x_lim[0] = x_data[0];
    do2D.x_lim[1] = x_data[nx-1];
    do2D.y_lim[0] = y_data[0];
    do2D.y_lim[1] = y_data[ny-1];
    do2D.x_data = malloc( sizeof(double) * nx );
    do2D.y_data = malloc( sizeof(double) * ny );
    do2D.z_data = malloc( sizeof(double) * nx * ny );

    for (i = 0; i < nx; i++)
    {
        do2D.x_data[i] = x_data[i];
    }
    for (i = 0; i < ny; i++)
    {
        do2D.y_data[i] = y_data[i];
    }
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            do2D.z_data[j * nx + i] = z_data[j * nx + i];
        }
    }

    return do2D;
}

data_object_1D init_do1D( size_t nx, double *x_data, double *z_data )
{
    data_object_1D do1D;
    unsigned int i;
    do1D.nx = nx;
    do1D.x_lim[0] = x_data[0];
    do1D.x_lim[1] = x_data[nx-1];
    do1D.x_data = malloc( sizeof(double) * nx );
    do1D.z_data = malloc( sizeof(double) * nx );

    for (i = 0; i < nx; i++)
    {
        do1D.x_data[i] = x_data[i];
    }
    for (i = 0; i < nx; i++)
    {
        do1D.z_data[i] = z_data[i];
    }

    return do1D;
}



















#endif
