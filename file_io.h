#ifndef file_io_h
#define file_io_h

//#include <unistd.h>
#include <stdio.h>
//#include <stdlib.h>
#include <math.h>

//#include "astro_const.h"
//#include "gal_rad.h"


char *string_cat( const char *str1, const char *str2 )
{
    char * string_out = malloc(strlen(str1) + strlen(str2) + 1);
    string_out = strcpy(string_out, str1);
    string_out = strcat(string_out, str2);
    return string_out;
}




void write_1D_file( unsigned int nx, double *xdata, char *label, char *filepath )
{
    FILE *outfile;
    unsigned int i;
     
    // open file for writing
    outfile = fopen( filepath, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filepath);
    }
    else
    {
        if (label != NULL)
        {
            fprintf( outfile, "%s\n", label );
        }

        for (i = 0; i < nx; i++)
        {
            fprintf( outfile, "%le ", xdata[i] );
        }
        fprintf( outfile, "\n" );
        fclose( outfile );
        printf("Successfully written file %s\n", filepath);
    }
    free( filepath );
    return;
}

void write_2D_file( unsigned int nx, unsigned int ny, double **data, char *label, char *filepath )
{
    FILE *outfile;
    unsigned int i, j;
     
    // open file for writing
    outfile = fopen( filepath, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filepath);
    }
    else
    {
        if (label != NULL)
        {
            fprintf( outfile, "%s\n", label );
        }

        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                fprintf( outfile, "%le ", data[i][j] );
            }
            fprintf( outfile, "\n" );
        }
        fclose( outfile );
        printf("Successfully written file %s\n", filepath);
    }
    free( filepath );
    return;
}


void write_2D_spec_file( unsigned int n_gal, unsigned int n_E, double **data, double *E__GeV, double **tau_gg, double **tau_EBL, double *distmod, char *label, char *filepath )
{
    FILE *outfile;
    unsigned int i, j;
     
    // open file for writing
    outfile = fopen( filepath, "w" );
    if (outfile == NULL)
    {
        printf("Error writing file %s: can't open output file\n", filepath);
    }
    else
    {
        if (label != NULL)
        {
            fprintf( outfile, "%s\n", label );
        }

        for (i = 0; i < n_gal; i++)
        {
            for (j = 0; j < n_E; j++)
            {
                fprintf( outfile, "%le ", data[i][j] * exp(-tau_gg[i][j]) * exp(-tau_EBL[i][j]) * distmod[i] * pow( E__GeV[j], 2 ) );
            }
            fprintf( outfile, "\n" );
        }
        fclose( outfile );
        printf("Successfully written file %s\n", filepath);
    }
    free( filepath );
    return;
}







#endif
