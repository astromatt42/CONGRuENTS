#ifndef gen_funcs_h
#define gen_funcs_h

void free2D( int n, double **array )
{
    for (int i = 0; i < n; i++)
    {
        free(array[i]);
    }
    free(array);
}






#endif
