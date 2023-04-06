# CONGRuENTS
COsmic-ray, Neutrino, Gamma-ray and Radio Non-Thermal Spectra

This code requires:

GNU Scientific Library - GSL 2.7.1 at time of writing

Cubature available from https://github.com/stevengj/cubature

Compile create_interp_objects.c to compute the interpolation objects separately - this is recommended if running multiple instances as it avoids having to recalculate them every time, however if making changes to the underlying code, ensure that they are deleted so they are recalculated on the next run.

For the main code compile spectra.c and link to the necessary libraries

