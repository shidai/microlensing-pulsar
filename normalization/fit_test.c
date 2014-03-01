#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
     
#include "func_test.c"

int integral (const gsl_vector * x);

int main (void)
{
    const size_t p=5;
    
    double x_init[5]={1.5,0.6,0.4,13.4*3.1415926/180.0,2.0*pow(10,9.0)};

    gsl_vector_view x=gsl_vector_view_array (x_init, p);
    integral (&x.vector);

    return 0;
}

   
     
