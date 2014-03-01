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
     
#include "func.c"
     
#define N 3000
     
void print_state (size_t iter, gsl_multifit_fdfsolver * s);
     
int main (void)
{
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    unsigned int i, iter = 0;
    size_t n;
    const size_t p = 5;
     
    double y[N];

    /* This is the data to be fitted */
    double l,b;
    double xl[1]={3.0};
    double xu[1]={13.0};

    const gsl_rng_type *T1;
    gsl_rng *r1;
    size_t calls=5000000;

    gsl_rng_env_setup ();

    T1=gsl_rng_default;
    r1=gsl_rng_alloc (T1);

    gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc (1);

    gsl_monte_function G;

    double res,err;

    for (i=0, l=-10.0*3.1415926/180.0; l<10.0*3.1415926/180.0; l+=(20.0*3.1415926/180.0)/10.0)
    {
        for (b=-10.0*3.1415926/180.0; b<10.0*3.1415926/180.0; b+=(20.0*3.1415926/180.0)/10.0)
        {
            struct data_params data_params = {l, b};

            {
                G.f=&data;
                G.dim=1;
                G.params=&data_params;

                gsl_monte_vegas_integrate (&G, xl, xu, 1, 10000, r1, s1, &res, &err);

                do
                {
                    gsl_monte_vegas_integrate (&G, xl, xu, 1, calls/5, r1, s1, &res, &err);
                }
                while (fabs (gsl_monte_vegas_chisq (s1) - 1.0) > 0.5);
            }

            y[i]=res;
            printf ("%e\n", y[i]);
            i++;
        }
    }
    gsl_monte_vegas_free (s1);

    gsl_rng_free (r1);

    
    n=i;

    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    struct data d = {n, y};
    gsl_multifit_function_fdf f;
    double x_init[5] = { 1.5, 0.6, 0.4, 100, 20.0*3.1415926/180.0 };
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    const gsl_rng_type * type;
    gsl_rng * r;
     
    gsl_rng_env_setup();
     
    type = gsl_rng_default;
    r = gsl_rng_alloc (type);
     
    f.f = &func_f;
    f.df = &func_df;
    f.fdf = &func_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;
     
     
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);
     
    print_state (iter, s);
     
    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);
     
        printf ("status = %s\n", gsl_strerror (status));
     
        print_state (iter, s);
     
        if (status)
            break;
     
        status = gsl_multifit_test_delta (s->dx, s->x,
                                             1e-2, 1e-2);
    }
    while (status == GSL_CONTINUE && iter < 500);
     
    gsl_multifit_covar (s->J, 0.0, covar);
     
    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
     
    { 
        double chi = gsl_blas_dnrm2(s->f);
        double dof = n - p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
     
        printf("chisq = %g\n",  pow(chi, 2.0));
        printf("dof = %g\n",  dof);
        printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
     
        printf ("x0      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        printf ("y0 = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
        printf ("z0      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
        printf ("theta      = %.5f +/- %.5f\n", FIT(3), c*ERR(3));
        printf ("p0      = %.5f +/- %.5f\n", FIT(4), c*ERR(4));
    }
     
    printf ("status = %s\n", gsl_strerror (status));
    
    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_rng_free (r);
    return 0;
}
     
void print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f % 15.8f "
               "|f(x)| = %g\n",
               iter,
               gsl_vector_get (s->x, 0), 
               gsl_vector_get (s->x, 1),
               gsl_vector_get (s->x, 2), 
               gsl_vector_get (s->x, 3), 
               gsl_vector_get (s->x, 4), 
               gsl_blas_dnrm2 (s->f));
}

