// calculate the mass of disk stars
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

void display_results (char *title, double result, double error)
{
    printf ("%s ==================\n", title);
    printf ("result = % .6e\n", result);
    printf ("sigma = % .6e\n", error);
    //printf ("error = % .6f = %.2g sigma\n", result - exact,fabs (result - exact) / error);
}

double g (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,p0,beta,ita,h1,h2,H;
    double r0;

    r0=8.0;

    // lense distribution: disk
    rdd=k[0]*sin(k[1]);
    zdd=k[0]*cos(k[1]);

    p0=0.0493*pow(10,9); // M*kpc^(-3)
    //p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=1.3888*0.0493*pow(10.0,9.0); // kpc^(-3)
    //p0=0.0493*pow(10.0,9.0)/3.1; // kpc^(-3)
    beta=0.565;
    h1=0.270; // kpc
    h2=0.440; // kpc
    H=2.75; // kpc

    if (((rdd/9.025)+0.114)<=0.670)
    {
        ita=0.670;
    }
    else
    {
        ita=(rdd/9.025)+0.114;
    }

    p_disk=(p0/ita)*exp(-(rdd-r0)/H)*((1.0-beta)*pow(cosh(zdd/(ita*h1)),-2.0)+beta*exp(-fabs(zdd)/(ita*h2)));  // :disk

    F=2.0*3.1415926*k[0]*k[0]*sin(k[1])*p_disk;

    return F;
}

int main (int argc, char *argv[])
{
    double res,err;

    double xl[2]={0.0,0.0};
    double xu[2]={10.0,3.1415926};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,2,0};

    size_t calls=50000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,&res, &err);

    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,&res, &err);
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    printf ("%e\n", res);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    return 0;
}
