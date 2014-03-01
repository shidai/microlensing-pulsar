// normalize psr distribution 
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
    double rd,zd;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double F;

    r0=8.0;

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    B=4.01;
    E=0.05; // kpc

    rd=k[0]*sin(k[1]);
    zd=k[0]*cos(k[1]);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    F=2.0*3.1415926*k[0]*k[0]*sin(k[1])*p_pulsar;

    return F;
}

double g_r (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double p_pulsar;
    double A,a,B;
    double r0,r1;
    double F;

    r0=8.0;

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    A=1.37*2000.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    B=4.01;
    //E=0.05; // kpc

    rd=k[0]*sin(k[1]);
    //zd=k[0]*cos(k[1]);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)));  // pulsar radial density
    //p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    F=p_pulsar;
    //F=2.0*3.1415926*k[0]*k[0]*sin(k[1])*p_pulsar;

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

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // calculation
    
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
