// calculate optical depth of pulsars by bulge at BW
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
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,rds,r0,pnum,p_bulge;
    // Dd: k[1]
    double rd,zd,p_pulsar;
    int ifunc;
    double Deff;
    double A,a,B,D,E;

    r0=8.0; // kpc

    x0=1.58;
    y0=0.62;
    z0=0.43;  // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;

    // source distribution: bulge

    //pnum=0.003869*pow(10,9);  // kpc^(-3)
    pnum=0.8549*1.388789*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=k[0]*cos(b)*cos(l)-r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density

    p_bulge=pnum*exp(-0.5*rds*rds);

    if (k[0]<=k[1])
    {
        ifunc=0.0;
    }
    else
    {
        ifunc=1.0;
    }

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    // pulsar distribution
    //A=41.0; // kpc^(-2)
    A=2000.0; // kpc^(-2)
    a=1.9;
    B=5.0;
    D=0.39; // kpc^(-1)
    E=0.33; // kpc

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);

    p_pulsar=A*pow(rd/r0,a)*exp(-B*((rd-r0)/r0))*D*exp(-fabs(zd)/E);  // pulsar density
 
    F=6.029*pow(10.0,-16.0)*1.4*k[0]*k[0]*cos(b)*p_bulge*p_pulsar*Deff*Deff*ifunc;

    return F;
}

double g2 (double *k, size_t dim, void *params)
{
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,rds,r0,pnum,p_bulge;
    // Dd: k[3]

    r0=8.0; // kpc

    x0=1.58;
    y0=0.62;
    z0=0.43;  // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;

    // source distribution: bulge 

    //pnum=0.003869*pow(10,9);  // kpc^(-3)
    pnum=0.8549*1.388789*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=k[0]*cos(b)*cos(l)-r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density

    p_bulge=pnum*exp(-0.5*rds*rds);

    F=k[0]*k[0]*cos(b)*p_bulge;

    return F;
}

int main (void)
{
    double res,err;

    double xl[2]={0.00001,0.00001};
    double xu[2]={20.0,20.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,2,0};

    size_t calls=50000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    double res2,err2;

    double xl2[1]={0.00001};
    double xu2[1]={20.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,1,0};

    size_t calls2=50000000;

    gsl_rng_env_setup ();

    T2=gsl_rng_default;

    // calculation
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,&res, &err);

    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,&res, &err);
        //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    //display_results ("vegas final", res, err);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    r2 = gsl_rng_alloc (T2);

    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G2, xl2, xu2, 1, 10000, r2, s2, &res2, &err2);

    do
    {
        gsl_monte_vegas_integrate (&G2, xl2, xu2, 1, calls2/5, r2, s2, &res2, &err2);
    }
    while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);

    printf ("Result: %e %e\n", res, res/res2);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    return 0;
}
