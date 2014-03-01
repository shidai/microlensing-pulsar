// bulge number 
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
    double F;
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,pmass,p_bulge;
    // Dd: k[1]
    //double M,mfunc;
    double theta;

    //r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=-2.68*3.1415926/180.0;
    //l=1.5*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=4.3656*9.0*pow(10.0,9.0)/3.14346;  // M*kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    //xs=-k[0]*cos(b)*cos(l)+r0;
    xs=k[0]*sin(k[2])*cos(k[1]);
    ys=k[0]*sin(k[2])*sin(k[1]);
    zs=k[0]*cos(k[2]);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    //F=k[0]*sin(k[2])*p_bulge;
    F=k[0]*k[0]*sin(k[2])*p_bulge;

    return F;
}

int main (void)
{
    double res,err;

    double xl[3]={0.0,0.0,0.0};
    double xu[3]={5.0,2.0*3.1415928,3.1415926};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,3,0};

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    // calculation
    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,&res, &err);

    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");
    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,&res, &err);
        //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    //display_results ("vegas final", res, err);

    printf ("Result: %e \n", res);

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    return 0;
}
