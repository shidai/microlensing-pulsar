// normalize bulge distribution by optical depth 
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
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,r0,pmass,p_bulge;
    // Dd: k[1]
    double p_disk,rdd,zdd,p0,beta,ita,h1,h2,H;
    double Deff;
    //double M,mfunc;
    double theta;

    r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=-2.82*3.1415926/180.0;
    l=1.55*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    // lense distribution: disk
    rdd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zdd=k[1]*sin(b);

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

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);
 
    // mass function of stars
    //if (k[2]>0.7)
    //{
    //    mfunc=pow(k[2]/0.7,-2.0);
    //}
    //else
    //{
    //    mfunc=pow(k[2]/0.7,-1.3);
    //}

    //if (k[2]<=1.0)
    //{
    //    M=k[2];
    //}
    //else if (k[2]>1.0 && k[2]<=8.0)
    //{
    //    M=0.6;
    //}
    //else if (k[2]>8.0 && k[2]<=40.0)
    //{
    //    M=1.35; 
    //}
    //else 
    //{
    //    M=5.0;
    //}    

    if (k[0] > k[1])
    {
        F=6.029*pow(10.0,-16.0)*k[0]*cos(b)*p_bulge*p_disk*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_disk*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_disk*M*mfunc*Deff*Deff/4.3656;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double gb (double *k, size_t dim, void *params)
{
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,r0,pmass,p_bulge;
    // Dd: k[1]
    double xs2,ys2,zs2,xb2,yb2,zb2,rs2,pmass2,p_bulge2;
    double Deff;
    //double M,mfunc;
    double theta;

    r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=-2.82*3.1415926/180.0;
    l=1.55*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    // lenses distribution: bulge

    pmass2=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass2=1.388789*9.6*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //xb2=xs2*cos(theta)-ys2*sin(theta);
    //yb2=xs2*sin(theta)+ys2*cos(theta);
    zb2=zs2;

    rs2=pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0),0.5);  // source density

    p_bulge2=pmass2*exp(-rs2);

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);
 
    // mass function of stars
    //if (k[2]>0.7)
    //{
    //    mfunc=pow(k[2]/0.7,-2.0);
    //}
    //else
    //{
    //    mfunc=pow(k[2]/0.7,-1.3);
    //}

    //if (k[2]<=1.0)
    //{
    //    M=k[2];
    //}
    //else if (k[2]>1.0 && k[2]<=8.0)
    //{
    //    M=0.6;
    //}
    //else if (k[2]>8.0 && k[2]<=40.0)
    //{
    //    M=1.35; 
    //}
    //else 
    //{
    //    M=5.0;
    //}    

    if (k[0] > k[1])
    {
        F=6.029*pow(10.0,-16.0)*k[0]*cos(b)*p_bulge*p_bulge2*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_bulge2*Deff*Deff;
        //F=6.029*pow(10.0,-16.0)*k[0]*k[0]*cos(b)*p_bulge*p_bulge2*M*mfunc*Deff*Deff/4.3656;
    }
    else
    {
        F=0.0;
    }

    return F;
}

double g2 (double *k, size_t dim, void *params)
{
    // Ds, l, b: k[0]
    double F,l,b;
    double x0,y0,z0,xs,ys,zs,xb,yb,zb,rs,r0,pmass,p_bulge;
    double theta;

    r0=8.0; // kpc

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=0.89;
    //y0=0.3827;
    //z0=0.2492;  // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    b=-2.82*3.1415926/180.0;
    l=1.55*3.1415926/180.0;
    //theta=23.8*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;

    // source distribution: bulge

    pmass=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
    
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rs=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    p_bulge=pmass*exp(-rs);

    F=k[0]*cos(b)*p_bulge;
    //F=k[0]*k[0]*cos(b)*p_bulge;

    return F;
}

int main (void)
{
    double res,err;

    double xl[2]={4.5,0.0};
    double xu[2]={11.5,11.5};
    //double xl[3]={4.5,0.0,0.03};
    //double xu[3]={11.5,11.5,120.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,2,0};

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    double resb,errb;

    double xlb[2]={4.5,4.5};
    double xub[2]={11.5,11.5};
    //double xlb[3]={4.5,4.5,0.03};
    //double xub[3]={11.5,11.5,120.0};

    const gsl_rng_type *Tb;
    gsl_rng *rb;

    gsl_monte_function Gb={&gb,2,0};

    size_t callsb=5000000;

    Tb=gsl_rng_default;

    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////

    double res2,err2;

    double xl2[1]={4.5};
    double xu2[1]={11.5};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,1,0};

    size_t calls2=5000000;

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

    rb = gsl_rng_alloc (Tb);

    gsl_monte_vegas_state *sb = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&Gb, xlb, xub, 2, 10000, rb, sb, &resb, &errb);

    //display_results ("vegas warm-up", res, err);
    //printf ("converging...\n");
    do
    {
        gsl_monte_vegas_integrate (&Gb, xlb, xub, 2, callsb/5, rb, sb, &resb, &errb);
        //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    }
    while (fabs (gsl_monte_vegas_chisq (sb) - 1.0) > 0.5);
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

    printf ("Result: %e %e\n", res2, (fabs(res)+fabs(resb))/fabs(res2));

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (sb);

    gsl_rng_free (rb);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    return 0;
}
