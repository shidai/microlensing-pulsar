// func for test

struct params
{
    double l;
    double b;
    double x0;
    double y0;
    double z0;
    double theta;
    double p0;
};
     
double source (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    double theta,xb,yb,zb;

    struct params *fp = (struct params *)p;

    // source distribution: bulge
    x0=fp->x0;
    //x0=1.58;
    y0=fp->y0;
    z0=fp->z0;  // kpc
    //z0=0.43;  // kpc

    r0=8.0; // kpc

    l=fp->l;
    b=fp->b;
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=fp->theta;

    pnum=fp->p0;
    //pnum=pow(8.5/8.0,2.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    F=pnum*exp(-0.5*rds*rds);

    return F;
}

int integral (const gsl_vector *x)
{
    double res,err;
    double l,b,x0,y0,z0,theta,p0;

    x0=gsl_vector_get (x,0);
    y0=gsl_vector_get (x,1);
    z0=gsl_vector_get (x,2);
    theta=gsl_vector_get (x,3);
    p0=gsl_vector_get (x,4);

    double xl[1]={3.0};
    double xu[1]={13.0};

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G;
 
    struct params par = {l, b, x0, y0, z0, theta, p0};
    G.f=&source;
    G.dim=1;
    G.params=&par;

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G, xl, xu, 1, 10000, r, s, &res, &err);

    do
    {
        gsl_monte_vegas_integrate (&G, xl, xu, 1, calls/5, r, s, &res, &err);
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    printf ("%e %e %e\n", l, b, res);
    fflush(stdout);        

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    return 0;
}

