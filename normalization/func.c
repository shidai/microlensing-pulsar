// functions for fit 
     
struct data 
{
    size_t n;
    double * y;
    //double * sigma;
};

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
     
struct data_params
{
    double l;
    double b;
};

double data (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pmass,xs,ys,zs;
    double theta,xb,yb,zb;

    struct data_params *fp = (struct data_params *)p;

    // source distribution: bulge
    //x0=fp->x0;
    x0=1.58;
    y0=0.62;
    //z0=fp->z0;  // kpc
    z0=0.43;  // kpc

    r0=8.5; // kpc

    l=fp->l;
    b=fp->b;
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=13.4*3.1415926/180.0;

    //pnum=fp->p0;
    pmass=0.8549; // M/kpc^3
  
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

    F=pmass*exp(-0.5*rds*rds);

    return F;
}

double lum (double *k, size_t dim, void *p)
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

double lum_dx0 (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    double theta,xb,yb,zb;
    double rds_dx0;

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

    rds_dx0=-pow(rds,-3.0)*(pow(xb/x0,2.0)+pow(yb/y0,2.0))*((xb*xb)/(x0*x0*x0));

    F=-pnum*exp(-0.5*rds*rds)*rds*rds_dx0;

    return F;
}

double lum_dy0 (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    double theta,xb,yb,zb;
    double rds_dy0;

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

    rds_dy0=-pow(rds,-3.0)*(pow(xb/x0,2.0)+pow(yb/y0,2.0))*((yb*yb)/(y0*y0*y0));

    F=-pnum*exp(-0.5*rds*rds)*rds*rds_dy0;

    return F;
}

double lum_dz0 (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    double theta,xb,yb,zb;
    double rds_dz0;

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

    rds_dz0=-pow(rds,-3.0)*(pow(zb,4.0)/pow(z0,5.0));

    F=-pnum*exp(-0.5*rds*rds)*rds*rds_dz0;

    return F;
}

double lum_dtheta (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    double theta,xb,yb,zb;
    double rds_dtheta,xb_dtheta,yb_dtheta;

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

    xb_dtheta=-xs*sin(theta)+ys*cos(theta);
    yb_dtheta=-xs*cos(theta)-ys*sin(theta);

    rds_dtheta=pow(rds,-3.0)*(pow(xb/x0,2.0)+pow(yb/y0,2.0))*(xb_dtheta*xb/(x0*x0)+yb_dtheta*yb/(y0*y0));

    F=-pnum*exp(-0.5*rds*rds)*rds*rds_dtheta;

    return F;
}

double lum_dp0 (double *k, size_t dim, void *p)
{
    double F;
    double r0,rds;
    double x0,y0,z0,b,l,xs,ys,zs;
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

    //pnum=fp->p0;
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

    F=exp(-0.5*rds*rds);

    return F;
}

int func_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    double res,err;
    double l,b,x0,y0,z0,theta,p0;

    size_t n=((struct data *)data)->n;
    double *y=((struct data *)data)->y;

    x0=gsl_vector_get (x,0);
    y0=gsl_vector_get (x,1);
    z0=gsl_vector_get (x,2);
    theta=gsl_vector_get (x,3);
    p0=gsl_vector_get (x,4);

    double xl[1]={3.0};
    double xu[1]={13.0};

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;

    const gsl_rng_type *T;
    gsl_rng *r;
    size_t calls=5000000;

    gsl_monte_function G;
    gsl_rng_env_setup ();
    T=gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);

    size_t i;
 
    for (i=0, l=-10.0*3.1415926/180.0; l<10.0*3.1415926/180.0; l+=(20.0*3.1415926/180.0)/10.0)
    {
        for (b=-10.0*3.1415926/180.0; b<10.0*3.1415926/180.0; b+=(20.0*3.1415926/180.0)/10.0)
        {
            struct params par = {l, b, x0, y0, z0, theta, p0};
            G.f=&lum;
            G.dim=1;
            G.params=&par;

            gsl_monte_vegas_integrate (&G, xl, xu, 1, 10000, r, s, &res, &err);

            do
            {
                gsl_monte_vegas_integrate (&G, xl, xu, 1, calls/5, r, s, &res, &err);
            }
            while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            gsl_vector_set (f, i, res-y[i]);
            i++;
        }
    }

    //printf ("%e %e %e\n", l, b, res);
    //fflush(stdout);        

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    return GSL_SUCCESS;
}

int func_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    double l,b,x0,y0,z0,theta,p0;

    size_t n=((struct data *)data)->n;

    x0=gsl_vector_get (x,0);
    y0=gsl_vector_get (x,1);
    z0=gsl_vector_get (x,2);
    theta=gsl_vector_get (x,3);
    p0=gsl_vector_get (x,4);

    double xl[1]={3.0};
    double xu[1]={13.0};

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;

    const gsl_rng_type *T;
    gsl_rng *r;
    size_t calls=5000000;

    gsl_rng_env_setup ();
    T=gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);

    gsl_monte_function G_dx0;
    gsl_monte_function G_dy0;
    gsl_monte_function G_dz0;
    gsl_monte_function G_dtheta;
    gsl_monte_function G_dp0;
    double res_dx0,err_dx0;
    double res_dy0,err_dy0;
    double res_dz0,err_dz0;
    double res_dtheta,err_dtheta;
    double res_dp0,err_dp0;

    size_t i;
 
    for (i=0, l=-10.0*3.1415926/180.0; l<10.0*3.1415926/180.0; l+=(20.0*3.1415926/180.0)/10.0)
    {
        for (b=-10.0*3.1415926/180.0; b<10.0*3.1415926/180.0; b+=(20.0*3.1415926/180.0)/10.0)
        {
            struct params par = {l, b, x0, y0, z0, theta, p0};

            {
                G_dx0.f=&lum_dx0;
                G_dx0.dim=1;
                G_dx0.params=&par;

                gsl_monte_vegas_integrate (&G_dx0, xl, xu, 1, 10000, r, s, &res_dx0, &err_dx0);

                do
                {
                    gsl_monte_vegas_integrate (&G_dx0, xl, xu, 1, calls/5, r, s, &res_dx0, &err_dx0);
                }
                while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            }
          
            /////////////////////////////////////////////////////////////////////////////////

            {
                G_dy0.f=&lum_dy0;
                G_dy0.dim=1;
                G_dy0.params=&par;

                gsl_monte_vegas_integrate (&G_dy0, xl, xu, 1, 10000, r, s, &res_dy0, &err_dy0);

                do
                {
                    gsl_monte_vegas_integrate (&G_dy0, xl, xu, 1, calls/5, r, s, &res_dy0, &err_dy0);
                }
                while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            }
          
            /////////////////////////////////////////////////////////////////////////////////

            {
                G_dz0.f=&lum_dz0;
                G_dz0.dim=1;
                G_dz0.params=&par;

                gsl_monte_vegas_integrate (&G_dz0, xl, xu, 1, 10000, r, s, &res_dz0, &err_dz0);

                do
                {
                    gsl_monte_vegas_integrate (&G_dz0, xl, xu, 1, calls/5, r, s, &res_dz0, &err_dz0);
                }
                while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            }
          
            /////////////////////////////////////////////////////////////////////////////////

            {
                G_dtheta.f=&lum_dtheta;
                G_dtheta.dim=1;
                G_dtheta.params=&par;

                gsl_monte_vegas_integrate (&G_dtheta, xl, xu, 1, 10000, r, s, &res_dtheta, &err_dtheta);

                do
                {
                    gsl_monte_vegas_integrate (&G_dtheta, xl, xu, 1, calls/5, r, s, &res_dtheta, &err_dtheta);
                }
                while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            }
          
            /////////////////////////////////////////////////////////////////////////////////

            {
                G_dp0.f=&lum_dp0;
                G_dp0.dim=1;
                G_dp0.params=&par;

                gsl_monte_vegas_integrate (&G_dp0, xl, xu, 1, 10000, r, s, &res_dp0, &err_dp0);

                do
                {
                    gsl_monte_vegas_integrate (&G_dp0, xl, xu, 1, calls/5, r, s, &res_dp0, &err_dp0);
                }
                while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            }

            gsl_matrix_set (J, i, 0, res_dx0);
            gsl_matrix_set (J, i, 1, res_dy0);
            gsl_matrix_set (J, i, 2, res_dz0);
            gsl_matrix_set (J, i, 3, res_dtheta);
            gsl_matrix_set (J, i, 4, res_dp0);
            i++;
        }
    }

    //printf ("%e %e %e\n", l, b, res);
    //fflush(stdout);        

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    return GSL_SUCCESS;
}

int func_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
    func_f (x, data, f);
    func_df (x, data, J);
     
    return GSL_SUCCESS;
}    
