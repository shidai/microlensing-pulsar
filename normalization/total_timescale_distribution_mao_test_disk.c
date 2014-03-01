// calculate event rate and time scale distribution for all objects, Mao's model, with monte calo
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <mpi.h>

double t;

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
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double rds,mfunc;
    //double mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double vl,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    double vl,vb,vrot,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,ve;
    double theta,xb,yb,zb;
    //double s,abar,albar,rbar,pmass,q,p_bar;

    // source distribution: bulge
    //x0=1.58*8.5/8.0;
    x0=1.58;
    //y0=0.2;
    y0=0.62;
    //y0=0.62;
    //z0=0.43*8.5/8.0;  // kpc
    z0=0.43;  // kpc
    //abar=1.0;  // kpc
    //albar=1.8;
    //q=0.35;
    //rbar=3.0; // kpc
    //pmass=1.22*pow(10.0,9.0);  // M*kpc^(-3)

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=13.4*3.1415926/180.0;

    //pnum=0.003869*pow(10,9);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*0.8549*pow(10.0,9.0);  // M*kpc^(-3)
    pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    //p_bar=pmass*pow(s/abar,-albar)*exp(-(s*s)/(rbar*rbar));

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

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

    // mass function of stars
    if (k[2]>0.7)
    {
        mfunc=pow(k[2]/0.7,-2.0);
    }
    else
    {
        mfunc=pow(k[2]/0.7,-1.3);
    }

    if (k[2]<=1.0)
    {
        M=k[2];
    }
    else if (k[2]>1.0 && k[2]<=8.0)
    {
        M=0.6;
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=1.35; 
    }
    else 
    {
        M=5.0;
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        //vrot=100.0*pow(xs*xs+ys*ys,0.5); // 100km/s,x/kpc
        vrot=100.0*xs; // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc
    }
    else
    {
        //vrot=100.0;
        vrot=100.0*xs/pow(xs*xs+ys*ys,0.5);
        //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
    }

    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vlb=k[4];
    //vbb=k[5];
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    ////////////////////////////////////////////////////////////////////
    vl=k[3];
    //vb;
    mean_vl=(220.0-vrot)*(k[1]/k[0]);
    mean_vb=0.0;
    sigma_vl=pow(30.0*30.0+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    sigma_vb=pow(20.0*20.0+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //vb=vbd-vbb*(k[1]/k[0]);

    //if ((184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl)<=0)
    //{
    //    ifunc_vbd=0.0;
    //}
    //else
    //{
    //    ifunc_vbd=1.0;
    //}

    // vbd as a function of vlb,vx,vy,vz,M
    //vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*k[3]*k[3]/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(k[3]-mean_vl)*(k[3]-mean_vl)/(sigma_vl*sigma_vl));
    vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    ve=365.0*13.57*Deff*pow(M,0.5)/t;

    //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);

    if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    {
        F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vlfunc*vbfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbfunc*vlfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-14)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc)*ifunc*ifunc_vb;
        //F=pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*(Deff*Deff*365.0/t)*(mfunc*pow(k[2]*M,0.5))*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(2.0*3.1415926,0.5)*110.0*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*p_disk*ifunc*ifunc_vbd*(Deff*Deff*365.0/t)*(mfunc)*M*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*pnum*exp(-0.5*rds*rds)*p_disk*ifunc*ifunc_vbd*(Deff*Deff/t)*(mfunc)*M*(vbdfunc*vldfunc*vxfunc*vyfunc*vzfunc);
    }
    else
    {
        F=0.0;
    }

    return F;
}

double source (double *k, size_t dim, void *params)
{
    double F;
    double r0,rds;
    //double r0;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double theta,xb,yb,zb;
    //double s,abar,albar,q,rbar,pmass,p_bar;
    // source distribution: bulge
    //x0=1.58*8.5/8.0;
    x0=1.58;
    //y0=0.2;
    y0=0.62;
    //y0=0.62;
    //z0=0.43*8.5/8.0;  // kpc
    z0=0.43;  // kpc
    //abar=1.0;  // kpc
    //albar=1.8;
    //q=0.35;
    //rbar=3.0; // kpc
    //pmass=1.22*pow(10.0,9.0);  // M*kpc^(-3)

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=13.4*3.1415926/180.0;

    //pnum=0.003869*pow(10,9);  // kpc^(-3)
    pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
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

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    //p_bar=pmass*pow(s/abar,-albar)*exp(-(s*s)/(rbar*rbar));
    F=k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds);
    //F=k[0]*k[0]*cos(b)*p_bar;

    return F;
}

int main (int argc, char *argv[])
{
    double total;
    int i;
    int id;  //  process rank
    int p;   //  number of processes

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    // all stars, disk to bulge
    // positive
    double res,err;

    //double xl[6]={3.0,0.0,0.03,-150.0,-400.0,-300.0};
    //double xu[6]={13.0,7.6,120.0,150.0,400.0,300.0};
    double xl[4]={3.0,0.0,0.03,-500.0};
    double xu[4]={13.0,8.0,120.0,500.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,4,0};

    size_t calls=500000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    // positive
    //double resb,errb;

    //double xlb[6]={3.0,0.0,0.03,0.0,0.0,0.0};
    //double xub[6]={13.0,7.6,120.0,600.0,600.0,600.0};
    //double xlb[4]={3.0,3.0,0.03,-300.0};
    //double xub[4]={8.0,8.0,120.0,300.0};

    //const gsl_rng_type *Tb;
    //gsl_rng *rb;

    //gsl_monte_function Gb={&g2,4,0};

    //size_t callsb=5000000;

    //gsl_rng_env_setup ();

    //Tb=gsl_rng_default;

    // source
    double res_s,err_s;

    double xl_s[1]={3.0};
    double xu_s[1]={13.0};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,1,0};

    size_t calls_s=500000000;

    T_s=gsl_rng_default;

    // calculation
    for (i=id;i<=200.0;i+=p)
    {
        t=pow(10.0,-0.5+0.018*i);    
    
        r = gsl_rng_alloc (T);

        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (4);
        gsl_monte_vegas_integrate (&G, xl, xu, 4, 10000, r, s,&res, &err);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 4, calls/5, r, s,&res, &err);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
        //display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        //rb = gsl_rng_alloc (Tb);

        //gsl_monte_vegas_state *sb = gsl_monte_vegas_alloc (4);
        //gsl_monte_vegas_integrate (&Gb, xlb, xub, 4, 10000, rb, sb, &resb, &errb);

        //do
        //{
        //    gsl_monte_vegas_integrate (&Gb, xlb, xub, 4, callsb/5, rb, sb, &resb, &errb);
        //}
        //while (fabs (gsl_monte_vegas_chisq (sb) - 1.0) > 0.5);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        r_s = gsl_rng_alloc (T_s);

        gsl_monte_vegas_state *s_s = gsl_monte_vegas_alloc (1);
        gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, 10000, r_s, s_s, &res_s, &err_s);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, calls_s/5, r_s, s_s, &res_s, &err_s);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s_s) - 1.0) > 0.5);
        //display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        //total = (2.0*fabs(resb))/res_s;
        total = (2.0*fabs(res))/res_s;
        //total = (fabs(res)+fabs(resb))/res_s;
        //total = (2.0*fabs(res)+2.0*fabs(resb))/res_s;
        printf ("%e %e\n", t,total);
        fflush(stdout);        

        gsl_monte_vegas_free (s);

        gsl_rng_free (r);

        //gsl_monte_vegas_free (sb);

        //gsl_rng_free (rb);

        gsl_monte_vegas_free (s_s);

        gsl_rng_free (r_s);

        //gsl_monte_vegas_free (s_neg);

        //gsl_rng_free (r_neg);

        //gsl_monte_vegas_free (sb_neg);

        //gsl_rng_free (rb_neg);

    }

    //printf ("Process %d is done\n", id);
    fflush (stdout);
    MPI_Finalize ();
    return 0;
}
