// calculate event rate and time scale relation of pulsars over bulge star, with pulsar distribution, with monte calo
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
    double vl,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double vl,vb,vrot,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,ve;
    double theta,xb,yb,zb;
    //double s,abar,albar,rbar,pmass,q,p_bar;

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

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
        M=0.0; 
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
    vlb=k[4];
    vbb=k[5];
    vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    ////////////////////////////////////////////////////////////////////
    //vl=k[3];
    //vb;
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(30.0*30.0+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(20.0*20.0+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
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
    vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////
    // vb as a function of vl
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*k[3]*k[3]/900.0);
    vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(k[3]-mean_vl)*(k[3]-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    ve=365.0*13.57*Deff*pow(M,0.5)/t;

    //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);

    if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    {
        F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff*Deff*365.0/t)*(M*mfunc)*(vldfunc*vyfunc*vzfunc*vbdfunc)/4.3656;
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

double g2 (double *k, size_t dim, void *params)
{
    double F;
    double x0,y0,z0,b,l,pnum,xs,ys,zs,r0;
    double xs2,ys2,zs2,pnum2;
    double rds,rds2,mfunc;
    double Deff;
    double vl,vbd,vlb,vld,vbb,vrot,vrot2;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vxfunc;
    double vyfunc,vzfunc,vldfunc,vbdfunc;
    double M,ve;
    double theta,xb,yb,zb;
    double xb2,yb2,zb2;

    // source distribution: bulge
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // M*kpc^(-3)
    pnum2=1.3888*9.0*pow(10.0,9.0);  // M*kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density

    // lense distribution: bulge

    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb2=zs2;

    rds2=pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0),0.5);  // source density

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
        M=0.0; // use psr distribution instead of star distribution
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
        vrot=100.0*xs; // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0*xs/pow(xs*xs+ys*ys,0.5);
    }

    if ((xs2*xs2+ys2*ys2)<1.0)
    {
        vrot2=100.0*xs2; // 100km/s,x/kpc
    }
    else
    {
        vrot2=100.0*xs2/pow(xs2*xs2+ys2*ys2,0.5);
    }

    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    vlb=k[4];
    vbb=k[5];
    vld=vrot2+k[3];
    vl=vld-220.0+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);


    // vbd as a function of vlb,vx,vy,vz,M
    vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[3]*k[3]/(82.5*82.5));
    vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbd*vbd/(66.3*66.3));

    ve=365.0*13.57*Deff*pow(M,0.5)/t;

    //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);

    if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    {
        F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff*Deff*365.0/t)*(M*mfunc)*(vldfunc*vbdfunc*vyfunc*vzfunc)/4.3656;
        //F=2.3*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc)*ifunc*ifunc_vb;
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*M*(vbfunc*vlfunc);
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8275*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*(mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc)*ifunc*ifunc_vbd;
        //F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8275*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum2*exp(-0.5*rds2*rds2)*(Deff*Deff*365.0/t)*M*(mfunc)*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(2.0*3.1415926,0.5)*110.0*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum*exp(-0.5*rds2*rds2)*ifunc*(Deff*Deff*365.0/t)*(mfunc)*M*(vbdfunc*vldfunc*vyfunc*vzfunc);
        //F=pow(2.0*3.1415926,0.5)*110.0*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*1.0726*pow(10,-26)*k[0]*k[0]*cos(b)*pnum*exp(-0.5*rds*rds)*pnum*exp(-0.5*rds2*rds2)*ifunc*ifunc_vbd*(Deff*Deff*365.0/t)*(mfunc)*M*(vbdfunc*vldfunc*vyfunc*vzfunc);
    }
    else
    {
        F=0.0;
    }

    return F;
}

double g_psr (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double x0,y0,z0,rds,pnum,xs,ys,zs,p_bulge;
    double xb,yb,zb;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double b,l;
    double Deff;
    double F;
    double vl,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    double vldfunc,vbdfunc,vyfunc,vzfunc;
    double M,ve;
    double theta;

    M=1.4;

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.0;
    //l=0.0;
    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
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

    rds=pow(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0),0.5);  // source density
    p_bulge=pnum*exp(-rds);

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/5.6; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    B=4.01;
    E=0.05; // kpc

    rd=pow(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0,0.5);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=pow(k[1]*fabs(k[0]-k[1])/k[0],0.5);

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
    vlb=k[4];
    vbb=k[5];
    vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    ////////////////////////////////////////////////////////////////////

    // vbd as a function of vlb,vx,vy,vz,M
    vbd=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5)+vbb*(k[1]/k[0]);
    /////////////////////////////////////////////////////////////////////////////////////

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*290.0))*exp(-0.5*k[3]*k[3]/(290.0*290.0));
    //vldfunc=(1.0/2.0*180.0)*exp(-fabs(k[3])/180.0);
    vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*k[4]*k[4]/(82.5*82.5));
    vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*k[5]*k[5]/(66.3*66.3));
    //vbdfunc=(1.0/2.0*180.0)*exp(-fabs(vbd)/180.0);
    vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*290.0))*exp(-0.5*vbd*vbd/(290.0*290.0));
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(k[3]-mean_vl)*(k[3]-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    ve=365.0*13.57*Deff*pow(M,0.5)/t;

    //vfunc=pow(59.6*Deff*Deff*k[2]/(t*t)-k[3]*k[3],0.5)*exp(0.5*0.002377*k[3]*k[3]);

    if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( k[0] > k[1] && (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    //if ( (184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl) > 0 )
    {
        F=2.3*pow(fabs(ve*ve-vl*vl),-0.5)*ve*ve*3.8375*pow(10,-16)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff*Deff*365.0/t)*(M)*(vldfunc*vyfunc*vzfunc*vbdfunc);
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

int main (int argc, char *argv[])
{
    double ratio;
    int i;
    int id;  //  process rank
    int p;   //  number of processes

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);
    
    // all stars, disk to bulge
    // positive
    double res,err;

    //double xl[7]={0.0,0.0,0.03,-150.0,-550.0,-400.0,-350.0};
    //double xu[7]={50.0,50.0,120.0,150.0,550.0,400.0,350.0};
    double xl[6]={4.5,0.0,0.03,-150.0,-400.0,-350.0};
    double xu[6]={11.5,11.5,120.0,150.0,400.0,350.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g,6,0};

    size_t calls=5000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    double res2,err2;

    //double xl[7]={0.0,0.0,0.03,-150.0,-550.0,-400.0,-350.0};
    //double xu[7]={50.0,50.0,120.0,150.0,550.0,400.0,350.0};
    double xl2[6]={4.5,4.5,0.03,-400.0,-400.0,-350.0};
    double xu2[6]={11.5,11.5,120.0,400.0,400.0,350.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,6,0};

    size_t calls2=5000000;

    T2=gsl_rng_default;

    // negtive
    //double res_neg, err_neg;

    //double xl_neg[6]={0.0,0.0,0.03,-150.0,-400.0,-350.0};
    //double xu_neg[6]={50.0,50.0,120.0,150.0,400.0,350.0};

    //const gsl_rng_type *T_neg;
    //gsl_rng *r_neg;

    //gsl_monte_function G_neg={&g_neg,6,0};

    //size_t calls_neg=500000000;

    //T_neg=gsl_rng_default;

    // all stars, bulge to bulge
    //double res2_neg,err2_neg;

    //double xl[7]={0.0,0.0,0.03,-150.0,-550.0,-400.0,-350.0};
    //double xu[7]={50.0,50.0,120.0,150.0,550.0,400.0,350.0};
    //double xl2_neg[6]={0.0,0.0,0.03,-150.0,-400.0,-350.0};
    //double xu2_neg[6]={50.0,50.0,120.0,150.0,400.0,350.0};

    //const gsl_rng_type *T2_neg;
    //gsl_rng *r2_neg;

    //gsl_monte_function G2_neg={&g2_neg,6,0};

    //size_t calls2_neg=500000000;

    //T2_neg=gsl_rng_default;

    // neutron stars, pulsar to disk
    double res_x,err_x;

    double xl_x[5]={4.5,0.0,-1000.0,-400.0,-350.0};
    double xu_x[5]={11.5,11.5,1000.0,400.0,350.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,5,0};

    size_t calls_x=5000000;

    T_x=gsl_rng_default;
    
    // calculation
    for (i=id;i<=200.0;i+=p)
    {
        t=pow(10.0,0.015*i);    
    
        r = gsl_rng_alloc (T);

        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (6);
        gsl_monte_vegas_integrate (&G, xl, xu, 6, 10000, r, s,&res, &err);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 6, calls/5, r, s,&res, &err);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
        //display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        r2 = gsl_rng_alloc (T2);

        gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (6);
        gsl_monte_vegas_integrate (&G2, xl2, xu2, 6, 10000, r2, s2, &res2, &err2);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G2, xl2, xu2, 6, calls2/5, r2, s2, &res2, &err2);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);
        //display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        r_x=gsl_rng_alloc(T_x);
        
        gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (5);
        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 5, 10000, r_x, s_x,&res_x, &err_x);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        do
        {
            gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 5, calls_x/5, r_x, s_x,&res_x, &err_x);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        }
        while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);
        // display_results ("vegas final", res, err);

        /////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////

        //r_neg=gsl_rng_alloc(T_neg);

        //gsl_monte_vegas_state *s_neg = gsl_monte_vegas_alloc (6);
        //gsl_monte_vegas_integrate (&G_neg, xl_neg, xu_neg, 6, 10000, r_neg, s_neg, &res_neg, &err_neg);

        //do
        //{
        //    gsl_monte_vegas_integrate (&G_neg, xl_neg, xu_neg, 6, calls_neg/5, r_neg, s_neg, &res_neg, &err_neg);
        //}
        //while (fabs (gsl_monte_vegas_chisq (s_neg) - 1.0) > 0.5);

        //////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////

        //r2_neg = gsl_rng_alloc (T2_neg);

        //gsl_monte_vegas_state *s2_neg = gsl_monte_vegas_alloc (6);
        //gsl_monte_vegas_integrate (&G2_neg, xl2_neg, xu2_neg, 6, 10000, r2_neg, s2_neg, &res2_neg, &err2_neg);

        // display_results ("vegas warm-up", res, err);
        // printf ("converging...\n");
        //do
        //{
        //    gsl_monte_vegas_integrate (&G2_neg, xl2_neg, xu2_neg, 6, calls2_neg/5, r2_neg, s2_neg, &res2_neg, &err2_neg);
          //printf ("result = % .6ef sigma = % .6ef chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
        //}
        //while (fabs (gsl_monte_vegas_chisq (s2_neg) - 1.0) > 0.5);
        //display_results ("vegas final", res, err);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        //ratio = fabs(res2)/(fabs(res)+fabs(res_x));
        //ratio = fabs(res_x);
        ratio = fabs(res_x)/(fabs(res)+fabs(res_x)+fabs(res2));
        printf ("%e %e\n", t,ratio);
        fflush(stdout);        

        gsl_monte_vegas_free (s);

        gsl_rng_free (r);

        gsl_monte_vegas_free (s2);

        gsl_rng_free (r2);

        gsl_monte_vegas_free (s_x);

        //gsl_rng_free (r_neg);
    
        //gsl_monte_vegas_free (s_neg);

        //gsl_rng_free (r2_neg);
    
        //gsl_monte_vegas_free (s2_neg);

        gsl_rng_free (r_x);
    
    }

    //printf ("Process %d is done\n", id);
    fflush (stdout);
    MPI_Finalize ();
    return 0;
}
