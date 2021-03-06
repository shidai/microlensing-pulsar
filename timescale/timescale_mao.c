// calculate percentage contribution to event rate by different objects
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

double g1 (double *k, size_t dim, void *params)
{
    double F,p_disk;
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double x0,y0,z0,b,l,pnum,xs,ys,zs;
    //double b,l,xs,ys,zs;
    double rds,mfunc;
    //double mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double v,vl,vb,vrot_disk,vld,vbd,vrot,vbb,vlb;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha,xb,yb,zb;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=0.000001;
    //l=-20.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    alpha=3.1415926-acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    alphas=3.1415926-acos( (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>0.0 ? (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)-1e-8 : (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)+1e-8);
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: bulge
    //x0=1.58;
    //y0=0.62;
    //z0=0.43;  // kpc
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    //theta=13.4*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;

    pnum=9.0*pow(10,9);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*0.8549*pow(10.0,9.0);  // M*kpc^(-3)
    //pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=-k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge

    //p_bar=pmass*pow(s/abar,-albar)*exp(-(s*s)/(rbar*rbar));

    // lense distribution: disk
    rdd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[1]*sin(b);

    //p0=0.0493*pow(10,9); // M*kpc^(-3)
    p0=1.388789*0.0493*pow(10.0,9.0); // kpc^(-3)
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

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

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
    //else if (k[2]>40.0) 
    {
        M=5.0;
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=-100.0*xs; // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=-100.0*xs/pow(xs*xs+ys*ys,0.5);
    }

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
          //vrot_disk=-220.0*xl/pow(xl*xl+yl*yl,0.5);
    //}

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=20.0;
    sigma_vdy=30.0;
    sigma_vdz=20.0;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*sqrt(M);   // km
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        //F=2.77*pow(10.0,-8.0)*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_disk*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
    }

    return F;
}

double g2 (double *k, size_t dim, void *params)
{
    double F;
    double x0,y0,z0,b,l,pnum,pnum2,xs,ys,zs;
    double xs2,ys2,zs2,r0;
    double rds,rds2,mfunc;
    //double ifunc,ifunc_vbd;
    double Deff;
    //double v,vl,vb,vld,vbd,vlb,vbb,vrot,vrot2;
    //double v,vrot,vrot2;
    // vld,vx,vy,vz are variables: k[3],k[4],k[5],k[6]
    //double vldfunc,vbdfunc,vyfunc,vzfunc;
    //double vl,vb,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M,re;
    double delta_func;
    double theta,alpha,xb,yb,zb,xb2,yb2,zb2;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot2;
    double d,R;

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    alpha=3.1415926-acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    alphas=3.1415926-acos( (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>0.0 ? (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)-1e-8 : (ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)+1e-8);
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    // source distribution: bulge
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc
    //x0=1.58;
    //y0=0.62;
    //z0=0.43;  // kpc

    //theta=13.4*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=(90.0-24.56)*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    pnum2=1.3888*9.0*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=pow(8.5/8.0,3.0)*0.8549*pow(10.0,9.0);  // M*kpc^(-3)
    //pnum2=pow(8.5/8.0,2.0)*1.3888*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum=1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=1.3888*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pnum2=0.8549*pow(10.0,9.0)/3.1;  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=-k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    // lense distribution: bulge
  
    //ys2=k[1]*cos(b)*cos(l)-r0;
    //xs2=-k[1]*cos(b)*sin(l);
    xs2=-k[1]*cos(b)*cos(l)+r0;
    ys2=-k[1]*cos(b)*sin(l);
    zs2=k[1]*sin(b);

    xb2=xs2*cos(theta)+ys2*sin(theta);
    //xb2=xs2*cos(theta)-ys2*sin(theta);
    yb2=-xs2*sin(theta)+ys2*cos(theta);
    //yb2=xs2*sin(theta)+ys2*cos(theta);
    zb2=zs2;

    rds2=sqrt(pow(xb2/x0,2.0)+pow(yb2/y0,2.0)+pow(zb2/z0,2.0));  // source density
    //rds2=pow(pow(pow(xb2/x0,2.0)+pow(yb2/y0,2.0),2.0)+pow(zb2/z0,4.0),0.25);  // lense density: bulge
    //rds2=pow(pow(pow(xs2/x0,2.0)+pow(ys2/y0,2.0),2.0)+pow(zs2/z0,4.0),0.25);  // lense density: bulge

    //if (k[0]<=k[1])
    //{
    //    ifunc=0.0;
    //}
    //else
    //{
    //    ifunc=1.0;
    //}

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

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
        //M=k[2];
    }
    else if (k[2]>8.0 && k[2]<=40.0)
    {
        M=1.35; 
        //M=k[2];
    }
    else 
    {
        M=5.0;
        //M=k[2];
    }    

    //A=0.23;
    //mfunc=fm*k[2]*exp(-A*Deff*Deff*k[2]/(t*t));

    // bulge rotational velocity function: source
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
    }

    // bulge rotational velocity function: lens
    if ((xs2*xs2+ys2*ys2)<1.0)
    {
        vrot2=100.0*sqrt(xs2*xs2+ys2*ys2); // 100km/s,x/kpc
        //vrot2=100.0*fabs(xs2); // 100km/s,x/kpc
    }
    else
    {
        vrot2=100.0;
        //vrot2=100.0*fabs(xs2)/pow(xs2*xs2+ys2*ys2,0.5);
    }

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: bulge velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[3];
    vdy=k[4];
    vdz=k[5];
    vdl=vrot2*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot2*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=110.0;
    sigma_vdy=82.5;
    sigma_vdz=66.3;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[6];
    vsy=k[7];
    vsz=k[8];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km
    // delta_func
    delta_func=(1.0/(sqrt(2.0*3.1415926)*0.01))*exp(-0.5*(log10((re/v)/(24.0*60.0*60.0))-t)*(log10((re/v)/(24.0*60.0*60.0))-t)/(0.01*0.01));

        //F=2.77*pow(10.0,-8.0)*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        //F=2.77*pow(10.0,-8.0)*(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*pnum2*exp(-rds2)*(Deff)*(sqrt(M)*mfunc)*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0)*delta_func/4.3656;
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

    // source distribution: bulge
    //x0=1.58;
    //y0=0.62;
    //z0=0.43*8.5/8.0;  // kpc
    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    r0=8.0; // kpc

    b=-2.75*3.1415926/180.0;
    l=1.16*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    //b=0.000001;
    //l=0.000001;
    //theta=13.4*3.1415926/180.0;
    theta=24.56*3.1415926/180.0;
    //theta=4.56*3.1415926/180.0;

    pnum=9.0*pow(10,9);  // kpc^(-3)
    //pnum=pow(8.5/8.0,3.0)*1.388789*0.8549*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    xs=-k[0]*cos(b)*cos(l)+r0;
    ys=-k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    xb=xs*cos(theta)+ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density
    //rds=pow(pow(pow(xb/x0,2.0)+pow(yb/y0,2.0),2.0)+pow(zb/z0,4.0),0.25);  // source density: bulge
    //rds=pow(pow(pow(xs/x0,2.0)+pow(ys/y0,2.0),2.0)+pow(zs/z0,4.0),0.25);  // source density: bulge

    //F=cos(b)*pnum*exp(-rds);
    //F=(1.0/(k[0]*k[0]))*cos(b)*pnum*exp(-rds);
    F=k[0]*k[0]*cos(b)*pnum*exp(-rds);
    //F=k[0]*k[0]*cos(b)*p_bar;

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

    double xl[9]={4.5,0.0,8.0,-100.0,-150.0,-100.0,-550.0,-450.0,-350.0};
    double xu[9]={11.5,11.5,40.0,100.0,150.0,100.0,550.0,450.0,350.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G={&g1,9,0};

    size_t calls=500000000;

    gsl_rng_env_setup ();

    T=gsl_rng_default;

    // all stars, bulge to bulge
    double res2,err2;

    //double xl2[5]={3.0,3.0,0.03,-800.0,-800.0};
    //double xu2[5]={13.0,13.0,120.0,800.0,800.0};
    double xl2[9]={4.5,4.5,8.0,-550.0,-450.0,-350.0,-550.0,-450.0,-350.0};
    double xu2[9]={11.5,11.5,40.0,550.0,450.0,350.0,550.0,450.0,350.0};

    const gsl_rng_type *T2;
    gsl_rng *r2;

    gsl_monte_function G2={&g2,9,0};

    size_t calls2=500000000;

    T2=gsl_rng_default;

    ////////////////////////////////////////////////////////////////////////////////////////////
    // source
    double res_s,err_s;

    double xl_s[1]={4.5};
    double xu_s[1]={11.5};

    const gsl_rng_type *T_s;
    gsl_rng *r_s;

    gsl_monte_function G_s={&source,1,0};

    size_t calls_s=500000000;

    T_s=gsl_rng_default;

    r_s = gsl_rng_alloc (T_s);

    gsl_monte_vegas_state *s_s = gsl_monte_vegas_alloc (1);
    gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, 10000, r_s, s_s, &res_s, &err_s);

    do
    {
        gsl_monte_vegas_integrate (&G_s, xl_s, xu_s, 1, calls_s/5, r_s, s_s, &res_s, &err_s);
    }
    while (fabs (gsl_monte_vegas_chisq (s_s) - 1.0) > 0.5);

    gsl_monte_vegas_free (s_s);

    gsl_rng_free (r_s);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

    r = gsl_rng_alloc (T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (9);

    r2 = gsl_rng_alloc (T2);

    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (9);

    // calculation
    for (i=id;i<=200.0;i+=p)
    {
        //t=pow(10.0,0.015*i);    
        t=-0.5+0.0175*i; 
    
        gsl_monte_vegas_integrate (&G, xl, xu, 9, 10000, r, s,&res, &err);

        do
        {
            gsl_monte_vegas_integrate (&G, xl, xu, 9, calls/5, r, s,&res, &err);
        }
        while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G2, xl2, xu2, 9, 10000, r2, s2, &res2, &err2);

        do
        {
            gsl_monte_vegas_integrate (&G2, xl2, xu2, 9, calls2/5, r2, s2, &res2, &err2);
        }
        while (fabs (gsl_monte_vegas_chisq (s2) - 1.0) > 0.5);

    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////

        ratio = (fabs(res)+fabs(res2))/(fabs(res_s));
        printf ("%e %e\n", t, ratio);
        fflush(stdout);        
    }

    gsl_monte_vegas_free (s);

    gsl_rng_free (r);

    gsl_monte_vegas_free (s2);

    gsl_rng_free (r2);

    fflush (stdout);
    MPI_Finalize ();
    return 0;
}
