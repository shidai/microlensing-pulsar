// calculate the timescale of event rate due to NS 
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

double g_psr (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double x0,y0,z0,rds,pnum,xs,ys,zs;
    double xb,yb,zb;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double theta,alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;
    double l=k[8];
    double b=k[9];

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=3.1415926-acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=3.1415926-acos(-1);
    else alpha=3.1415926-acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>1.0)         alphas=3.1415926-acos(1);
    else if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)<-1.0)    alphas=3.1415926-acos(-1);
    else alphas=3.1415926-acos((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs));
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=-k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=0.58*2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/7.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.33; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=-100.0*xs; // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=-100.0*xs/pow(xs*xs+ys*ys,0.5);
        //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
    }

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
    //}

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[2];
    vdy=k[3];
    vdz=k[4];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=290.0;
    sigma_vdy=290.0;
    sigma_vdz=290.0;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[5];
    vsy=k[6];
    vsz=k[7];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

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
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
    }

    return F;
}

double g_psr_disk (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
    double p_disk,p_pulsar;
    double A,a,B,E;
    double r1;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //l=k[8];
    //b=k[9];
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    double l=k[8];
    double b=k[9];

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=3.1415926-acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=3.1415926-acos(-1);
    else alpha=3.1415926-acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>1.0)         alphas=3.1415926-acos(1);
    else if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)<-1.0)    alphas=3.1415926-acos(-1);
    else alphas=3.1415926-acos((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs));
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);

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

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/7.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.33; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // disk rotational velocity function
    vrot=220.0;

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
    //}

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[2];
    vdy=k[3];
    vdz=k[4];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=290.0;
    sigma_vdy=290.0;
    sigma_vdz=290.0;

    // source: disk velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[5];
    vsy=k[6];
    vsz=k[7];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=20.0;
    sigma_vsy=30.0;
    sigma_vsz=20.0;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

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
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));


    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

        F=2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
    }

    return F;
}

double g_psrr (double *k, size_t dim, void *params)
{   
    double rd,zd;
    double x0,y0,z0,rds,pnum,xs,ys,zs;
    double xb,yb,zb;
    double p_pulsar;
    double A,a,B,E;
    double r0,r1;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double theta,alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    //l=5.0*3.1415926/180.0;
    //b=5.0*3.1415926/180.0;
    double l=k[8];
    double b=k[9];

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=3.1415926-acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=3.1415926-acos(-1);
    else alpha=3.1415926-acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>1.0)         alphas=3.1415926-acos(1);
    else if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)<-1.0)    alphas=3.1415926-acos(-1);
    else alphas=3.1415926-acos((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs));
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: bulge

    x0=0.97473;
    y0=0.35107;
    z0=0.2644;  // kpc

    theta=24.56*3.1415926/180.0;

    pnum=9.0*pow(10.0,9.0);  // kpc^(-3)
    //pmass=0.8549*pow(10,9); // M/kpc^3
  
    xs=-k[0]*cos(b)*cos(l)+r0;
    //ys=k[0]*cos(b)*cos(l)-r0;
    //xs=-k[0]*cos(b)*sin(l);
    ys=-k[0]*cos(b)*sin(l);
    zs=k[0]*sin(b);

    //s=pow(xs*xs+(ys*ys+zs*zs)/(q*q),0.5);
    xb=xs*cos(theta)+ys*sin(theta);
    //xb=xs*cos(theta)-ys*sin(theta);
    yb=-xs*sin(theta)+ys*cos(theta);
    //yb=xs*sin(theta)+ys*cos(theta);
    zb=zs;

    rds=sqrt(pow(xb/x0,2.0)+pow(yb/y0,2.0)+pow(zb/z0,2.0));  // source density

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=0.58*2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/7.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.33; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // bulge rotational velocity function
    if ((xs*xs+ys*ys)<1.0)
    {
        vrot=100.0*sqrt(xs*xs+ys*ys); // 100km/s,x/kpc
        //vrot=-100.0*xs; // 100km/s,x/kpc
        //vrot=100.0*fabs(xs); // 100km/s,x/kpc
    }
    else
    {
        vrot=100.0;
        //vrot=-100.0*xs/pow(xs*xs+ys*ys,0.5);
        //vrot=100.0*fabs(xs)/pow(xs*xs+ys*ys,0.5);
    }

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
    //}

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[2];
    vdy=k[3];
    vdz=k[4];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=290.0;
    sigma_vdy=290.0;
    sigma_vdz=290.0;

    // source: bulge velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[5];
    vsy=k[6];
    vsz=k[7];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=110.0;
    sigma_vsy=82.5;
    sigma_vsz=66.3;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

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
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));

    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    double re;
    re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km

        F=(re/v)*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*pnum*exp(-rds)*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
    }

    return F;
}

double g_psr_diskr (double *k, size_t dim, void *params)
{   
    double rdd,zdd,r0,p0,beta,ita,h1,h2,H;
    double rd,zd;
    double p_disk,p_pulsar;
    double A,a,B,E;
    double r1;
    double Deff;
    double F;
    //double v,vl,vb,vrot,vrot_disk,sigma_y,sigma_z,mean_vl,mean_vb,sigma_vl,sigma_vb,vlfunc,vbfunc;
    double M;
    double alpha;
    double v,vl,vb,vol,vob,vdl,vdb,vsl,vsb,vdx,vdy,vdz,vsx,vsy,vsz,sigma_vdx,sigma_vdy,sigma_vdz,sigma_vsx,sigma_vsy,sigma_vsz,vdxfunc,vdyfunc,vdzfunc,vsxfunc,vsyfunc,vszfunc,vrot,vrot_disk;
    double d,R;

    r0=8.0; // kpc

    //l=k[8];
    //b=k[9];
    //b=-2.75*3.1415926/180.0;
    //l=1.16*3.1415926/180.0;
    double l=k[8];
    double b=k[9];

    if ( k[0] <= k[1] )
    {
        F=0.0;
    }
    else
    {
    d=k[1]*cos(b);
    R=sqrt(r0*r0+d*d-2.0*r0*d*cos(l));
    if ((d*d+R*R-r0*r0)/(2.0*d*R)>1.0)         alpha=3.1415926-acos(1);
    else if ((d*d+R*R-r0*r0)/(2.0*d*R)<-1.0)    alpha=3.1415926-acos(-1);
    else alpha=3.1415926-acos((d*d+R*R-r0*r0)/(2.0*d*R));
    if(sin(l)<0) alpha=-1.0*alpha;

    double ds,Rs,alphas;
    ds=k[0]*cos(b);
    Rs=sqrt(r0*r0+ds*ds-2.0*r0*ds*cos(l));
    if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)>1.0)         alphas=3.1415926-acos(1);
    else if ((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs)<-1.0)    alphas=3.1415926-acos(-1);
    else alphas=3.1415926-acos((ds*ds+Rs*Rs-r0*r0)/(2.0*ds*Rs));
    if(sin(l)<0) alphas=-1.0*alphas;
    //alpha=acos( (d*d+R*R-r0*r0)/(2.0*d*R)>0.0 ? (d*d+R*R-r0*r0)/(2.0*d*R)-1e-8 : (d*d+R*R-r0*r0)/(2.0*d*R)+1e-8);
    //alpha=asin(r0*sin(l)/R);
    //alpha=2.0*3.1415926-asin(r0*sin(l)/R);

    M=1.4;

    // source distribution: disk
    rdd=sqrt(k[0]*k[0]*cos(b)*cos(b)-2.0*k[0]*r0*cos(b)*cos(l)+r0*r0);
    zdd=k[0]*sin(b);

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

    // lens distribution: psr

    //A=41.0*1.6*pow(10.0,6.0); // kpc^(-2)
    //A=41.0*8.61*pow(10.0,6.0); // kpc^(-2)  by Lorimer 2006
    //A=2000.0*4.5*pow(10.0,5.0); // kpc^(-2)    by Kaspi 2006
    //A=2000.0*10000/1.18; // kpc^(-2)    by Kaspi 2006
    A=2000.0*10000/7.0; // kpc^(-2)    by Kaspi 2006
    r1=0.55;  // kpc
    a=1.64;
    //B=9.01;
    B=4.01;
    E=0.33; // kpc

    rd=sqrt(k[1]*k[1]*cos(b)*cos(b)-2.0*k[1]*r0*cos(b)*cos(l)+r0*r0);
    zd=k[1]*sin(b);

    p_pulsar=A*pow((rd+r1)/(r0+r1),a)*exp(-B*((rd-r0)/(r0+r1)))*exp(-fabs(zd)/E);  // pulsar density

    Deff=sqrt(k[1]*fabs(k[0]-k[1])/k[0]);

    // disk rotational velocity function
    vrot=220.0;

    // disk rotational velocity function
    //if ((xl*xl+yl*yl)<4.0)
    //{
    //    vrot_disk=220.0*xl; // 100km/s,x/kpc
    //}
    //else
    //{
          vrot_disk=220.0;
    //}

    // observor velocity
    vol=220.0*cos(l);
    vob=220.0*sin(l)*sin(b);

    // deflector: disk velocity
    // vdx,vdy,vdz,sigma_dx=20.0,sigma_dy=30.0,sigma_dz=20
    vdx=k[2];
    vdy=k[3];
    vdz=k[4];
    vdl=vrot_disk*cos(alpha)+vdx*sin(l)-vdy*cos(l);
    vdb=-vrot_disk*sin(alpha)*sin(b)+vdx*cos(l)*sin(b)+vdy*sin(l)*sin(b)+vdz*cos(b);
    sigma_vdx=290.0;
    sigma_vdy=290.0;
    sigma_vdz=290.0;

    // source: disk velocity
    // vsx,vsy,vsz,sigma_sx=110.0,sigma_sy=82.5,sigma_sz=66.3
    vsx=k[5];
    vsy=k[6];
    vsz=k[7];
    vsl=vrot*cos(alphas)+vsx*sin(l)-vsy*cos(l);
    vsb=-vrot*sin(alphas)*sin(b)+vsx*cos(l)*sin(b)+vsy*sin(l)*sin(b)+vsz*cos(b);
    sigma_vsx=20.0;
    sigma_vsy=30.0;
    sigma_vsz=20.0;

    vl=vdl-vol+(vol-vsl)*(k[1]/k[0]);
    vb=vdb-vob+(vob-vsb)*(k[1]/k[0]);
    v=sqrt(vl*vl+vb*vb);

    //sigma_y=290.0;
    //sigma_z=290.0;
    //sigma_y=-5.625*pow(xl*xl+yl*yl,0.5)+75.0;
    //sigma_z=-3.75*pow(xl*xl+yl*yl,0.5)+50.0;
    // lens-source relative transverse velocity
    //vlb=(k[5]*cos(l)-k[4]*sin(l))/cos(b);
    //vlb=k[5]*cos(l)-k[4]*sin(l);
    //vbb=k[6]*cos(b)-tan(b)*(k[4]*cos(b)*cos(l)+k[5]*cos(b)*sin(l));
    //vlb=-k[4];
    ///////////////////////////////////////////////////////////////////
    //vld=k[3];
    //vbd=k[6];
    //vlb=k[4];
    //vbb=k[5];
    //vl=(vld+vrot_disk-220.0)+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vl=k[3]+(220.0-vrot-vlb)*(k[1]/k[0]);
    //vb=vbd-vbb*(k[1]/k[0]);
    //v=pow(vl*vl+vb*vb,0.5);
    ////////////////////////////////////////////////////////////////////
    //vl=k[2];
    //vb=k[3];
    //mean_vl=(vrot_disk-220.0)+(220.0-vrot)*(k[1]/k[0]);
    //mean_vl=(220.0-vrot)*(k[1]/k[0]);
    //mean_vb=0.0;
    //sigma_vl=pow(sigma_y*sigma_y+82.5*82.5*k[1]*k[1]/(k[0]*k[0]),0.5);
    //sigma_vb=pow(sigma_z*sigma_z+66.3*66.3*k[1]*k[1]/(k[0]*k[0]),0.5);
    //v=pow(vl*vl+vb*vb,0.5);

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
    //vb=pow(fabs(184.2*365.0*365.0*Deff*Deff*M/(t*t)-vl*vl),0.5);

    // f(vld),f(vbd),f(vx),f(vy),f(vx)
    //vldfunc=(1.0/(pow(2.0*3.1415926,0.5)*30.0))*exp(-0.5*vld*vld/900.0);
    //vyfunc=(1.0/(pow(2.0*3.1415926,0.5)*82.5))*exp(-0.5*vlb*vlb/(82.5*82.5));
    //vzfunc=(1.0/(pow(2.0*3.1415926,0.5)*66.3))*exp(-0.5*vbb*vbb/(66.3*66.3));
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    //vbdfunc=(1.0/(pow(2.0*3.1415926,0.5)*20.0))*exp(-0.5*vbd*vbd/400.0);
    ////////////////////////////////////////////////////////////////////////////////////////
    //vlfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vl))*exp(-0.5*(vl-mean_vl)*(vl-mean_vl)/(sigma_vl*sigma_vl));
    //vbfunc=(1.0/(pow(2.0*3.1415926,0.5)*sigma_vb))*exp(-0.5*(vb-mean_vb)*(vb-mean_vb)/(sigma_vb*sigma_vb));


    vdxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdx))*exp(-0.5*(vdx*vdx)/(sigma_vdx*sigma_vdx));
    vdyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdy))*exp(-0.5*(vdy*vdy)/(sigma_vdy*sigma_vdy));
    vdzfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vdz))*exp(-0.5*(vdz*vdz)/(sigma_vdz*sigma_vdz));
    vsxfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsx))*exp(-0.5*(vsx*vsx)/(sigma_vsx*sigma_vsx));
    vsyfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsy))*exp(-0.5*(vsy*vsy)/(sigma_vsy*sigma_vsy));
    vszfunc=(1.0/(sqrt(2.0*3.1415926)*sigma_vsz))*exp(-0.5*(vsz*vsz)/(sigma_vsz*sigma_vsz));

    double re;
    re=4.28*pow(10.0,8.0)*Deff*pow(M,0.5);   // km

        F=(re/v)*2.77*pow(10.0,-8.0)*k[0]*k[0]*cos(b)*p_disk*p_pulsar*(Deff)*(sqrt(M))*(v*vdxfunc*vdyfunc*vdzfunc*vsxfunc*vsyfunc*vszfunc)*1.02*pow(10.0,-9.0);
    }

    return F;
}

int main (int argc, char *argv[])
{
    double timescale,re,event;

    // neutron stars, pulsar to bulge
    double res_x,err_x;

    double xl_x[10]={4.5,0.0,-2000.0,-2000.0,-2000.0,-550.0,-450.0,-350.0,-3.1415926,-3.1415926/2.0};
    double xu_x[10]={11.5,11.5,2000.0,2000.0,2000.0,550.0,450.0,350.0,3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_x;
    gsl_rng *r_x;

    gsl_monte_function G_x={&g_psr,10,0};

    size_t calls_x=500000000;

    gsl_rng_env_setup ();

    T_x=gsl_rng_default;

    // neutron stars, pulsar to disk 
    double res_xd,err_xd;

    double xl_xd[10]={0.0,0.0,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100.0,-3.1415926,-3.1415926/2.0};
    double xu_xd[10]={11.5,11.5,2000.0,2000.0,2000.0,100.0,150.0,100.0,3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_xd;
    gsl_rng *r_xd;

    gsl_monte_function G_xd={&g_psr_disk,10,0};

    size_t calls_xd=500000000;

    T_xd=gsl_rng_default;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    r_x=gsl_rng_alloc(T_x);
        
    gsl_monte_vegas_state *s_x = gsl_monte_vegas_alloc (10);

    r_xd=gsl_rng_alloc(T_xd);
        
    gsl_monte_vegas_state *s_xd = gsl_monte_vegas_alloc (10);

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    // timescale, re/v

    // neutron stars, pulsar to bulge
    double res_xr,err_xr;

    double xl_xr[10]={4.5,0.0,-2000.0,-2000.0,-2000.0,-550.0,-450.0,-350.0,-3.1415926,-3.1415926/2.0};
    double xu_xr[10]={11.5,11.5,2000.0,2000.0,2000.0,550.0,450.0,350.0,3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_xr;
    gsl_rng *r_xr;

    gsl_monte_function G_xr={&g_psrr,10,0};

    size_t calls_xr=500000000;

    T_xr=gsl_rng_default;

    // neutron stars, pulsar to disk 
    double res_xdr,err_xdr;

    double xl_xdr[10]={0.0,0.0,-2000.0,-2000.0,-2000.0,-100.0,-150.0,-100.0,-3.1415926,-3.1415926/2.0};
    double xu_xdr[10]={11.5,11.5,2000.0,2000.0,2000.0,100.0,150.0,100.0,3.1415926,3.1415926/2.0};

    const gsl_rng_type *T_xdr;
    gsl_rng *r_xdr;

    gsl_monte_function G_xdr={&g_psr_diskr,10,0};

    size_t calls_xdr=500000000;

    T_xdr=gsl_rng_default;
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    r_xr=gsl_rng_alloc(T_xr);
        
    gsl_monte_vegas_state *s_xr = gsl_monte_vegas_alloc (10);

    r_xdr=gsl_rng_alloc(T_xdr);
        
    gsl_monte_vegas_state *s_xdr = gsl_monte_vegas_alloc (10);

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    // calculation
    
        gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 10, 10000, r_x, s_x,&res_x, &err_x);

        do
        {
            gsl_monte_vegas_integrate (&G_x, xl_x, xu_x, 10, calls_x/5, r_x, s_x,&res_x, &err_x);
        }
        while (fabs (gsl_monte_vegas_chisq (s_x) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 10, 10000, r_xd, s_xd, &res_xd, &err_xd);

        do
        {
            gsl_monte_vegas_integrate (&G_xd, xl_xd, xu_xd, 10, calls_xd/5, r_xd, s_xd, &res_xd, &err_xd);
        }
        while (fabs (gsl_monte_vegas_chisq (s_xd) - 1.0) > 0.5);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

	// timescale, re/v

        gsl_monte_vegas_integrate (&G_xr, xl_xr, xu_xr, 10, 10000, r_xr, s_xr, &res_xr, &err_xr);

        do
        {
            gsl_monte_vegas_integrate (&G_xr, xl_xr, xu_xr, 10, calls_xr/5, r_xr, s_xr,&res_xr, &err_xr);
        }
        while (fabs (gsl_monte_vegas_chisq (s_xr) - 1.0) > 0.5);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

        gsl_monte_vegas_integrate (&G_xdr, xl_xdr, xu_xdr, 10, 10000, r_xdr, s_xdr, &res_xdr, &err_xdr);

        do
        {
            gsl_monte_vegas_integrate (&G_xdr, xl_xdr, xu_xdr, 10, calls_xdr/5, r_xdr, s_xdr, &res_xdr, &err_xdr);
        }
        while (fabs (gsl_monte_vegas_chisq (s_xdr) - 1.0) > 0.5);

        //////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////

        re = res_xr+res_xdr;
        event = res_x+res_xd;
        timescale = re/event;
        printf ("%e\n", timescale);

    gsl_monte_vegas_free (s_x);

    gsl_rng_free (r_x);

    gsl_monte_vegas_free (s_xd);

    gsl_rng_free (r_xd);
     
    /////////////////////////////////////////////////////////////////////////////////////////////////

    gsl_monte_vegas_free (s_xr);

    gsl_rng_free (r_xr);

    gsl_monte_vegas_free (s_xdr);

    gsl_rng_free (r_xdr);
     
    return 0;
}
