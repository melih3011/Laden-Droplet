#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

#define etam 1.e-3

double pressure(double eta,double rho0,double leps,double llam)
{
  double hc;

  hc=1.0-eta;
  return rho0*(1+eta*(1+eta*(1-eta)))/(hc*hc*hc)-4.*leps*pow(llam,3)*rho0*eta;
}

double sk0_bulk(double eta,double leps,double llam)
{
  double hc;

  hc=1.0-eta;
  return pow(hc,4)/(1.+4.*eta+4.*eta*eta-eta*eta*eta*(4.-eta)+leps*pow(llam,3)*(-8.*eta+32.*eta*eta-48.*pow(eta,3)+32.*pow(eta,4)-8.*pow(eta,5)));
}

double pressure_mixture(double n0,double n1,double n2,double n3)
{
  double hc;

  hc=1/(1.-n3);

  return hc*(n0+hc*(n1*n2+hc*(n2*n2*n2*(1.-n3/3.))/(12.*pi)));
}

double chemical_potential(double eta,double leps,double llam)
{
  double hc;

  hc=1.-eta;
  return eta*(8.+3.*(eta-3.)*eta)/(hc*hc*hc)-8.*leps*pow(llam,3)*eta;
}

double detadmu(double eta, double leps,double llam) {

    return eta*sk0_bulk(eta,leps,llam);
}

double chemical_potential_mixture(double n0,double n1,double n2,double n3,
				  double hR)
{
  double hc;

  hc=1.-n3;

  return -log(hc)+hR*n2/hc+pi_4*hR*hR*(n1/hc+n2*n2*(n3+hc*hc*log(hc))/
				       (12.*pi*hc*hc*n3*n3))+
    pi_4*hR*hR*hR/3.0*(n0/hc+n1*n2/(hc*hc)-
		       n2*n2*n2*(n3*(2.-5.*n3+n3*n3)+2*hc*hc*hc*log(hc))
		       /(36.*pi*hc*hc*hc*n3*n3*n3));
}

double chemical_potential_mixture_sw(double n0,double n1,double n2,double n3,double leps,double hR) {


  return -n3*leps+n0*leps*pi_4/3.*hR*hR*hR-n2*leps*hR+n1*leps*pi_4*hR*hR;

}

double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2)
{
  return (-n0)*log(1-n3);
}

double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2)
{
    return (n1*n2-nv2*nv1)/(1-n3);
}

double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2)
{
  double n13,hc;

  if (n3==0) return 0;

  n13=1.-n3;
  hc=n13*n13;

  return n2*(n2*n2-3*nv2*nv2)*(n3+hc*log(n13))/(36.0*pi*n3*n3*hc); 
}

double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,
        double hR)
{
  double d21,d22,d23,hc,ln3;
  
  ln3=1.-n3;
  hc=ln3*n3;
  
  if (n3<etam) {

    d21=n1/ln3+(n2*n2-nv2*nv2)/pi*(1./8.+2.*n3/9.+5.*n3*n3/16.+2.*n3*n3*n3/5.);
    d22=n2/ln3;
    d23=-log(ln3);
    
    return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);

  }

  d21=n1/ln3+(n2*n2-nv2*nv2)*(n3+ln3*ln3*log(ln3))/(12.*pi*hc*hc);
  d22=n2/ln3;
  d23=(-log(ln3));

  return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);
}

double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,
            double hR)
{
  double d13,d23,ln3,hc;

  ln3=1.-n3;
  hc=ln3*n3;

  d13=n0/ln3+(n1*n2-nv1*nv2)/(ln3*ln3);

  if (n3<etam) {
    d23=n2*(n2*n2-3.*nv2*nv2)/pi*(2./27.+5.*n3/48.+2.*n3*n3/15.+35.*n3*n3*n3/54.);
    return d13+d23;
  }

  d23=-n2*(n2*n2-3.*nv2*nv2)*(n3*(2.-5.*n3+n3*n3)+2.*ln3*ln3*ln3*log(ln3))/(36.*pi*hc*hc*hc);

  return d13+d23;

}

double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,
             double hR)
{
  double dv21,dv22,dv2,ln3,hc;

  ln3=1.-n3;

  hc=ln3*n3;

  if (n3<etam) {

    dv21=nv1/ln3+2.*n2*nv2/pi*(1./8.+2.*n3/9.+5.*n3*n3/16.+2.*n3*n3*n3/5.);
    dv22=nv2/ln3;

    return -(dv21+dv22/(pi_4*hR));

  }

  dv21=nv1/ln3+n2*nv2*(n3+ln3*ln3*log(ln3))/(6.*pi*hc*hc);

  dv22=nv2/ln3;

  return -(dv21+dv22/(pi_4*hR));
}

double myd3_sw(double n0,double n1,double n2,double n3,double nv1,double nv2,double hR) {

  return n0;
} 

double myd2_sw(double n0,double n1,double n2,double n3,double nv1,double nv2,double hR) {

  return n1-n2/(pi_4*hR)-n3/(pi_4*hR*hR);
} 

double mydv2_sw(double n0,double n1,double n2,double n3,double nv1,double nv2,double hR) {

  return nv1-nv2/(pi_4*hR);
} 


void filename(char *fln)
{
  sprintf(fln,"%s","wb");
}
