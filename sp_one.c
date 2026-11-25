#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769  /* 2*pi */
#define pi_4 12.5663706143591729539 /* 4*pi */

#define R 50
#define R2 (int)(1 * R)
#define N (800 * R)

#define lam1 1.5
#define eps1 1.2

#define lam2 1.5
#define eps2 0.0

#define eps12 1.72518//sqrt(eps1*eps2)

#define lo1 (int)(2*R*lam1)
#define lo2 (int)(2*R2*lam2)

#define lo12 (int)((double)(lo1+lo2)*0.5)

//eps1 1.5
//#define rho1l (3./pi_4*0.44272121)
//#define rho1v (3./pi_4*0.00097083)

//eps1 1.2
#define rho1l (3./pi_4*0.38022500)
#define rho1v (3./pi_4*0.00586050)

//eps1 1.0
//#define rho1l (3./pi_4*0.3150217436)
//#define rho1v (3./pi_4*0.0193657384)

//eps1 5/6, i.e. closer to critical point
//#define rho1l (3./pi_4*0.21679094074)
//#define rho1v (3./pi_4*0.06258312516)

//eps1 1/1.2235
//#define rho1l (3./pi_4*0.20032127644643)
//#define rho1v (3./pi_4*0.0730433213868317)


#define Nn 992513

double rho02l;

int Rb, Rn;
int Rmax = lo2 > lo1 ? lo2 : lo1;

double n3[N], n2[N], n1[N], n0[N], n2v[N], n1v[N], nsw1[N], nsw2[N], nsw12[N], nsw21[N];

double rho1[N];
double w3[2*R+1][2*R+1], w2[2*R+1][2*R+1], w2v[2*R+1][2*R+1];

double rho2[N];
double w32[2*R2+1][2*R2+1], w22[2*R2+1][2*R2+1], w2v2[2*R2+1][2*R2+1];

double wsw1[2*lo1+1][2*lo1+1], wsw2[2*lo2+1][2*lo2+1], wsw12[2*lo12+1][2*lo12+1];

double phi[N], phisw[N], phiid[N], phiid2[N], dphi[N], dphi2[N], dep[N], dep2[N], ads[N], Vext[N];

double muid, mu1l_coex, mu1l;
double dr = 1.0/(1.0*R), a1, a2, fac;

FILE *fp;

double get_magnitude(const double *arr) {

  int i;
  double max = fabs(arr[0]);

  for (i = 1; i < N-5*Rmax; ++i) {
    if (fabs (arr[i]) > max) max = fabs(arr[i]);
  }

  return max;
}

int get_radius(const double *arr, int mode) {

  double mid, tol;
  int rad = 0;

  if (mode == 0) {

    tol = 1.e-3;
    mid = 1./3.*arr[0]+2./3.*arr[N-5*Rmax];

    while (fabs(arr[rad]-mid) > tol) rad++;
  
  } else if (mode == 1) {

    tol = 1.e-6;
    while (fabs(arr[rad]-arr[N-5*Rmax]) > tol) rad++;
  
  } else exit(-1);

  return rad;

}

double fr2int(const double *f, double ldr, int n, double a) {
  double val=0.0, hr1=6.*a*a, hr2=4.*a*ldr, hr3=ldr*ldr;
  register int j;

  for (j = 0; j < n; j++) {
    val += f[j]*(hr1+hr2*(1.+3.*j)+hr3*(1.+4.*j+6.*j*j)) +
        f[j+1]*(hr1+hr2*(2.+3.*j)+hr3*(3.+8.*j+6.*j*j));
  }
  
  return (ldr*val/12.0);
}

double uni_random() {

  double sign = pow(-1., rand()%2);

  return sign*(double)rand()/RAND_MAX;
}

void print_results() {

  void filename(char *fln);

  char fln[250], fname[50];
  
  filename(fname);

  sprintf(fln, "ld_%.2f_%.2f_%.4f", (double)(R*dr/(R2*dr)), N*dr/2., mu1l+muid);
  fp = fopen(fln, "w");
  
  fprintf(fp, "# eps1 = %.3f, R11 = %.3f\n", eps1, lo1*dr);
  fprintf(fp, "# eps2 = %.3f, R22 = %.3f\n", eps2, lo2*dr);
  fprintf(fp, "# eps12 = %.3f, R12 = %.3f\n", eps12, lo12*dr);
  fprintf(fp, "# a1 = %.2e, a2 = %.2e, fac = %.2e\n", a1, a2, fac);
  fprintf(fp, "# Nn = %d, mu = %.5f\n", Nn, mu1l+muid);
  fprintf(fp, "# xm = %.3f, xl = %.3f, xn = %.3f\n", (N-3*Rmax)*dr/2., get_radius(rho1, 1)*dr/2., get_radius(rho2, 1)*dr/2.);
  fprintf(fp, "# xi = %.6f\n", pi_4*fr2int(rho2, dr, get_radius(rho2, 1), 0.)/Nn);

  for (int i = 0; i < N-5*Rmax; i++) {
    fprintf(fp, "%10.5lf  %.10e %.10e %.10e %.10e\n", i*dr/2., rho1[i], rho2[i], dphi[i], n3[i]);
  }
  fclose(fp);

}

void print_file(int j) {

  void filename(char *fln);
  char fln[250], command[300];

  sprintf(fln, "temp.dat");
  fp = fopen(fln, "w");

  for (int i = 0; i < N-5*Rmax; i++) {
    fprintf(fp, "%10.5lf  %.10e %.10e %.10e %.10e\n", i*dr/2., rho1[i], rho2[i], dphi[i], exp(-dep2[i]));
  }

  fclose(fp);

  sprintf(command, "mv %s ./temp", fln);
  system(command);

}

double sprod(const double *arr1, const double *arr2) {

  int i;
  double prod[N];

  for (i=0; i<N-3*Rmax; i++) prod[i] = arr1[i]*arr2[i];

  return pi_4*fr2int(prod, dr, N-4*Rmax, 0.0);

}

void param_init() {

  double q = 1.*R/R2, eta20;
  rho02l = (double)(3.*Nn/pi_4/pow(Rn*dr, 3));
  eta20 = rho02l*pi_4/3.*pow(R2*dr, 3);
  printf("\n***************\n");
  printf("Initialized %d particles with packing fraction %.4e\n", Nn, eta20);
  printf("Initial size Rb=%.2f and Rn=%.2f\n", Rb*dr/2., Rn*dr/2.);
  printf("Square-well energies: eps1=%.3f, eps2=%.3f, eps12=%.3f\n", eps1, eps2, eps12);
  printf("Square-well ranges: R11=%.3f, R22=%.3f, R12=%.3f\n", lo1*dr, lo2*dr, lo12*dr);
  printf("***************\n\n");
  
  muid = log(pi_4*rho1l/3.);

  //eps1 1.0
  //mu1l=-4.30639-muid;
  //mu1l=-4.34874-muid; //95%
  //mu1l=-4.39396-muid; // 90%
  //mu1l=-4.44237-muid; //85%
  //mu1l=-4.4731-muid; //82%
  //mu1l=-4.49434-muid; //80%
  //mu1l=-4.51622-muid; //78%
  //mu1l=-4.55032-muid; //75%
  //mu1l=-4.61083-muid; //70%
  //mu1l=-4.67654-muid; //65%
  mu1l=-4.91418-muid; //50%
  //mu1l=-5.01144-muid; //45%
  //mu1l=-5.16929-muid; //38%
  //mu1l=-5.78284-muid; //20%
  //mu1l=-6.46058-muid; //10%

  //eps1 5/6
  //mu1l=-3.612-muid; //>100%
  //mu1l=-3.6135-muid; //100%
  //mu1l=-3.64487-muid; //95%
  //mu1l=-3.76173-muid; //80%

  //eps1 1/1.2235
  //mu1l=-3.55377-muid; //100%
  //mu1l=-3.61681-muid; //90%
  //mu1l=-3.69603-muid; //80%


}

void init_matrix(int n, double wm[][2 * n + 1], const double *lwm) {

  int i, j;
  double p1[2], p2[3], p3[4], p4[5];

  p1[0] = 0.5; p1[1] = 0.5;
  p2[0] = 1./3.; p2[1] = 4./3.; p2[2] = 1./3.;
  p3[0] = 3./8.; p3[1] = 9./8.; p3[2] = 9./8.; p3[3] = 3./8.;
  p4[0] = 14./45.; p4[1] = 64./45.; p4[2] = 24./45.; p4[3] = 64./45.; p4[4] = 14./45.;

  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 2 * n; ++j) {
      wm[i][j] = 0.0;
    }
  }

  for (j = 0; j < 2; j++) wm[1][j] = lwm[j]*p1[j];

  for (j = 0; j < 3; j++) wm[2][j] = lwm[j]*p2[j];

  for (j = 0; j < 4; j++) wm[3][j] = lwm[j]*p3[j];
  
  for (j = 0; j < 5; j++) wm[4][j] = lwm[j]*p4[j];

  for (i = 5; i <= 2 * n; i++) {
    for (j = 0; j < 2 * n - i; j++) {
      wm[i][2 * n - j] = 0.0;
    }

    wm[i][i] = 3./8.*lwm[i];
    wm[i][i-1] = 7./6.*lwm[i-1];
    wm[i][i-2] = 23./24.*lwm[i-2];

    for (j = 3; j <= i-3; j++) wm[i][j] = lwm[j];

    wm[i][2] = 23./24.*lwm[2];
    wm[i][1] = 7./6.*lwm[1];
    wm[i][0] = 3./8.*lwm[0];
  }
}

void data_init() {

  double lw3[2*R+1], lw2[2*R+1], lw2v[2*R+1];
  double lw32[2*R2+1], lw22[2*R2+1], lw2v2[2*R2+1];
  double lwsw1[2*lo1+1], lwsw2[2*lo2+1], lwsw12[2*lo12+1];

  double z, lR1 = R*dr, lR2 = R2*dr, llo1 = 2.*lR1*lam1, llo2 = 2.*lR2*lam2, llo12 = lo12*dr;
  int i;

  param_init();

  for (i = 0; i <= 2*R; i++) {
    z = (i-R)*dr;
    lw3[i] = pi*dr*(lR1*lR1-z*z);
    lw2[i] = pi_2*dr*lR1;
    lw2v[i] = pi_2*dr*z;
  }

  for (i = 0; i <= 2*R2; i++) {
    z = (i-R2)*dr;
    lw32[i] = pi*dr*(lR2*lR2-z*z);
    lw22[i] = pi_2*dr*lR2;
    lw2v2[i] = pi_2*dr*z;
  }

  for (i=0; i <= 2*lo1; ++i) {
    z = (i-lo1)*dr;
    lwsw1[i] = -eps1*dr*pi*(llo1*llo1-z*z);

  }

  for (i=0; i <= 2*lo2; ++i) {
    z = (i-lo2)*dr;
    lwsw2[i] = -eps2*dr*pi*(llo2*llo2-z*z);
  }

  for (i=0; i <= 2*lo12; ++i) {
    z = (i-lo12)*dr;
    lwsw12[i] = -eps12*dr*pi*(llo12*llo12-z*z);
  }


  init_matrix(R, w3, lw3);
  init_matrix(R, w2, lw2);
  init_matrix(R, w2v, lw2v);

  init_matrix(R2, w32, lw32);
  init_matrix(R2, w22, lw22);
  init_matrix(R2, w2v2, lw2v2);

  init_matrix(lo1, wsw1, lwsw1);
  init_matrix(lo2, wsw2, lwsw2);
  init_matrix(lo12, wsw12, lwsw12);

}

void weighted_densities() {

  int i, j, jp, jm;

  double lR1 = R*dr, lR2 = R2*dr;

  for (i=0; i<N; i++) {

    n3[i]=0.0, n2[i]=0.0, n1[i]=0.0, n0[i]=0.0, n1v[i]=0.0, n2v[i]=0.0;
    nsw1[i]=0.0, nsw2[i]=0.0, nsw12[i]=0.0, nsw21[i]=0.0;

  }

  n3[0] = pi_4*fr2int(rho1, dr, R, 0.)+pi_4*fr2int(rho2, dr, R2, 0.);
  n2[0] = pi_4*pow(lR1, 2)*rho1[R]+pi_4*pow(lR2, 2)*rho2[R2];
  n1[0] = lR1*rho1[R]+lR2*rho2[R2];
  n0[0] = rho1[R]+rho2[R2];
  n2v[0] = 0.0;
  n1v[0] = 0.0;

  nsw1[0] = -eps1*pi_4*fr2int(rho1, dr, lo1, 0.);
  nsw2[0] = -eps2*pi_4*fr2int(rho2, dr, lo2, 0.);
  nsw12[0] = -eps12*pi_4*fr2int(rho2, dr, lo12, 0.);
  nsw21[0] = -eps12*pi_4*fr2int(rho1, dr, lo12, 0.);


  //1st component HS
  //=======================================================================
  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = 1; i < R; i++) {

    double hr = i*dr, hrho;

    double n3i = 0.0, n2i = 0.0, nv2i = 0.0, n0i, n1i, nv1i;

    jp = R+i;
    jm = R-i;

    for (j = jm; j <= jp; j++) {
      hrho = rho1[j]*(j*dr);
      n3i += hrho*w3[jp-jm][R+i-j];
      n2i += hrho*w2[jp-jm][R+i-j];
      nv2i += hrho*w2v[jp-jm][R+i-j];
    }

    n1i = n2i/(pi_4*lR1);
    n0i = n1i/lR1;
    nv1i = nv2i/(pi_4*lR1);

    n3[i] += pi_4*fr2int(rho1,dr,R-i,0.0)+n3i/hr;
    n2[i] += n2i/hr;
    n1[i] += n1i/hr;
    n0[i] += n0i/hr;
    n2v[i] += (n3i/hr+nv2i)/hr;
    n1v[i] += (n3i/lR1/(hr*pi_4)+nv1i)/hr;    

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = R; i < N-Rmax; i++) {

    double hr = i*dr, hrho;

    double n3i = 0.0, n2i = 0.0, nv2i = 0.0, n0i, n1i, nv1i;

    jp = i+R;
    jm = i-R;

    for (j = jm; j <= jp; j++) {
      hrho = rho1[j]*(j*dr);
      n3i += hrho*w3[jp-jm][R+i-j];
      n2i += hrho*w2[jp-jm][R+i-j];
      nv2i += hrho*w2v[jp-jm][R+i-j];
    }

    n1i = n2i/(pi_4*lR1);
    n0i = n1i/lR1;
    nv1i = nv2i/(pi_4*lR1);

    n3[i] += n3i/hr;
    n2[i] += n2i/hr;
    n1[i] += n1i/hr;
    n0[i] += n0i/hr;
    n2v[i] += (n3i/hr+nv2i)/hr;
    n1v[i] += (n3i/lR1/(hr*pi_4)+nv1i)/hr;    

  }
  //=======================================================================

  //2nd component HS
  //=======================================================================
  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = 1; i < R2; i++) {

    double hr = i*dr, hrho;
    double m3i = 0.0, m2i = 0.0, mv2i = 0.0, m0i, m1i, mv1i;

    jp = R2+i;
    jm = R2-i;

    for (j = jm; j <= jp; j++) {
      hrho = rho2[j]*(j*dr);
      m3i += hrho*w32[jp-jm][R2+i-j];
      m2i += hrho*w22[jp-jm][R2+i-j];
      mv2i += hrho*w2v2[jp-jm][R2+i-j];
    }

    m1i = m2i/(pi_4*lR2);
    m0i = m1i/lR2;
    mv1i = mv2i/(pi_4*lR2);

    n3[i] += pi_4*fr2int(rho2,dr,R2-i,0.0)+m3i/hr;
    n2[i] += m2i/hr;
    n1[i] += m1i/hr;
    n0[i] += m0i/hr;
    n2v[i] += (m3i/hr+mv2i)/hr;
    n1v[i] += (m3i/lR2/(hr*pi_4)+mv1i)/hr;    

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = R2; i < N-Rmax; i++) {

    double hr = i*dr, hrho;
    double m3i = 0.0, m2i = 0.0, mv2i = 0.0, m0i, m1i, mv1i;

    jp = i+R2;
    jm = i-R2;

    for (j = jm; j <= jp; j++) {
      hrho = rho2[j]*(j*dr);
      m3i += hrho*w32[jp-jm][R2+i-j];
      m2i += hrho*w22[jp-jm][R2+i-j];
      mv2i += hrho*w2v2[jp-jm][R2+i-j];
    }

    m1i = m2i/(pi_4*lR2);
    m0i = m1i/lR2;
    mv1i = mv2i/(pi_4*lR2);

    n3[i] += m3i/hr;
    n2[i] += m2i/hr;
    n1[i] += m1i/hr;
    n0[i] += m0i/hr;
    n2v[i] += (m3i/hr+mv2i)/hr;
    n1v[i] += (m3i/lR2/(hr*pi_4)+mv1i)/hr;    

  }
  //=======================================================================

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=1; i<lo1; i++) {

    double sw1 = 0.0, hr=i*dr, hrho;

    jp = lo1+i;
    jm = lo1-i;

    for (j = jm; j <= jp; j++) {
      hrho = rho1[j]*(j*dr);
      sw1 += hrho*wsw1[jp-jm][lo1+i-j];
    }

    nsw1[i] += -eps1*pi_4*fr2int(rho1,dr,lo1-i,0.0)+sw1/hr;

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=lo1; i<N-Rmax; i++) {

    double sw1 = 0.0, hr=i*dr, hrho;

    jp = i+lo1;
    jm = i-lo1;

    for (j = jm; j <= jp; j++) {
      hrho = rho1[j]*(j*dr);
      sw1 += hrho*wsw1[jp-jm][lo1+i-j];
    }

    nsw1[i] += sw1/hr;

  }


  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=1; i<lo2; i++) {

    double sw2 = 0.0, hr=i*dr, hrho;

    jp = lo2+i;
    jm = lo2-i;

    for (j = jm; j <= jp; j++) {
      hrho = rho2[j]*(j*dr);
      sw2 += hrho*wsw2[jp-jm][lo2+i-j];
    }

    nsw2[i] += -eps2*pi_4*fr2int(rho2,dr,lo2-i,0.0)+sw2/hr;

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=lo2; i<N-Rmax; i++) {

    double sw2 = 0.0, hr=i*dr, hrho;

    jp = i+lo2;
    jm = i-lo2;

    for (j = jm; j <= jp; j++) {
      hrho = rho2[j]*(j*dr);
      sw2 += hrho*wsw2[jp-jm][lo2+i-j];
    }

    nsw2[i] += sw2/hr;

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=1; i<lo12; i++) {

    double sw12 = 0.0, sw21 = 0.0, hr=i*dr;

    jp = lo12+i;
    jm = lo12-i;

    for (j = jm; j <= jp; j++) {
      sw12 += rho2[j]*j*dr*wsw12[jp-jm][lo12+i-j];
      sw21 += rho1[j]*j*dr*wsw12[jp-jm][lo12+i-j];
    }

    nsw12[i] += -eps12*pi_4*fr2int(rho2,dr,lo12-i,0.0)+sw12/hr;
    nsw21[i] += -eps12*pi_4*fr2int(rho1,dr,lo12-i,0.0)+sw21/hr;

  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i=lo12; i<N-Rmax; i++) {

    double sw12 = 0.0, sw21 = 0.0, hr=i*dr;

    jp = i+lo12;
    jm = i-lo12;

    for (j = jm; j <= jp; j++) {
      sw12 += rho2[j]*j*dr*wsw12[jp-jm][lo12+i-j];
      sw21 += rho1[j]*j*dr*wsw12[jp-jm][lo12+i-j];
    }

    nsw12[i] += sw12/hr;
    nsw21[i] += sw21/hr;

  }
	  
}

void derivative() {
  
  double myd2(double n0, double n1, double n2, double n3, double nv1, double nv2, double hR);
  double myd3(double n0, double n1, double n2, double n3, double nv1, double nv2, double hR);
  double mydv2(double n0, double n1, double n2, double n3, double nv1, double nv2, double hR);

  int i, j, jp, jm;
  double d2[N], d22[N], d3[N], dv2[N], dv22[N];
  double lR1 = R*dr, lR2 = R2*dr;

  #pragma omp parallel for private(i) num_threads(4)
  for (i = 0; i < N-2*Rmax; i++) {

    double ln3 = n3[i], ln2 = n2[i], ln1 = n1[i];
    double ln0 = n0[i], lnv2 = n2v[i], lnv1 = n1v[i];

    d2[i] = myd2(ln0, ln1, ln2, ln3, lnv1, lnv2, lR1);
    d3[i] = myd3(ln0, ln1, ln2, ln3, lnv1, lnv2, lR1);
    dv2[i] = mydv2(ln0, ln1, ln2, ln3, lnv1, lnv2, lR1);

    d22[i] = myd2(ln0, ln1, ln2, ln3, lnv1, lnv2, lR2);
    dv22[i] = mydv2(ln0, ln1, ln2, ln3, lnv1, lnv2, lR2);

  }

  dep[0] = pi_4*fr2int(d3, dr, R, 0.)+pi_4*pow(lR1, 2)*(d2[R]+dv2[R])+
           nsw1[0]+nsw12[0]-mu1l;

  dep2[0] = pi_4*fr2int(d3, dr, R2, 0.)+pi_4*pow(lR2, 2)*(d22[R2]+dv22[R2])+
            nsw2[0]+nsw21[0];

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = 1; i < R; i++) {

    double hr = i*dr, hrs;
    double c3i = 0.0, c2i = 0.0, cv2i = 0.0;

    jp = R+i;
    jm = R-i;

    for (j = jm; j <= jp; j++) {
      hrs = j*dr;
      c3i += d3[j]*hrs*w3[jp-jm][(R+i)-j];
      c2i += d2[j]*hrs*w2[jp-jm][(R+i)-j];
      cv2i += dv2[j]*(w3[jp-jm][(R+i)-j]-hrs*w2v[jp-jm][(R+i)-j]);
    }

    dep[i] = pi_4*fr2int(d3,dr,R-i,0.0)+(c3i+c2i+cv2i)/hr+nsw1[i]+nsw12[i]-mu1l;
    dphi[i] = dep[i]+log(rho1[i]/rho1l);

  }
    
  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = R; i < N-2*Rmax; i++) {

    double hr = i*dr, hrs;
    double c3i = 0.0, c2i = 0.0, cv2i = 0.0;

    jp = i+R;
    jm = i-R;

    for (j = jm; j <= jp; j++) {
      hrs = j*dr;
      c3i += d3[j]*hrs*w3[jp-jm][(R+i)-j];
      c2i += d2[j]*hrs*w2[jp-jm][(R+i)-j];
      cv2i += dv2[j]*(w3[jp-jm][(R+i)-j]-hrs*w2v[jp-jm][(R+i)-j]);
    }

    dep[i] = (c3i+c2i+cv2i)/hr+nsw1[i]+nsw12[i]-mu1l;
    dphi[i] = dep[i]+log(rho1[i]/rho1l);

  }    
    
  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = 1; i < R2; i++) {

    double hr = i*dr, hrs;
    double e3i = 0.0, e2i = 0.0, ev2i = 0.0;

    jp = R2+i;
    jm = R2-i;

    for (j = jm; j <= jp; j++) {
      hrs = j*dr;
      e3i += d3[j]*hrs*w32[jp-jm][(R2+i)-j];
      e2i += d22[j]*hrs*w22[jp-jm][(R2+i)-j];
      ev2i += dv22[j]*(w32[jp-jm][(R2+i)-j]-hrs*w2v2[jp-jm][(R2+i)-j]);
    }

    dep2[i] = pi_4*fr2int(d3,dr,R2-i,0.0)+(e3i+e2i+ev2i)/hr+nsw2[i]+nsw21[i];
    dphi2[i] = dep2[i]+log(rho2[i]/rho02l);
    
  }

  #pragma omp parallel for private(i, j, jm, jp) num_threads(4)
  for (i = R2; i < N-Rmax; i++) {

    double hr = i*dr, hrs;
    double e3i = 0.0, e2i = 0.0, ev2i = 0.0;

    jp = i+R2;
    jm = i-R2;

    for (j = jm; j <= jp; j++) {
      hrs = j*dr;
      e3i += d3[j]*hrs*w32[jp-jm][(R2+i)-j];
      e2i += d22[j]*hrs*w22[jp-jm][(R2+i)-j];
      ev2i += dv22[j]*(w32[jp-jm][(R2+i)-j]-hrs*w2v2[jp-jm][(R2+i)-j]);
    }

    dep2[i] = (e3i+e2i+ev2i)/hr+nsw2[i]+nsw21[i];
    dphi2[i] = dep2[i]+log(rho2[i]/rho02l);
  }

}

void get_external_potential(double lrho) {

  int i = 0;

  for (i=0; i<N-3*Rmax; i++) rho1[i] = 0.0;

  weighted_densities();
  derivative();

  for (i=0; i<N-3*Rmax; i++) Vext[i] = -log(rho2[i]/lrho)-dep2[i];

}

void minimize() {
  
  double hw = 0.0, hw_old = 1.0, hw_eps, c, dN, magn;
  double rho1_old[N], rho2_old[N];
  int i, j = 0, ITS = 250000;

  hw_eps = 1.e-18;

  system("rm ./temp/*.*");

  for (i = 0; i < Rb; i++) rho1[i] = rho1l;
  for (i = Rb; i < N-2*Rmax; i++) rho1[i] = rho1v;

  for (i = 0; i <= Rn; i++) rho2[i] = rho02l;
  for (i = Rn+1; i < N-2*Rmax; i++) rho2[i] = 0.0;
  
  while ((fabs(hw-hw_old) > hw_eps) && (j < ITS)) {

    hw_old = hw;

    weighted_densities();

    derivative();

    a2 = fac*exp(-get_magnitude(dep2));
    
    #pragma omp parallel for private(i)
    for (i = 0; i < N-3*Rmax; i++) {
      
      rho1_old[i] = rho1[i];
      rho2_old[i] = rho2[i];

      rho1[i] = (1.-a1)*rho1[i]+a1*rho1l*exp(-dep[i]);
      rho2[i] = (1.-a2)*rho2[i]+a2*rho02l*exp(-dep2[i]);
      
    }

    c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0);

    #pragma parallel for reduction (+:hw)
    for (i = 0; i < N-3*Rmax; ++i) {
      rho2[i] *= (double)Nn/c;
      hw += 0.5*(rho1_old[i]-rho1[i])*(rho1_old[i]-rho1[i])+0.5*(rho2_old[i]-rho2[i])*(rho2_old[i]-rho2[i]);
    }

    dN = (double)Nn-pi_4*fr2int(rho2, dr, N-3*Rmax, 0.);
    magn = pi_4*get_magnitude(rho1)/3.;

    print_file(j);

    printf("\r%.2f %% ... %.2e ... %.2e ... %.4e ...  %.2f ... %.2f ... %.4f", 100.*j/ITS, get_magnitude(dphi),
     fabs(hw-hw_old), magn, get_radius(rho1, 1)*dr/2., get_radius(rho2, 1)*dr/2., pi_4*fr2int(rho2, dr, get_radius(rho2, 1), 0.)/Nn);
    fflush(stdout);
    
    j++;
  }

  printf("\ndOm/drho=%.3e, dN=%.3e\n", get_magnitude(dphi), dN);

  if (j >= ITS) printf("Exceeded maximum number of iterations.\n");
  else printf("Did %d iterations\n", j);
}

void minimize_eff(double x, double b) {

  double hw = 0.0, hw_old = 1.0, hw_eps, c, dN, magn, a, q, rl=Rb*dr, rn=Rn*dr;
  double c1, c2, a11, a22, a12, a01, a02, hr, dphi_max, dphi2_max;
  double rho1_old[N], rho2_old[N], rhon1[3][N], rhon2[3][N], dn1[3][N], dn2[3][N], d01[N], d02[N];
  int i, j = 0, k, ITS = 10000, ITS_picard = 1000, xi = (int)(x/dr);
  
  q = (double)(R*dr/(R2*dr));
  a = -eps12*rho02l*rho1l*pow(lam1*q+lam2, 3)*pi_4/3.*pow(R*dr, 3)*pi_4/3.*pow(R2*dr, 3);

  for (i=0; i<N-3*Rmax; i++) Vext[i] = a + (b-a)/(1.+exp(-1.1*(i*dr-x)));

  hw_eps = 1.e-19;

  system("rm ./temp/*.*");

  //Stage 1: Minimize 2nd component first using an external potential best suited
  // In general, Vext resembles a pot being attractive inside and repulsive outside the droplet 
  //============================================================
  printf("Beginning to minimize 2.component with Vext\n");
  for (i=0; i<N-2*Rmax; i++) rho1[i] = 0.0;
  for (i = 0; i <= Rn; i++) rho2[i] = rho02l;
  for (i = Rn+1; i < N-2*Rmax; i++) rho2[i] = 0.0; 
  
  while ((fabs(hw-hw_old) > hw_eps) && (j < ITS)) {
    hw_old = hw;
    weighted_densities();
    derivative();
    #pragma omp parallel for private(i)
    for (i = 0; i < N-3*Rmax; i++) {
      rho2_old[i] = rho2[i];
      rho2[i] = (1.-0.1)*rho2[i]+0.1*rho02l*exp(-dep2[i]-Vext[i]);
    }
    c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0);
    #pragma parallel for reduction (+:hw)
    for (i = 0; i < N-3*Rmax; ++i) {
      rho2[i] *= (double)Nn/c;
      hw += (rho2_old[i]-rho2[i])*(rho2_old[i]-rho2[i]);
    }    
    j++;
  }
  //============================================================

  hw = 0.0;
  hw_old = 1.0;
  j = 0;

  //Stage 2: Now, minimize only the 1st component by keeping the density profile of the 2nd component constant. 
  //============================================================
  printf("Beginning to minimize 1.component\n");
  for (i=0; i<Rb; i++) rho1[i] = rho1l;  
  for (i=Rb; i<N-2*Rmax; i++) rho1[i] = rho1v;  

  while ((fabs(hw-hw_old) > hw_eps) && (j < ITS)) {
    hw_old = hw;
    weighted_densities();
    derivative();
    #pragma omp parallel for private(i)
    for (i = 0; i < N-3*Rmax; i++) {
      rho1_old[i] = rho1[i];
      rho1[i] = (1.-0.1)*rho1[i]+0.1*rho1l*exp(-dep[i]);
    }
    #pragma parallel for reduction (+:hw)
    for (i = 0; i < N-3*Rmax; ++i) {
      hw += (rho1_old[i]-rho1[i])*(rho1_old[i]-rho1[i]);
    }    
    j++;
  }
//============================================================

  hw = 0.0;
  hw_old = 1.0;
  j = 0;

  //Stage 3: Perform usual Picard-Iteration (or Ng version) of the pre-modified density profiles.
  // Herewith, the overall convergence against an equilibrium solution should be accelerated.
  //============================================================
  printf("Starting Picard-Iteration\n"); 
  while ((fabs(hw-hw_old) > hw_eps) && (j < ITS_picard)) {
  //while (get_magnitude(dphi) > 5.e-5) {

    hw_old = hw;

    weighted_densities();

    derivative();

    //a2 = fac;
    a2 = fac*exp(-get_magnitude(dep2));
    //rho02l = 3./pow(Rb*dr,3)*fr2int(rho2, dr, Rb, 0.0);

    #pragma omp parallel for private(i)
    for (i = 0; i < N-3*Rmax; i++) {  
      rho1_old[i] = rho1[i];
      rho2_old[i] = rho2[i];
      
      rho1[i] = (1.-a1)*rho1[i]+a1*rho1l*exp(-dep[i]);
      rho2[i] = (1.-a2)*rho2[i]+a2*rho02l*exp(-dep2[i]);
      
    }

    c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0);

    #pragma parallel for reduction (+:hw)
    for (i = 0; i < N-3*Rmax; ++i) {
      rho2[i] *= (double)Nn/c;
      hw += 0.5*(rho1_old[i]-rho1[i])*(rho1_old[i]-rho1[i])+0.5*(rho2_old[i]-rho2[i])*(rho2_old[i]-rho2[i]);
    }

    rl = get_radius(rho1, 1)*dr;
    rn = get_radius(rho2, 1)*dr;

    dN = (double)Nn-pi_4*fr2int(rho2, dr, N-3*Rmax, 0.);
    magn = get_magnitude(n3);

    j++;
  }
  //============================================================

  hw = 0.0;
  hw_old = 1.0;
  j = 0;

  //Stage 4: Finally, perform Ng of the premodified density profiles.
  // Herewith, the overall convergence against an equilibrium solution should be accelerated.
  //============================================================
  printf("Starting Ng-algorithm\n"); 

  while ((fabs(hw-hw_old) > hw_eps) && (j < 250000)) {

    
    hw_old = hw;
    
    a2 = fac*exp(-get_magnitude(dep2));

    for (k=2; k>=0; k--) {

      weighted_densities();

      derivative();

      dphi_max = get_magnitude(dphi);
      dphi2_max = get_magnitude(dphi2);

      #pragma omp parallel for private(i)
      for (i=0; i<N-3*Rmax; i++) {
        rhon1[k][i] = (1.-a1)*rho1[i]+a1*rho1l*exp(-dep[i]);
        rhon2[k][i] = (1.-a2)*rho2[i]+a2*rho02l*exp(-dep2[i]);
      }
      
      c = pi_4*fr2int(rhon2[k], dr, N-3*Rmax, 0.0); 
      #pragma parallel for num_threads(4)
      for (i = 0; i < N-3*Rmax; ++i) rhon2[k][i] *= (double)Nn/c;

      #pragma omp parallel for private(i)
      for (i=0; i<N-3*Rmax; i++) {
        dn1[k][i] = rhon1[k][i] - rho1[i];
        dn2[k][i] = rhon2[k][i] - rho2[i];

        rho1[i] = rhon1[k][i];
        rho2[i] = rhon2[k][i];
      }
    }

      for (i=0; i<N-3*Rmax; i++) {
        d01[i] = dn1[0][i]-dn1[1][i];
        d02[i] = dn1[0][i]-dn1[2][i];
      }
      
      a11 = sprod(d01, d01);
      a22 = sprod(d02, d02);
      a12 = sprod(d01, d02);
      a01 = sprod(dn1[0], d01);
      a02 = sprod(dn1[0], d02);
      hw += 0.5*sprod(dn1[0], dn1[0]);

      hr = a12*a12 - a11*a22;
      c1 = (a02*a12 - a01*a22)/hr;
      c2 = (a01*a12 - a02*a11)/hr;

      for (i=0; i<N-3*Rmax; i++) {
        rho1[i] = (1.-c1-c2)*rhon1[0][i] + c1*rhon1[1][i] + c2*rhon1[2][i];
      }

//=============================================================================
      for (i=0; i<N-3*Rmax; i++) {
        d01[i] = dn2[0][i]-dn2[1][i];
        d02[i] = dn2[0][i]-dn2[2][i];
      }
      
      a11 = sprod(d01, d01);
      a22 = sprod(d02, d02);
      a12 = sprod(d01, d02);
      a01 = sprod(dn2[0], d01);
      a02 = sprod(dn2[0], d02);
      hw += 0.5*sprod(dn2[0], dn2[0]);
      
      hr = a12*a12 - a11*a22;
      c1 = (a02*a12 - a01*a22)/hr;
      c2 = (a01*a12 - a02*a11)/hr;

      for (i=0; i<N-3*Rmax; i++) {
        rho2[i] = (1.-c1-c2)*rhon2[0][i] + c1*rhon2[1][i] + c2*rhon2[2][i];
      }
      
      c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0); 
      #pragma parallel for
      for (i = 0; i < N-3*Rmax; ++i) rho2[i] *= (double)Nn/c;

      magn = get_magnitude(n3);
 
      rl = get_radius(rho1, 1)*dr;
      rn = get_radius(rho2, 1)*dr;

      if (j%100 == 0) print_file(j);

      ++j;

      printf("\r%.2f %% ... %.2e ... %.2e ... %.2e ... %.4e ...  %.2f ... %.2f ... %.4f", 
	100.*j/250000, dphi_max, dphi2_max, (hw-hw_old), magn, 0.5*rl, 0.5*rn, pi_4*fr2int(rho2, dr, (int)(rn/dr), 0.)/Nn);
      fflush(stdout);
    
  }


  printf("\ndOm/drho=%.3e, dN=%.3e\n", dphi_max, dN);

  if (j >= ITS_picard) printf("Exceeded maximum number of iterations.\n");
  else printf("Did %d iterations\n", j);

}

void minimize_ng(double accuracy) {
  
  double hw = 0.0, hw2 = 0.0, hw_old = 1.0, hw_old2 = 1.0, hw_eps, c, dN, magn;
  double rho1_old[N], rho2_old[N], rhon1[3][N], rhon2[3][N], dn1[3][N], dn2[3][N], d01[N], d02[N];
  double c1, c2, a11, a22, a12, a01, a02, hr, dphi_max, dphi2_max;
  int i, j = 0, k, ITS = 250000;

  hw_eps = accuracy;

  printf("\nBeginning Ng-Iteration\n");

  while (((fabs(hw-hw_old) > hw_eps) || (fabs(hw2-hw_old2) > hw_eps)) && (j < ITS)) {

    hw_old = hw;
    hw_old2 = hw2;

    weighted_densities();

    derivative();

    dphi_max = get_magnitude(dphi);
    dphi2_max = get_magnitude(dphi2);

    a2 = fac*exp(-get_magnitude(dep2));

    if (j < 500) {
      #pragma omp parallel for private(i)
      for (i = 1; i < N-3*Rmax; i++) {
        
        rho1_old[i] = rho1[i];
        rho2_old[i] = rho2[i];

        rho1[i] = (1.-a1)*rho1[i]+a1*rho1l*exp(-dep[i]);
        rho2[i] = (1.-a2)*rho2[i]+a2*rho02l*exp(-dep2[i]);
        
      }
      
      c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0);
      #pragma parallel for reduction (+:hw)
      for (i = 0; i < N-3*Rmax; ++i) {
        rho2[i] *= (double)Nn/c;
        hw += 0.5*(rho1_old[i]-rho1[i])*(rho1_old[i]-rho1[i])+0.5*(rho2_old[i]-rho2[i])*(rho2_old[i]-rho2[i]);
      }

    } else {

      for (k=2; k>=0; k--) {

        weighted_densities();

        derivative();

        #pragma omp parallel for private(i)
        for (i=0; i<N-3*Rmax; i++) {
          rhon1[k][i] = (1.-a1)*rho1[i]+a1*rho1l*exp(-dep[i]);
          rhon2[k][i] = (1.-a2)*rho2[i]+a2*rho02l*exp(-dep2[i]);
        }
        
        c = pi_4*fr2int(rhon2[k], dr, N-3*Rmax, 0.0); 
        #pragma parallel for num_threads(4)
        for (i = 0; i < N-3*Rmax; ++i) rhon2[k][i] *= (double)Nn/c;

        #pragma omp parallel for private(i)
        for (i=0; i<N-3*Rmax; i++) {
          dn1[k][i] = rhon1[k][i] - rho1[i];
          dn2[k][i] = rhon2[k][i] - rho2[i];

          rho1[i] = rhon1[k][i];
          rho2[i] = rhon2[k][i];
        }
      }

        for (i=0; i<N-3*Rmax; i++) {
          d01[i] = dn1[0][i]-dn1[1][i];
          d02[i] = dn1[0][i]-dn1[2][i];
        }
        
        a11 = sprod(d01, d01);
        a22 = sprod(d02, d02);
        a12 = sprod(d01, d02);
        a01 = sprod(dn1[0], d01);
        a02 = sprod(dn1[0], d02);
        hw2 += 0.5*sprod(dn1[0], dn1[0]);

        hr = a12*a12 - a11*a22;
        c1 = (a02*a12 - a01*a22)/hr;
        c2 = (a01*a12 - a02*a11)/hr;

        for (i=0; i<N-3*Rmax; i++) {
          rho1[i] = (1.-c1-c2)*rhon1[0][i] + c1*rhon1[1][i] + c2*rhon1[2][i];
        }

  //=============================================================================
        for (i=0; i<N-3*Rmax; i++) {
          d01[i] = dn2[0][i]-dn2[1][i];
          d02[i] = dn2[0][i]-dn2[2][i];
        }
        
        a11 = sprod(d01, d01);
        a22 = sprod(d02, d02);
        a12 = sprod(d01, d02);
        a01 = sprod(dn2[0], d01);
        a02 = sprod(dn2[0], d02);
        hw2 += 0.5*sprod(dn2[0], dn2[0]);
        
        hr = a12*a12 - a11*a22;
        c1 = (a02*a12 - a01*a22)/hr;
        c2 = (a01*a12 - a02*a11)/hr;

        for (i=0; i<N-3*Rmax; i++) {
          rho2[i] = (1.-c1-c2)*rhon2[0][i] + c1*rhon2[1][i] + c2*rhon2[2][i];
        }
        
        c = pi_4*fr2int(rho2, dr, N-3*Rmax, 0.0); 
        #pragma parallel for
        for (i = 0; i < N-3*Rmax; ++i) rho2[i] *= (double)Nn/c;

    }


    dN = ((double)Nn-pi_4*fr2int(rho2, dr, N-3*Rmax, 0.))/Nn;
    magn = pi_4*get_magnitude(rho1)/3.;

    print_file(j);

    printf("\r%.2f %% ... %.2e ... %.2e ... %.2e ... %.4e ...  %.2f ... %.2f ... %.4f", 100.*j/ITS, dphi_max, dphi2_max,
     fabs(hw-hw_old)+fabs(hw2-hw_old2), magn, get_radius(rho1, 1)*dr/2., get_radius(rho2, 1)*dr/2., pi_4*fr2int(rho2, dr, get_radius(rho2, 1), 0.)/Nn);
    fflush(stdout);

    j++;
  }

  printf("\nELG<=%.3e, dN/N=%.3e\n", get_magnitude(dphi), dN);

  if (j >= ITS) printf("Exceeded maximum number of iterations.\n");
  else printf("Did %d iterations\n", j);
}

void minimize_next(double lmu, double lfac, double accuracy) {

  fac = lfac;
  mu1l = lmu-log(pi_4/3.*rho1l);
  minimize_ng(accuracy);
  print_results();

}



int main(int argc, char *argv[]) {
 
  double Rb_r, Rn_r;

  if (argc >= 4) {

    Rb_r = 2.*atof(argv[1]);
    a1 = atof(argv[2]);
    fac = atof(argv[3]);
    
    Rb = (int)(Rb_r * R);
    Rn = Rb;

    data_init();
    
    if (argc >= 5) mu1l = atof(argv[4])-log(pi_4/3.*rho1l);
     
    //minimize_ng();
    minimize_eff(Rb_r, 10.);
    print_results();

    //minimize_next(-5.7648, fac, 1.e-16);
    //minimize_next(-5.6176, fac, 1.e-16);
    //minimize_next(-5.5521, fac, 1.e-16);
    //minimize_next(-5.4902, fac, 1.e-16);
  
  } else {
    printf("\033[0;31m");
    printf("Enter %s Rb a1 fac (mu_ex)\n", argv[0]);
  }

  return 0;
}
