/* prolate.c — negative spectrum of the Connes–Moscovici self-adjoint prolate operator W_sa.
   Problem on (lam, inf):  ((x^2-lam^2) u')' + (om^2 x^2 + E^2) u = 0,  om = 2 pi lam,  mu = -E^2,
   u regular (analytic) at x = lam (u(lam)=1).  Even eigenvalue <=> lim [u,s] = 0, s = sin(om x)/x;
   odd <=> lim [u,c] = 0, c = cos(om x)/x  (CM21 l.567-583).
   Method: Frobenius series at lam, then Taylor-series steps in x (polynomial coefficients, recurrence),
   Pruefer angle phi = atan2(k u, sqrt(x^2-lam^2) u'), k = sqrt(E^2+om^2 x^2), unwrapped;
   Theta(E) = lim (phi - om x) evaluated at x = X plus tail: int (k - kap sinh t) dt  +  [a cos2phi/(2 phidot)] - int a^2/(2k) dt
   (t = acosh(x/lam), kap = 2 pi lam^2, a = (1/2)coth t + kap^2 sinh t cosh t /(2k^2)), error O(1/k(X)^2).
   Usage:
     prolate lam theta E [omX]
     prolate lam grid E0 E1 dE [omX]
     prolate lam eig parity Emax [omX] [dEscan] [Emin]
     prolate lam wr E0 E1 dE X        (raw Wronskians [u,s](X), [u,c](X), no tail: independent sign-change count)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static double LAM, OM, KAP, E2v, Ev;
static double STEPPH=1.5, RADF=0.3; /* phase per Taylor step, fraction of distance to x=lam (env PF_STEP, PF_RAD) */
static long NSTEPS;

/* Gauss-Legendre nodes on [-1,1] */
#define NGL 10
static double glx[NGL], glw[NGL];
static void gl_init(void){
  int n=NGL;
  for(int i=0;i<n;i++){
    double x=cos(M_PI*(i+0.75)/(n+0.5)), dp=0;
    for(int it=0;it<100;it++){
      double p0=1,p1=x;
      for(int k=2;k<=n;k++){double p2=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p2;}
      dp=n*(x*p1-p0)/(x*x-1);
      double dx=p1/dp; x-=dx; if(fabs(dx)<1e-16)break;
    }
    glx[i]=x; glw[i]=2/((1-x*x)*dp*dp);
  }
}

/* Frobenius at lam: u = sum U_n h^n */
static void frob(double h,double *u,double *up){
  double V0=OM*OM*LAM*LAM+E2v;
  double Um2=0,Um1=0,Un=1,s=1,sp=0,hn=1;
  for(int n=0;n<300;n++){
    double Unp1=-((n*(n+1.0)+V0)*Un+2*OM*OM*LAM*Um1+OM*OM*Um2)/(2*LAM*(n+1.0)*(n+1.0));
    double t=Unp1*hn*h, tp=(n+1)*Unp1*hn;
    s+=t; sp+=tp;
    Um2=Um1;Um1=Un;Un=Unp1;hn*=h;
    if(n>20 && fabs(t)<1e-19*fabs(s) && fabs(tp*h)<1e-19*fabs(s)) break;
  }
  *u=s;*up=sp;
}

/* one Taylor step from xc with (u,up) to xc+h */
static void tstep(double xc,double h,double *u,double *up){
  double pc=xc*xc-LAM*LAM, Vc=OM*OM*xc*xc+E2v, w2=OM*OM;
  double c[400]; c[0]=*u; c[1]=*up;
  double su=c[0]+c[1]*h, sp=c[1], hn=h; /* hn = h^(n) for term n */
  double A=fabs(c[0])+fabs(c[1])*h+1e-300;
  int small=0;
  for(int n=0;n<396;n++){
    double cm1=n>=1?c[n-1]:0, cm2=n>=2?c[n-2]:0;
    c[n+2]=-(2*xc*(n+1.0)*(n+1.0)*c[n+1]+(n*(n+1.0)+Vc)*c[n]+2*w2*xc*cm1+w2*cm2)/(pc*(n+1.0)*(n+2.0));
    double hn2=hn*h; /* h^(n+2) */
    double t=c[n+2]*hn2, tp=(n+2)*c[n+2]*hn;
    su+=t; sp+=tp; hn=hn2;
    if(fabs(t)<1e-18*A && fabs(tp*h)<1e-18*A){ if(++small>=3) break; } else small=0;
  }
  *u=su;*up=sp;
}

typedef struct { double theta, raw, tail1, tail2, phiX, X; long steps; double ws, wc; } Res;

static double tail_int(double tX,int which){
  /* which=1: int (E^2+kap^2)/(k+kap sinh t) ; which=2: int a^2/(2k) ; over [tX, inf) via panels in s=t-tX */
  double tot=0;
  double edges[]={0,0.25,0.5,1,1.5,2,3,4,5,6,8,10,13,16,20,25,30,40,50,60};
  int ne=sizeof(edges)/sizeof(edges[0]);
  for(int p=0;p<ne-1;p++){
    double a0=edges[p],b0=edges[p+1],m=(a0+b0)/2,r=(b0-a0)/2;
    for(int i=0;i<NGL;i++){
      double t=tX+m+r*glx[i];
      double ch=cosh(t), sh=sinh(t);
      double k=sqrt(E2v+KAP*KAP*ch*ch);
      double v;
      if(which==1) v=(E2v+KAP*KAP)/(k+KAP*sh);
      else { double a=0.5*ch/sh+KAP*KAP*sh*ch/(2*k*k); v=a*a/(2*k); }
      tot+=r*glw[i]*v;
    }
  }
  return tot;
}

static double wrap(double d){ while(d>M_PI) d-=2*M_PI; while(d<=-M_PI) d+=2*M_PI; return d; }

/* integrate to X; returns Theta and raw Wronskians at X */
static Res theta_E(double E,double omX){
  Res R; Ev=E; E2v=E*E;
  double X=omX/OM; if(X<4*LAM) X=4*LAM;
  double V0=OM*OM*LAM*LAM+E2v;
  double h1=fmin(0.25*LAM,0.05*2*LAM/(1+V0));
  double u,up; frob(h1,&u,&up);
  double x=LAM+h1;
  double k=sqrt(E2v+OM*OM*x*x), sq=sqrt(x*x-LAM*LAM);
  double al=atan2(k*u,sq*up); /* in (0,pi) near pi/2 */
  long W=0; /* phi = 2 pi W + al */
  long steps=0;
  while(x<X){
    double pc=x*x-LAM*LAM;
    double freq=sqrt((E2v+OM*OM*x*x)/pc);
    double h=fmin(RADF*(x-LAM),STEPPH/freq);
    if(x+h>X) h=X-x;
    { volatile double xn=x+h; h=xn-x; } /* exact step: x+h representable, no drift */
    tstep(x,h,&u,&up); x+=h; steps++;
    /* renormalise amplitude to avoid under/overflow */
    double nrm=fabs(u)+fabs(up); if(nrm>1e100||nrm<1e-100){u/=nrm;up/=nrm;}
    k=sqrt(E2v+OM*OM*x*x); sq=sqrt(x*x-LAM*LAM);
    double aln=atan2(k*u,sq*up);
    double d=wrap(aln-al);
    double phinew=al+d; /* may leave (-pi,pi] */
    if(phinew>M_PI){W++; } else if(phinew<=-M_PI){W--;}
    al=aln;
  }
  /* phi(X) - om X */
  double lx=LAM*x; double fr=lx-floor(lx); /* om x/(2pi) = lam x */
  double raw=2*M_PI*((double)W-floor(lx))+al-2*M_PI*fr;
  double tX=acosh(x/LAM);
  double ch=cosh(tX), sh=sinh(tX);
  k=sqrt(E2v+KAP*KAP*ch*ch);
  double a=0.5*ch/sh+KAP*KAP*sh*ch/(2*k*k);
  double phi=al; double phid=k+a*sin(2*phi);
  double t1=tail_int(tX,1);
  double t2=a*cos(2*phi)/(2*phid)-tail_int(tX,2);
  R.theta=raw+t1+t2; R.raw=raw; R.tail1=t1; R.tail2=t2; R.X=x; R.steps=steps; R.phiX=al;
  /* raw Wronskians [u,s](X), [u,c](X), [f,g] = (x^2-lam^2)(f g' - f' g); normalised by amplitude rho */
  double s=sin(OM*x)/x, sp=OM*cos(OM*x)/x-sin(OM*x)/(x*x);
  double cc=cos(OM*x)/x, cp=-OM*sin(OM*x)/x-cos(OM*x)/(x*x);
  double pc=x*x-LAM*LAM;
  double m2=sq*k; /* m^2 = sqrt(x^2-lam^2) k */
  double rho=sqrt(m2*u*u+pc*pc*up*up/m2);
  R.ws=pc*(u*sp-up*s)/rho; R.wc=pc*(u*cp-up*cc)/rho;
  NSTEPS=steps;
  return R;
}

static double omX_default(double E){ double v=30*E; return v>1e5?v:1e5; }

int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage\n");return 1;}
  gl_init();
  if(getenv("PF_STEP")) STEPPH=atof(getenv("PF_STEP"));
  if(getenv("PF_RAD")) RADF=atof(getenv("PF_RAD"));
  LAM=atof(argv[1]); OM=2*M_PI*LAM; KAP=2*M_PI*LAM*LAM;
  const char*mode=argv[2];
  if(!strcmp(mode,"theta")){
    double E=atof(argv[3]); double omX=argc>4?atof(argv[4]):omX_default(E);
    Res R=theta_E(E,omX);
    double Phi=E>0?E*log(4*E/KAP)-E+M_PI/4:NAN;
    printf("E=%.15g Theta=%.15f raw=%.15f tail1=%.15g tail2=%.6e X=%.6g steps=%ld Theta-Phi=%.6e ws=%.6e wc=%.6e\n",E,R.theta,R.raw,R.tail1,R.tail2,R.X,R.steps,R.theta-Phi,R.ws,R.wc);
  } else if(!strcmp(mode,"grid")){
    double E0=atof(argv[3]),E1=atof(argv[4]),dE=atof(argv[5]); double omX=argc>6?atof(argv[6]):-1;
    for(double E=E0;E<=E1+1e-12;E+=dE){
      Res R=theta_E(E,omX>0?omX:omX_default(E));
      printf("%.15g %.15f\n",E,R.theta); fflush(stdout);
    }
  } else if(!strcmp(mode,"loggrid")){
    double E0=atof(argv[3]),E1=atof(argv[4]); int n=atoi(argv[5]); double omX=argc>6?atof(argv[6]):-1;
    for(int i=0;i<n;i++){
      double E=E0*pow(E1/E0,(double)i/(n-1));
      Res R=theta_E(E,omX>0?omX:omX_default(E));
      printf("%.15g %.15f %ld\n",E,R.theta,R.steps); fflush(stdout);
    }
  } else if(!strcmp(mode,"eig")){
    int par=atoi(argv[3]); double Emax=atof(argv[4]); double omX=argc>5?atof(argv[5]):-1; double dEs=argc>6?atof(argv[6]):0.2; double Emin=argc>7?atof(argv[7]):0;
    double off=par?0.5*M_PI:0.0;
    double Ea=Emin; Res Ra=theta_E(Emin>0?Emin:1e-9,omX>0?omX:omX_default(Emin)); double Ta=Ra.theta;
    printf("# lam=%.15g parity=%d Emin=%.15g Theta(Emin)=%.15f\n",LAM,par,Emin,Ta);
    for(double Eb=Emin+dEs;Eb<=Emax+1e-12;Eb+=dEs){
      Res Rb=theta_E(Eb,omX>0?omX:omX_default(Eb)); double Tb=Rb.theta;
      if(Tb<=Ta){printf("# NONMONOTONE at E=%.6f %.15f %.15f\n",Eb,Ta,Tb);}
      /* levels j*pi+off in (Ta,Tb] */
      long j0=(long)ceil((Ta-off)/M_PI-1e-15), j1=(long)floor((Tb-off)/M_PI);
      for(long j=j0;j<=j1;j++){
        double target=j*M_PI+off; if(target<=Ta) continue;
        double lo=Ea,hi=Eb,flo=Ta-target,fhi=Tb-target;
        double E=0;int side=0;
        for(int it=0;it<80;it++){
          E=(lo*fhi-hi*flo)/(fhi-flo);
          if(!(E>lo&&E<hi)) E=0.5*(lo+hi);
          Res Rm=theta_E(E,omX>0?omX:omX_default(E)); double fm=Rm.theta-target;
          if(fabs(fm)<1e-13||hi-lo<1e-14*hi) {lo=hi=E;break;}
          if((fm<0)){lo=E;flo=fm; if(side==-1) fhi*=0.5; side=-1;} else {hi=E;fhi=fm; if(side==1) flo*=0.5; side=1;}
        }
        Res Rf=theta_E(E,omX>0?omX:omX_default(E));
        printf("%ld %.15f %.3e %.6e %.6e\n",j,E,Rf.theta-target,Rf.ws,Rf.wc); fflush(stdout);
      }
      Ea=Eb;Ta=Tb;
    }
  } else if(!strcmp(mode,"wr")){
    double E0=atof(argv[3]),E1=atof(argv[4]),dE=atof(argv[5]),X=atof(argv[6]);
    for(double E=E0;E<=E1+1e-12;E+=dE){
      Res R=theta_E(E,OM*X);
      printf("%.10g %.12e %.12e\n",E,R.ws,R.wc); fflush(stdout);
    }
  }
  return 0;
}
