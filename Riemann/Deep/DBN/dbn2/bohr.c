/* bohr.c — heuristic Bohr-torus obstruction for the Polymath zero-free
   criterion at window N: minimise Phi(theta) = |sum w_n z_n| - g |sum u_n z_n|
   over completely multiplicative unimodular z_n = exp(-i theta(n)),
   w_n = bt(n) n^{-sigma}, u_n = n^y bt(n) n^{-sigma}, g = (x_N/4pi)^{-y/2}.
   min Phi <= 0  ==>  no phase-uniform lower bound (triangle/mollifier type)
   can certify H_t != 0 on the whole window (Kronecker heuristic).
   usage: bohr N t y2 iters seed init(0=Liouville,1=random,2=mixed)
   Exploration only. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;
int main(int argc,char**argv){
  u64 N=strtoull(argv[1],0,10); double t=atof(argv[2]), y2=atof(argv[3]); int iters=atoi(argv[4]); unsigned seed=atoi(argv[5]); int init=atoi(argv[6]); const char *loadf = argc>7? argv[7] : 0; const char *savef = argc>8? argv[8] : 0;
  double y=sqrt(y2), pi=M_PI;
  double xN=4*pi*(double)N*(double)N - pi*t/4, L=log(xN/(4*pi));
  double sig=(1+y)/2 + t/4*L; double g=exp(-y/2*L);
  unsigned *spf=calloc(N+1,sizeof(unsigned)); int *pidx=malloc((N+1)*sizeof(int));
  int np=0; for(u64 i=2;i<=N;i++){ if(!spf[i]){ for(u64 j=i;j<=N;j+=i) if(!spf[j]) spf[j]=(unsigned)i; pidx[i]=np++; } }
  double *w=malloc((N+1)*sizeof(double)), *u=malloc((N+1)*sizeof(double)), *th=malloc((N+1)*sizeof(double)), *cn=malloc((N+1)*sizeof(double));
  double *tp=malloc(np*sizeof(double)), *gr=malloc(np*sizeof(double)), *m1=calloc(np,sizeof(double)), *m2=calloc(np,sizeof(double));
  for(u64 n=1;n<=N;n++){ double l=log((double)n); double b=exp(0.25*t*l*l); w[n]=b*exp(-sig*l); u[n]=w[n]*exp(y*l); }
  srand(seed);
  { int k=0; for(u64 i=2;i<=N;i++) if(spf[i]==i){ double r=2*pi*rand()/(double)RAND_MAX; tp[k]= init==0? pi : init==1? r : (i<100? pi : r); k++; } }
  if(loadf && loadf[0]!='-'){ FILE *fp=fopen(loadf,"rb"); if(fp){ size_t rd=fread(tp,sizeof(double),np,fp); fclose(fp); fprintf(stderr,"loaded %zu\n",rd);} }
  double best=1e9, lr=0.05; 
  for(int it=0; it<=iters; it++){
    th[1]=0; for(u64 n=2;n<=N;n++){ unsigned p=spf[n]; th[n]=th[n/p]+tp[pidx[p]]; }
    double Fr=0,Fi=0,Gr=0,Gi=0;
    for(u64 n=1;n<=N;n++){ double c=cos(th[n]), s=sin(th[n]); Fr+=w[n]*c; Fi-=w[n]*s; Gr+=u[n]*c; Gi-=u[n]*s; }
    double aF=hypot(Fr,Fi), aG=hypot(Gr,Gi); double Phi=aF-g*aG;
    if(Phi<best) best=Phi;
    if(it%200==0 || it==iters) fprintf(stderr,"it %d Phi %.6f |F| %.6f g|G| %.6f\n",it,Phi,aF,g*aG);
    if(it==iters) break;
    /* dPhi/dtheta_p = sum_n v_p(n) c_n ; d z_n/d theta = -i z_n ; z_n=cos-i sin */
    for(u64 n=1;n<=N;n++){ double c=cos(th[n]), s=sin(th[n]);
      /* dF = -i w z = -i w (c - i s) = -w s - i w c ; d|F| = (Fr*dFr + Fi*dFi)/|F| */
      double dF = (Fr*(-w[n]*s) + Fi*(-w[n]*c))/aF;
      double dG = (Gr*(-u[n]*s) + Gi*(-u[n]*c))/aG;
      cn[n]= dF - g*dG; }
    for(int k=0;k<np;k++) gr[k]=0;
    for(u64 n=2;n<=N;n++){ u64 m=n; while(m>1){ unsigned p=spf[m]; gr[pidx[p]]+=cn[n]; m/=p; } }
    /* Adam */
    double b1=0.9,b2=0.999,eps=1e-12; double bc1=1-pow(b1,it+1), bc2=1-pow(b2,it+1);
    for(int k=0;k<np;k++){ m1[k]=b1*m1[k]+(1-b1)*gr[k]; m2[k]=b2*m2[k]+(1-b2)*gr[k]*gr[k];
      tp[k]-= lr*(m1[k]/bc1)/(sqrt(m2[k]/bc2)+eps); }
    if(it%300==299) lr*=0.6;
  }
  if(savef){ FILE *fp=fopen(savef,"wb"); fwrite(tp,sizeof(double),np,fp); fclose(fp); }
  printf("N=%llu t=%.6f y2=%.6f sigma=%.6f g=%.6f init=%d seed=%u bestPhi=%.6f Lambda=%.6f\n",N,t,y2,sig,g,init,seed,best,t+y2/2);
  return 0;
}
