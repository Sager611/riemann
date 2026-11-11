
// edesmo.c: exact optimal c2 = sup E - inf E for Q_q (no p^2 | n for p<=q), c1 = prod(1-1/p^2), E(t)=Q_q(t)-c1 t.
#include <stdio.h>
#include <stdint.h>
int main(){ int ps[6]={2,3,5,7,11,13};
  for(int k=1;k<=6;k++){ int64_t R=1,A=1; for(int i=0;i<k;i++){R*=ps[i]*ps[i]; A*=ps[i]*ps[i]-1;}
    int64_t Q=0, emax=0, emin=0; // R*E values
    for(int64_t n=1;n<=R;n++){ int ok=1; for(int i=0;i<k;i++){ int64_t p2=ps[i]*ps[i]; if(n%p2==0){ok=0;break;} }
      if(ok){ int64_t before=R*Q-A*n; if(before<emin)emin=before; Q++; int64_t at=R*Q-A*n; if(at>emax)emax=at; } }
    int64_t num=emax-emin; // c2 = num/R ; reduce
    int64_t a=num,b=R; while(b){int64_t t=a%b;a=b;b=t;} 
    int64_t g1=A,g2=R; int64_t x=g1,y=g2; while(y){int64_t t=x%y;x=y;y=t;}
    printf("q=%d: c1=%lld/%lld  sup E=%lld/%lld inf E=%lld/%lld  c2=sup-inf=%lld/%lld = %.6f\n",ps[k-1],(long long)(A/x),(long long)(R/x),(long long)emax,(long long)R,(long long)emin,(long long)R,(long long)(num/a),(long long)(R/a),(double)num/R);
  } return 0; }
