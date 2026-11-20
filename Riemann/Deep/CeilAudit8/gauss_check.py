import cmath, math, random
def G(a,k,q): return sum(cmath.exp(2j*math.pi*((a*r*r+k*r)%q)/q) for r in range(q))
def e(x): return cmath.exp(2j*math.pi*x)
random.seed(1); bad=0; tested=0; vanish_ok=0
for q in [6,10,14,18,22,30,8,12,16,20,24,28,36,40,44,9,15,21]:
  for a in range(1,3*q):
    if math.gcd(a,q)!=1: continue
    for k in range(0,2*q):
      g=G(a,k,q); par=k%2
      if q%2==0:
        should0 = (par==0) if q%4==2 else (par==1)
        if should0: vanish_ok += abs(g)<1e-8; tested+=1; bad += abs(g)>1e-8; continue
    # phase test in h' basis (even q) / h basis (odd q)
    for b in range(1,2*q):
      if q%4==2 and b%2==0: continue
      if q%4==0 and b%2==1: continue
      for h in range(0,6):
        if q%2:
          ia=pow(4*a,-1,q); ia2=pow(2*a,-1,q); pred=e((ia2*b*h-ia*h*h)/q); g1=G(a,b-h,q)/G(a,b,q)
        elif q%4==2:
          q1=q//2; i2=pow(2*a,-1,q1); pred=e((i2*b*h-i2*h*h)/q1); g1=G(a,b-2*h,q)/G(a,b,q)
        else:
          ia=pow(a,-1,q); pred=e((ia*b*h-ia*h*h)/q); g1=G(a,b-2*h,q)/G(a,b,q)
        tested+=1; bad += abs(g1-pred)>1e-7
  # |G| on allowed class
  a=1 if math.gcd(1,q)==1 else 1
print("tested",tested,"bad",bad,"vanish_ok",vanish_ok)
for q in [6,10,8,12,16]:
  a=1; k=1 if q%4==2 else 0
  print(q, round(abs(G(a,k,q)),6), round(math.sqrt(2*q),6))
