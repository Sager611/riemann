
from mpmath import mp, mpf, log, exp, findroot
mp.dps=30
def rh(L):
    ll=log(L); return (mpf('0.5')+mpf('6.992')/(L**mpf('0.252')*ll))*L/ll
def mk(c): return lambda L: log(mpf(c))+L/6+log(L)
def fodge(L): t6=exp(L/6); return log(mpf('0.3080532')*t6*L+mpf('2.405897')*t6+mpf('506.6'))
for name,f in [('ours0.568',mk('0.568')),('ours0.4821',mk('0.4821')),('fodge',fodge)]:
    r=findroot(lambda L: rh(L)-f(L), 110)
    print(name,'crossover log t =',mp.nstr(r,8),' log10 t =',mp.nstr(r/log(10),6))
# check sign after crossover (RH smaller) at a few points
for L10 in [45,47,49,51,60,100]:
    L=L10*log(10); print(L10, 'rh-fodge', mp.nstr(rh(L)-fodge(L),5), 'rh-0.4821', mp.nstr(rh(L)-mk('0.4821')(L),5))
