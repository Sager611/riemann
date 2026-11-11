import mpmath as mp
mp.mp.dps = 20
exec(open('t1_phase.py').read().split('ts = [')[0])
ts=[mp.mpf(7003.0)+k*mp.mpf('0.02') for k in range(301)]
v=[phidot(t) for t in ts]
segs=[];cur=None
for t,x in zip(ts,v):
    if x<0 and cur is None: cur=[t,t]
    elif x<0: cur[1]=t
    elif cur is not None: segs.append(cur); cur=None
if cur: segs.append(cur)
print('negative phase-velocity windows in [7003,7009]:', ', '.join('[%s,%s]'%(mp.nstr(a,7),mp.nstr(b,7)) for a,b in segs))
print('min %s' % mp.nstr(min(v),5))
zs=[mp.im(mp.zetazero(n)) for n in range(6706,6714)]
print('zeros 6706..6713:', ' '.join(mp.nstr(z,9) for z in zs))
print('phidot at zeros:', ' '.join(mp.nstr(phidot(z),3) for z in zs))
