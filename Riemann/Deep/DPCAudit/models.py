
import json, numpy as np, mpmath as mp
d=json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/dh_zeros.json'))
off=d['offline']; on=d['online']
print("DH: TMAX",d['TMAX'],"online",len(on),"offline",len(off), "sample offline:", off[:3])
# pairing check: functional equation pairs beta+i gamma with 1-beta+i gamma
gam=sorted(float(z[1]) if isinstance(z,(list,tuple)) else float(z) for z in on)
print("online gaps normalized: min %.3f"%min((g2-g1)*np.log(g1/(2*np.pi))/(2*np.pi) for g1,g2 in zip(gam,gam[1:]) if g1>20))
# picket-fence model: Lambda(s)=gamma(s)+gamma(1-s); on line = 2 Re gamma(1/2+it) = 2|gamma| cos(theta(t)); zeros theta = pi/2 mod pi
# V=(1/pi)Log(Z+i0), Z=cos(th+a/2)/cos(th-a/2), N = #zeros of cos in (th-a/2, th+a/2); th uniform mod pi
for phi in [0.3,0.6,0.9]:
    a=np.pi*phi
    th=(np.arange(2000000)+0.5)/2000000*np.pi
    Z=np.cos(th+a/2)/np.cos(th-a/2)
    # count k with pi/2+k pi in (th-a/2,th+a/2)
    N=np.floor((th+a/2-np.pi/2)/np.pi)-np.floor((th-a/2-np.pi/2)/np.pi)
    X=np.log(np.abs(Z))/np.pi
    V=X+1j*N
    m=[np.mean(V**j) for j in range(1,6)]
    print("picket phi=%.1f: P(N=1)=%.5f  E[V^j] j=1..5:"%(phi,np.mean(N)), np.round(m,4), " (i phi)^j:", np.round([(1j*phi)**j for j in range(1,6)],4))
    # sign check: sign Z = (-1)^N
    print("   sign(Z)==(-1)^N:", bool(np.all(np.sign(Z)==(-1.0)**N)))
