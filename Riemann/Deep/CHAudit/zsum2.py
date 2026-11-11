
import numpy as np
from flint import acb, acb_series, ctx
S='/Users/adriansager/playground/riemann/.scratch/chaudit/'
d=np.load(S+'zp.npy'); g0=d[0]; z0=d[1]
a=np.load(S+'z2a.npy'); b=np.load(S+'z2b.npy')
g1=np.concatenate([a[0],b[0]]); z1=np.concatenate([a[1],b[1]])
print('new zeros',len(g1),'sorted',bool(np.all(np.diff(g1)>0)),'min gap %.3e'%np.min(np.diff(g1)))
T=1893193.5
n_tot=len(g0)+np.sum(g1<T)
print('N(T=1893193.5) found = %d ; Rosser-Yohe-Schoenfeld (as quoted in 2511.14736): 3500000'%n_tot)
ctx.prec=96
for nn in [2001053,2500000,3000000,3499999,3500000]:
    k=nn-len(g0)-1
    zf=acb.zeta_zero(nn); gf=float(zf.imag.mid()); zpf=abs(complex(acb_series([acb(0.5,gf),1],prec=2).zeta()[1]))
    print('zero #%d: ours %.9f flint %.9f | |zeta prime| ours %.8f flint %.8f'%(nn,g1[k],gf,z1[k],zpf))
g=np.concatenate([g0,g1]); zp=np.concatenate([z0,z1]); m=g<T; dl=np.pi/(2*T)
Sum=dl*np.sum(np.abs(1/np.tanh(dl*(0.5+1j*g[m])))/zp[m])
print('PLATT T=1893193.5: sum over %d zeros = %.12f ; paper 2.85417779533422 ; diff %.2e'%(m.sum(),Sum,Sum-2.85417779533422))
