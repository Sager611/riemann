# mixed.py (explore23): DH on the critical line as the sum of two uncoupled 'cavities': Lambda_DH(1/2+it) = u (Z(t) + Z(-t)),
# Z(t) = eps^{1/2} Lambda(1/2+it, chibar) real (chi mod 5, chi(2) = i). Zeros of Z on t>0: L(s,chibar); on t<0 (mirrored): L(s,chi).
# Test: do DH's off-line ordinates sit where the two zero sets fail to interlace (Hermite-Kakeya-Obreschkoff)?
import sys, json, numpy as np, mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dh import hurwitz, vartheta, Zdh
PI = np.pi
chib = {1: 1.0, 2: -1j, 3: 1j, 4: -1.0}   # conj of chi(2) = i
TM = 420.0
ts = np.arange(-TM, TM, 0.01)
s = 0.5 + 1j*ts
L = 5.0**(-s)*sum(c*hurwitz(s, r/5.0, int(TM/2) + 60) for r, c in chib.items())
vt = np.sign(ts)*vartheta(np.abs(ts))                 # vartheta is odd in t
W = np.exp(1j*vt)*L                                    # = (real Z) * e^{-i phi0}
phi0 = np.angle(W[np.argmax(np.abs(W))])
Z = (W*np.exp(-1j*phi0)).real
out = {'imag_rms_after_phase': float(np.sqrt(np.mean((W*np.exp(-1j*phi0)).imag**2))/np.sqrt(np.mean(Z**2)))}
# DH = u (Z(t) + Z(-t)): compare with Zdh on t > 0 up to sign/scale
tp = ts[ts > 1]; Zp = Z[ts > 1]; Zm = Z[::-1][ts[::-1] < -1][::-1] if False else np.interp(-tp, ts, Z)
summ = Zp + Zm
zd = Zdh(tp[::50])
ratio = zd/summ[::50]
out['DH_over_(Z(t)+Z(-t))_ratio_[median, rel spread]'] = [float(np.median(ratio)), float(np.std(ratio[np.abs(summ[::50]) > 0.05])/abs(np.median(ratio)))]
def zeros_of(y, x):
    i = np.nonzero(np.sign(y[:-1])*np.sign(y[1:]) < 0)[0]
    return x[i] - y[i]*(x[i + 1] - x[i])/(y[i + 1] - y[i])
zb = zeros_of(Zp, tp)            # L(s, chibar) zeros at heights > 0
za = zeros_of(Zm, tp)            # L(s, chi) zeros at heights > 0 (from Z(-t))
D = json.load(open('dh_zeros.json')); off = np.array(D['offline'])
# merged sequence labels; a non-interlacing 'run' = two consecutive zeros (merged order) from the same set
lab = np.concatenate([np.zeros(len(za)), np.ones(len(zb))]); allz = np.concatenate([za, zb]); o = np.argsort(allz)
allz = allz[o]; lab = lab[o]
runs = [(allz[i], allz[i + 1]) for i in range(len(allz) - 1) if lab[i] == lab[i + 1]]
out['zeros_chi_chibar'] = [int(len(za)), int(len(zb))]
out['non_interlacing_runs'] = int(len(runs))
hit = []
for b_, g_ in off:
    inrun = any(a - 1e-9 <= g_ <= c + 1e-9 for (a, c) in runs)
    # nearest run distance in units of the local spacing
    dist = min(min(abs(g_ - a), abs(g_ - c)) if not (a <= g_ <= c) else 0.0 for (a, c) in runs)
    hit.append([round(float(g_), 3), round(float(b_), 4), bool(inrun), round(float(dist), 3)])
out['offline_in_run_[gamma, beta, inside, dist]'] = hit
# how many runs contain a DH off-line ordinate; local test: sign of Z1+Z2 at run ends
out['runs_with_offline'] = int(sum(any(a <= g_ <= c for g_ in off[:, 1]) for (a, c) in runs))
json.dump(out, open('mixed.json', 'w'), indent=1)
print(json.dumps(out))
