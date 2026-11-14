
# Falsification of both bounds against exact S(T), N(T) - M(T) on [1, 74920] from Odlyzko's first 1e5 zeros.
import math, json
import numpy as np
import mpmath as mp
g = np.loadtxt('.scratch/stcaudit/zeros1')
n = np.arange(1, len(g)+1, dtype=float)
def theta_as(t):
    return t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3) + 31/(80640*t**5)
# check asymptotic theta against mpmath at a few points
for t in [14.134725142, 100.0, 5000.0]:
    assert abs(theta_as(np.array([t]))[0] - float(mp.siegeltheta(t))) < 1e-9
th = theta_as(g)
Sp = n - th/np.pi - 1          # S just after gamma_n (simple zeros)
Sm = Sp - 1                    # S just before gamma_n
M = g/(2*np.pi)*np.log(g/(2*np.pi*np.e))
NMp = n - M; NMm = n - 1 - M
# [1, gamma_1): N = 0
tt = np.linspace(1, g[0]-1e-9, 2000)
S0 = np.array([-float(mp.siegeltheta(t))/math.pi - 1 for t in tt]); NM0 = -tt/(2*np.pi)*np.log(tt/(2*np.pi*np.e))
CG_S = [(0.20176, 0.844), (0.10823, 1.286), (0.09473, 1.946), (0.09349, 2.099), (0.08399, 4.318), (0.08382, 4.444)]
CG_N = [(0.20176, 1.369), (0.10823, 2.171), (0.09473, 2.834), (0.09349, 2.987), (0.08399, 5.202), (0.08382, 5.326)]
IS = [(0.140750, 0.150951, 1.377296), (0.134791, 0.159606, 1.505817), (0.129350, 0.169319, 1.643495), (0.122384, 0.186030, 1.859015),
      (0.118798, 0.197172, 1.996010), (0.115144, 0.211179, 2.163574), (0.111259, 0.154570, 2.730936), (0.108109, 0.207558, 2.795522), (0.105234, 0.214408, 3.175432)]
def env_cg(T, lines): L = np.log(T); return np.min([a*L + b for a, b in lines], axis=0)
def env_is(T, add=0.0): L = np.log(T); return np.min([a*L + b*np.log(L) + c + add for a, b, c in IS], axis=0)
absS = np.maximum(np.abs(Sp), np.abs(Sm)); absNM = np.maximum(np.abs(NMp), np.abs(NMm))
res = {}
res['max|S| on [gamma_1, 74920]'] = float(absS.max()); res['argmax T'] = float(g[absS.argmax()])
res['max|S| on [1, gamma_1)'] = float(np.abs(S0).max())
res['CG |S| min slack'] = float(min((env_cg(g, CG_S) - absS).min(), (env_cg(tt, CG_S) - np.abs(S0)).min()))
res['CG |N-M| min slack'] = float(min((env_cg(g, CG_N) - absNM).min(), (env_cg(tt, CG_N) - np.abs(NM0)).min()))
# each CG line separately (Theorem says each line holds for all T >= 1)
res['CG each line S min slack'] = [float(min((a*np.log(g)+b - absS).min(), (a*np.log(tt)+b - np.abs(S0)).min())) for a, b in CG_S]
res['CG each line N min slack'] = [float(min((a*np.log(g)+b - absNM).min(), (a*np.log(tt)+b - np.abs(NM0)).min())) for a, b in CG_N]
m = g >= math.e; mt = tt >= math.e
res['IS each row S min slack (T>=e)'] = [float(min((a*np.log(g)+b*np.log(np.log(g))+c - absS).min(), (a*np.log(tt[mt])+b*np.log(np.log(tt[mt]))+c - np.abs(S0[mt])).min())) for a, b, c in IS]
res['IS Cor1.4 |N-M| min slack (row 2, T>=e)'] = float(min((0.134791*np.log(g)+0.159606*np.log(np.log(g))+2.380817 - absNM).min(), (0.134791*np.log(tt[mt])+0.159606*np.log(np.log(tt[mt]))+2.380817 - np.abs(NM0[mt])).min()))
# computed-bound inputs used by Ismail below T0: |S|<=1 on [e,280], <=2 on (280, 6.8e6]
res['max|S| on [e,280]'] = float(max(np.abs(S0[mt]).max(), absS[g <= 280].max(), np.abs(Sp[g <= 280]).max()))
res['max|S| on (280, 74920]'] = float(absS[g > 280].max())
print(json.dumps(res, indent=1))
json.dump(res, open('Riemann/Deep/STCAudit/falsify_S.json', 'w'), indent=1)
