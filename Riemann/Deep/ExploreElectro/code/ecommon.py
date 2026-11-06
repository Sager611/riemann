# ecommon.py (explore22): completed functions with F(s) = F(1-s), real on the critical line (copied from explore18's topo_common).
#   xi : Riemann xi.  mv : xi with zeros 6, 7 replaced by 1/2 +- 0.15 + i*39.2524 (FE kept, no Euler product).
#   dh : Davenport-Heilbronn Lambda(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s), f = sum a_n n^-s, a = (1, k, -k, -1, 0) mod 5.
import mpmath as mp
mp.mp.dps = 25

def h(s):
    return s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)
def hlogd(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2
def xi(s):
    return h(s) * mp.zeta(s)
def xilogd(s):
    return hlogd(s) + mp.zeta(s, 1, 1) / mp.zeta(s)

GA = mp.mpf('37.586178158825671257217763481')   # 6th zero of zeta
GB = mp.mpf('40.918719012147495187398126914')   # 7th zero
GBAR = (GA + GB) / 2
DELTA = mp.mpf('0.15')
def mvnum():
    return [mp.mpc(0.5 + DELTA, GBAR), mp.mpc(0.5 + DELTA, -GBAR), mp.mpc(0.5 - DELTA, GBAR), mp.mpc(0.5 - DELTA, -GBAR)]
def mvden():
    return [mp.mpc(0.5, GA), mp.mpc(0.5, -GA), mp.mpc(0.5, GB), mp.mpc(0.5, -GB)]
def mvlogd(s):
    return xilogd(s) + sum(1 / (s - r) for r in mvnum()) - sum(1 / (s - r) for r in mvden())

def kappa():
    return (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
def acoef():
    k = kappa()
    return {1: mp.mpf(1), 2: k, 3: -k, 4: mp.mpf(-1)}
def dh_f(s):
    A = acoef()
    return mp.power(5, -s) * sum(A[r] * mp.zeta(s, mp.mpf(r) / 5) for r in A)
def dh_fp(s):
    A = acoef()
    return mp.power(5, -s) * sum(A[r] * (mp.zeta(s, mp.mpf(r) / 5, 1) - mp.log(5) * mp.zeta(s, mp.mpf(r) / 5)) for r in A)
def dh(s):
    return mp.power(5 / mp.pi, (s + 1) / 2) * mp.gamma((s + 1) / 2) * dh_f(s)
def dhlogd(s):
    return mp.log(5 / mp.pi) / 2 + mp.digamma((s + 1) / 2) / 2 + dh_fp(s) / dh_f(s)

LOGD = {'xi': xilogd, 'mv': mvlogd, 'dh': dhlogd}
# Known DH off-line zeros (Spira 1994; Balanzario-Sanchez-Ortiz 2007), re-checked by explore18/19.
DH_OFF = [(0.808517, 85.699348), (0.650830, 114.163343), (0.574356, 166.479306), (0.724258, 176.702461)]
