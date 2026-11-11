# PYAudit: independent mpmath re-implementation of Yang's A_k, B_k (k = 3, 4, 5) from the PROOF of Yang 2301.03165 Lem 2.4-2.5
# (v1 l. 330-401 and 403-536), compared with SubWeyl/certify.py's AkBk; and the effect of two misreadings of the recursion.
# uv run --no-project --with mpmath --with python-flint --with numpy python3 akbk_check.py
import os, sys
import mpmath as mp
mp.mp.dps = 50
pi = mp.pi
def A3B3(eta, h):
    lam0 = (1/eta + 32*mp.sqrt(eta)*h/(15*mp.sqrt(pi)))**(-3)
    kappa = mp.sqrt(1 + mp.mpf(3)/8*mp.sqrt(pi)*eta**1.5)/2 - mp.mpf(1)/2
    d3 = mp.sqrt(1 + kappa)                        # delta_3 of the third-derivative test (= sqrt(1/2 + 1/2 sqrt(...)))
    l03 = mp.cbrt(lam0)
    A3 = mp.sqrt(1/(eta*h) + 32/(15*mp.sqrt(pi))*mp.sqrt(eta + l03) + (eta + l03)*l03/3)*d3
    B3 = mp.sqrt(32)/(mp.sqrt(3)*pi**0.25*eta**0.25)*d3
    return A3, B3, d3
def step(A, B, j, eta, h, mode):
    J = 2**(j - 1)
    if mode == 'Kmisread': J = 2**j                 # misreading: K of the conclusion (A_{k+1}) used in the step
    d = mp.sqrt(1 + 2/mp.mpf(2337)**(1 - mp.mpf(2)/J)*(9*pi/1024*eta)**(mp.mpf(1)/J))   # bound for D_j = sqrt(1+q/N), eq. (k_deriv_test_D_j_bound)
    if mode == 'delta3alt' and j == 3: d = A3B3(eta, h)[2]   # misreading: delta_3 of Lemma 2.4 reused in the recursion
    alpha = mp.mpf(4)*(J - 1)**2/((2*J - 1)*(4*J - 3)); beta = mp.mpf(4)*(J - 1)**2/((2*J - 3)*(4*J - 5))
    A1 = d*(h**(-mp.mpf(1)/J) + 2**(mp.mpf(7)/12)*mp.sqrt(alpha*A))    # = delta_j (h^{-1/J} + 2^{19/12}(J-1)/sqrt((2J-1)(4J-3)) A^{1/2})
    B1 = d*mp.sqrt(2*beta*B)
    return A1, B1
def AkBk(k, eta, h, mode='proof'):
    A, B, _ = A3B3(eta, h)
    for j in range(3, k): A, B = step(A, B, j, eta, h, mode)
    return A, B
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SubWeyl'))
sys.dont_write_bytecode = True   # keep SubWeyl/ free of __pycache__
import certify as C
from flint import arb
H = mp.mpf('1.012'); amin = mp.mpf(2336)/(2*H*H*(H - 1)); kap = (amin/(amin - 1))**4.5
HC5 = 76545*mp.sqrt(2)/107264
cases = [(3, '1.12054', H**4), (4, '0.87152', H**4), (4, '0.5', mp.mpf('1.0513')**4), (5, '1.12054', HC5*H**9*kap), (5, '0.75', HC5*mp.mpf('1.0513')**9), (5, '2.0', HC5*H**9)]
worst = 0
for k, et, h in cases:
    eta = mp.mpf(et)
    A, B = AkBk(k, eta, h)
    Ac, Bc = C.AkBk(k, arb(et), arb(mp.nstr(h, 40)))
    ra = abs(A/mp.mpf(Ac.mid().str(40, radius=False)) - 1); rb = abs(B/mp.mpf(Bc.mid().str(40, radius=False)) - 1)
    worst = max(worst, ra, rb)
    alt = AkBk(k, eta, h, 'delta3alt'); km = AkBk(k, eta, h, 'Kmisread')
    print('k=%d eta=%s h=%s: A=%s B=%s | rel.diff vs certify %.1e %.1e | delta3-alt reading A x%s B x%s | K-misread A x%s B x%s' % (
        k, et, mp.nstr(h, 8), mp.nstr(A, 12), mp.nstr(B, 12), ra, rb, mp.nstr(alt[0]/A, 6), mp.nstr(alt[1]/B, 6), mp.nstr(km[0]/A, 6), mp.nstr(km[1]/B, 6)))
print('max relative difference mpmath vs certify.py: %.2e' % worst)
