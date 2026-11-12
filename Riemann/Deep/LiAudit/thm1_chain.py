# Theorem 1 at finite N with real zeros and real lambda_n: every inequality of P's steps 2-7, and the trend lemma (step 4).
# usage: python3 thm1_chain.py zeros1 R_fft.npy
import sys, numpy as np
from math import log, pi, sqrt
EG = 0.57721566490153286061; kap = 1 + log(2*pi) - EG
g = np.loadtxt(sys.argv[1])[:2000]; R = np.load(sys.argv[2])
th_all = 2*np.arctan(1/(2*g))
def smoothstep(y):
    y = np.clip(y, 0, 1); a = np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0); b = np.where(y < 1, np.exp(-1/np.maximum(1-y, 1e-300)), 0.0)
    return a/(a+b)
def make_phi(delta, M=2**15):
    x = np.linspace(-0.5, 0.5, M+1); h = x[1]-x[0]
    psi = smoothstep((0.5-np.abs(x))/delta); psi /= sqrt((psi**2).sum()*h)
    L = 2*M+1; F = np.fft.rfft(psi, 2*L); phi = np.fft.irfft(F*F, 2*L)[:2*M+1]*h   # phi on [-1,1], step h
    xs = np.linspace(-1, 1, 2*M+1); return xs, phi, (psi.sum()*h)**2
def chain(N, delta, u=0.6, Tover=None):
    xs, phi, intphi = make_phi(delta)
    n = np.arange(0, N+1); w = np.interp(n/N, xs, phi)          # phi(n/N), n >= 0 (phi even)
    Phi0 = w[0] + 2*w[1:].sum()
    if Tover is None:                                              # T^2 log(T/2pi) = u N int phi
        T = 30.0
        for _ in range(200): T = sqrt(u*N*intphi/log(T/(2*pi)))
    else: T = Tover
    th = th_all[g <= T]; NT = th.size
    nn = n[1:].astype(float); tn = 0.5*nn*(np.log(nn) - kap)
    E = np.exp(1j*np.outer(nn, th))                                # e^{i n theta_gamma}, n >= 1
    W = E.real.sum(axis=1); Z = E.sum(axis=1)
    RR = R[1:N+1]
    Rn2 = 2*(w[1:]*RR**2).sum()                                    # ||R~||^2
    RW = 2*(w[1:]*RR*W).sum()                                      # <R~, W>
    W2 = w[0]*NT**2 + 2*(w[1:]*W**2).sum()                         # ||W||^2
    Pm = w[0]*NT**2 + 2*(w[1:]*np.abs(Z)**2).sum()                 # sum_{g,g'<=T} Phi(th-th')
    Pp = w[0]*NT**2 + 2*(w[1:]*(Z**2).real).sum()                  # sum_{g,g'<=T} Phi(th+th')
    tau = 2*(w[1:,None]*tn[:,None]*E.real).sum(axis=0)             # tau(theta_gamma), exact
    tauT = -(np.log(2*pi*th)/th**2)                                 # P's main term per zero
    X = RW + tau.sum()                                              # = sum_g sum_g' B(th_g, th_g')
    taup = -tau.sum()                                               # tau' (exact)
    bound = Phi0*NT*(1 - taup/(Phi0*NT))**2 if 2*0.5*Phi0*NT > taup else float('nan')
    S = (RR**2).sum()
    print('N=%d delta=%.2f T=%.2f N(T)=%d int(phi)=%.4f Phi(0)/N=%.4f u_eff=%.3f' % (N, delta, T, NT, intphi, Phi0/N, T**2*log(T/(2*pi))/(N*intphi)))
    print('  ||W||^2 - Pm/2 = %.3e (Pp/2 = %.3e);  X + Pm = %.4e (claim <= tiny; Pm = %.4e)' % (W2 - Pm/2, Pp/2, X + Pm, Pm))
    print('  CS: <R,W>^2/||W||^2 = %.4e <= ||R||^2 = %.4e;  <R,W> = %.4e, -(2||W||^2 - tau\') = %.4e' % (RW**2/W2, Rn2, RW, -(2*W2 - taup)))
    print('  trend: max_g |tau - log(2pi th)/th^2| = %.3f (log N = %.2f) at th=%.4f (N^-1/2=%.4f); tau\'/tau_T = %.5f; tau\'/(Phi0 N(T)) = %.4f' % (np.abs(-tau - tauT).max(), log(N), th[np.abs(-tau-tauT).argmax()], N**-0.5, taup/tauT.sum(), taup/(Phi0*NT)))
    print('  bound Phi0 N(T)(1-tau\'/(Phi0 N(T)))^2 = %.4e;  sum_{n<=N} R^2 = %.4e (ratio %.3f);  c1-asymptotic %.4e' % (bound, S, S/bound, 0.05579*N**1.5*sqrt(log(N))))
for (N, delta) in ((1000, 0.1), (10**4, 0.1), (10**5, 0.1), (10**5, 0.02)):
    chain(N, delta)
chain(10**5, 0.1, Tover=0.5*sqrt(10**5))                            # corner: T = N^{1/2}/2, the largest P allows
chain(10**4, 0.1, Tover=0.5*sqrt(10**4))
