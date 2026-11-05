# t3c_smatrix.py — the T(2,5) vector (theta_{Re psi}, theta_{Im psi}) of weight 1/2 transforms under tau -> -1/tau by a real
# orthogonal reflection S (Lee–Yang modular data); its eigenvectors are the FE-eigencombinations G (+1) and H (-1) of t3_zeros.py.
import numpy as np
def psi20():
    a = np.zeros(20, complex)
    for n in range(20):
        if np.gcd(n, 20) != 1: continue
        p4 = 1 if n % 4 == 1 else -1
        k = [j for j in range(4) if pow(2, j, 5) == n % 5][0]
        a[n] = p4*(1j)**k
    return a
P = psi20(); n = np.arange(1, 400)
def th(y):  # theta_j(i y) = sum_{n>=1} chi_j(n) exp(-pi n^2 y / 20), j = Re, Im
    w = np.exp(-np.pi*n**2*y/20)
    return np.array([np.sum(P[n % 20].real*w), np.sum(P[n % 20].imag*w)])
ys = [0.3, 0.5, 0.7, 0.9, 1.3]
A = np.array([th(y) for y in ys]); Bm = np.array([th(1/y)/np.sqrt(y) for y in ys])
S, *_ = np.linalg.lstsq(A, Bm, rcond=None); S = S.T
print('S (basis Re psi, Im psi):\n', np.round(S, 12))
print('max residual:', np.max(np.abs(A @ S.T - Bm)))
print('Lee–Yang (2/sqrt5)[[sin72, sin36],[sin36, -sin72]]:\n', np.round(2/np.sqrt(5)*np.array([[np.sin(2*np.pi/5), np.sin(np.pi/5)], [np.sin(np.pi/5), -np.sin(2*np.pi/5)]]), 12))
w, V = np.linalg.eigh((S+S.T)/2)
th2 = np.angle(np.sum(P*np.exp(2j*np.pi*np.arange(20)/20))/np.sqrt(20))/2
print('eigenvalues:', np.round(w, 12), ' eigenvectors (columns):\n', np.round(V, 12))
print('(cos th/2, sin th/2) =', np.round([np.cos(th2), np.sin(th2)], 12), ' (-sin, cos) =', np.round([-np.sin(th2), np.cos(th2)], 12))
