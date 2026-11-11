# fourier.py -- lamz2audit: numerical check of the Fourier identities in the proof of Lemma 3.2 of 2609.02882 (TeX l.735-790)
# at complex arguments: FT(Q'') = (2 pi i z)^2 FT(Q), and K(z)^2 = FT(r_{T})(z) w(rho-rho') at z = i(rho-rho') log T/(2 pi).
# Test f = c (1-4x^2)^6 on [-1/2,1/2] (f = eta^2, eta = sqrt(c)(1-4x^2)^3, C^5 with 6 vanishing derivatives at the ends).
# Run: uv run --no-project --with mpmath python3 fourier.py > logs/fourier.log
import mpmath as mp, random
mp.mp.dps = 30
random.seed(3)
c = 1/mp.quad(lambda x: (1-4*x**2)**6, [-0.5, 0.5])
f   = lambda x: c*(1-4*x**2)**6
f2  = lambda x: mp.diff(lambda t: c*(1-4*t**2)**6, x, 2)
FT  = lambda g, z: mp.quad(lambda x: g(x)*mp.exp(-2j*mp.pi*z*x), [-0.5, 0, 0.5])
print("int f =", mp.nstr(mp.quad(f, [-0.5, 0.5]), 20), " (K(0)=1)")
worst1 = 0; worst2 = 0
for trial in range(8):
    z = mp.mpc(random.uniform(-6, 6), random.uniform(-2, 2))
    K = FT(f, z); F2 = FT(f2, z)
    e1 = abs(F2 - (2j*mp.pi*z)**2*K)/max(1, abs(F2)); worst1 = max(worst1, e1)
    # Q = f*f, Q'' = f''*f  =>  FT(Q) = K^2, FT(Q'') = F2*K;  r_T = Q - Q''/(4 log^2 T)
    T = mp.mpf(10)**random.uniform(2, 30); L = mp.log(T)
    d = mp.mpc(random.uniform(-0.9, 0.9), random.uniform(-40, 40))          # rho - rho' : |Re| < 1
    zz = 1j*d*L/(2*mp.pi)
    Kz = FT(f, zz); F2z = FT(f2, zz)
    rhat = Kz**2 - F2z*Kz/(4*L**2)
    w = 4/(4 - d**2)
    e2 = abs(rhat*w - Kz**2)/max(1, abs(Kz**2)); worst2 = max(worst2, e2)
    print(" z=%s  |FT(f'')-(2 pi i z)^2 FT(f)|/max(1,|.|)=%s ;  rho-rho'=%s T=1e%.1f  |rhat*w - K^2|/max(1,|K^2|)=%s" % (
        mp.nstr(z, 6), mp.nstr(e1, 3), mp.nstr(d, 6), float(mp.log10(T)), mp.nstr(e2, 3)))
print("worst relative errors:", mp.nstr(worst1, 3), mp.nstr(worst2, 3))
