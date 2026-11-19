"""explore51: exact heat-flow zero dynamics for periodic configurations (beta = infinity Dyson model).

F(z) = prod_j 2 sin(pi (z - z_j)/P), N zeros per period P (N even), is a trigonometric polynomial
F = sum_m c_m e^{2 pi i m z/P}, |m| <= N/2.  The de Bruijn-Newman flow d_t F = -F'' acts by
c_m(t) = c_m e^{(2 pi m/P)^2 t}; its zeros obey z_k' = 2 sum_j (pi/P) cot(pi (z_k - z_j)/P), i.e. the
zero ODE of H_t for the P-periodic continuation of the configuration.  So the flow is exact (no ODE step):
hitting time = first t with N sign changes of F_t on one period (bisection).
Density: rho = log(X/4 pi)/(4 pi) at Raymond's wall X = 6000000185827 (z = 2 gamma units), spacing s = 1/rho.
"""
import numpy as np, math, json, sys, time

X = 6000000185827.0
LOGT = math.log(X / (4 * math.pi))          # = log(T0/2pi), T0 = X/2
RHO = LOGT / (4 * math.pi)
S = 1 / RHO
TH141 = 0.8733458561067185768
Y0 = 2 * TH141 - 1


def t_row(rho, y0):
    c = math.pi * rho
    return math.log(math.cosh(c * y0)) / c ** 2


def coeffs(zs, P):
    """Fourier coefficients c_m, m = -N/2..N/2 (array index m + N/2), of prod 2 sin(pi(x - z_j)/P)."""
    N = len(zs)
    M = 4 * N
    x = np.arange(M) * P / M
    lf = np.zeros(M, dtype=complex)
    for zj in zs:
        lf += np.log(2 * np.sin(np.pi * (x - zj) / P))
    lf -= lf.real.max()
    f = np.exp(lf)
    f = f * np.exp(-1j * np.angle(f[np.argmax(abs(f))]))   # make real up to rounding
    c = np.fft.fft(f) / M
    m = np.arange(-N // 2, N // 2 + 1)
    return m, c[m % M]


def evolve(m, c, P, t):
    return c * np.exp((2 * np.pi * m / P) ** 2 * t)


def n_real(m, ct, P, M=2 ** 17):
    """number of sign changes of F_t on one period (= number of simple real zeros, if resolved)."""
    full = np.zeros(M, dtype=complex)
    full[m % M] = ct
    f = np.fft.ifft(full).real * M
    sg = np.sign(f)
    return int(np.sum(sg != np.roll(sg, 1)))


def roots(m, ct, P):
    """zeros z (mod P) of the trig polynomial: polynomial in w = e^{2 pi i z/P}, degree N."""
    w = np.roots(ct[::-1])
    return P / (2j * np.pi) * np.log(w)


def hit_time(zs, P, thi=None, tol=1e-7, M=2 ** 17):
    m, c = coeffs(zs, P)
    N = len(zs)
    ymax = max(abs(np.imag(zs)))
    lo, hi = 0.0, (thi if thi else ymax ** 2 / 2 + 1e-6)
    if n_real(m, evolve(m, c, P, hi), P, M) < N:
        return float('nan')
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if n_real(m, evolve(m, c, P, mid), P, M) >= N:
            hi = mid
        else:
            lo = mid
    return hi


def ode_rk4(zs, P, tend, dt):
    """direct N-body check: z_k' = 2 sum_j (pi/P) cot(pi (z_k - z_j)/P)."""
    z = np.array(zs, dtype=complex)
    def f(z):
        d = z[:, None] - z[None, :]
        np.fill_diagonal(d, 1.0)
        k = 2 * (np.pi / P) / np.tan(np.pi * d / P)
        np.fill_diagonal(k, 0.0)
        return k.sum(axis=1)
    t = 0.0
    out = [(t, max(abs(z.imag)))]
    while t < tend - 1e-15:
        h = min(dt, tend - t)
        k1 = f(z); k2 = f(z + h / 2 * k1); k3 = f(z + h / 2 * k2); k4 = f(z + h * k3)
        z = z + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
        out.append((t, max(abs(z.imag))))
    return z, out


def env_lattice(N):
    return (np.arange(N) + 0.5) * S


def env_cue(N, seed):
    rng = np.random.default_rng(seed)
    A = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / math.sqrt(2)
    Q, R = np.linalg.qr(A)
    Q = Q * (np.diag(R) / abs(np.diag(R)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)) % (2 * np.pi))
    return th / (2 * np.pi) * N * S


def insert_row(env, k, y):
    """replace 2k consecutive real zeros (centred) by k conjugate pairs at the midpoints of consecutive zeros."""
    N = len(env)
    i0 = N // 2 - k
    rem = env[i0:i0 + 2 * k]
    mids = (rem[0::2] + rem[1::2]) / 2
    keep = np.concatenate([env[:i0], env[i0 + 2 * k:]])
    return np.concatenate([keep.astype(complex), mids + 1j * y, mids - 1j * y]), mids


if __name__ == '__main__':
    t0 = time.time()
    res = {'rho': RHO, 's': S, 'logT0over2pi': LOGT, 'y0': Y0, 't_row_formula': t_row(RHO, Y0),
           'deBruijn': Y0 ** 2 / 2}
    print(f"rho={RHO:.6f} s={S:.6f} L={LOGT:.6f} y0={Y0:.12f}")
    print(f"formula t_row={t_row(RHO, Y0):.9f}  de Bruijn y0^2/2={Y0**2/2:.9f}")
    # (A) smallest case: one pair per period a = 2s (N = 2): exact formula
    a = 2 * S
    tA = hit_time(np.array([1j * Y0, -1j * Y0]), a, tol=1e-9, M=2 ** 14)
    print(f"(A) N=2 row, trig flow: {tA:.9f}")
    res['A_N2'] = tA
    # (B) N = 4: one pair + 2 real zeros per period 4s; trig flow vs RK4 ODE
    P4 = 4 * S
    z4 = np.array([0.5 * S, 2.5 * S, 1.5 * S + 1j * Y0, 1.5 * S - 1j * Y0])
    tB4 = hit_time(z4, P4, tol=1e-9, M=2 ** 14)
    m4, c4 = coeffs(z4, P4)
    zr, hist = ode_rk4(z4, P4, 0.02, 1e-6)
    zt = roots(m4, evolve(m4, c4, P4, 0.02), P4)
    print(f"(B) N=4 hit (trig) {tB4:.9f}; max Im at t=0.02: ODE {hist[-1][1]:.10f} trig {max(abs(zt.imag)):.10f}")
    res['B_N4'] = {'hit': tB4, 'ymax_t002_ode': hist[-1][1], 'ymax_t002_trig': float(max(abs(zt.imag)))}
    # (C) N = 400, environments x k pairs at height y0
    N = 400
    P = N * S
    envs = {'lattice': env_lattice(N), 'cue1': env_cue(N, 1), 'cue2': env_cue(N, 2)}
    try:
        zz = json.load(open('Riemann/Deep/ExploreRelax/logs/zeta_unfolded.json'))
        envs['zeta'] = np.array(zz['x'])
    except Exception as e:
        print('no zeta env:', e)
    ks = [1, 2, 4, 8, 16, 32, 64, 100, 150, 200]
    res['C'] = {}
    for name, env in envs.items():
        row = {}
        for k in ks:
            zs, mids = insert_row(env, k, Y0)
            th = hit_time(zs, P)
            row[k] = th
            print(f"(C) {name:8s} k={k:3d}: hit {th:.6f}", flush=True)
        res['C'][name] = row
    # (D) decay curves y_max(t) for k = 1 and the full row, CUE env
    ts = [0, 0.0125, 0.025, 0.05, 0.075, 0.1]
    res['D'] = {}
    for k in [1, 200]:
        zs, mids = insert_row(envs['cue1'], k, Y0)
        m, c = coeffs(zs, P)
        cur = []
        for t in ts:
            zt = roots(m, evolve(m, c, P, t), P)
            cur.append(float(max(abs(zt.imag))))
        res['D'][k] = cur
        print(f"(D) cue1 k={k}: y_max(t) at {ts}: " + ' '.join(f'{v:.5f}' for v in cur))
    json.dump(res, open('Riemann/Deep/ExploreRelax/logs/row_exact.json', 'w'), indent=1)
    print(f"done {time.time()-t0:.1f}s")
