"""explore51: 400 consecutive zeta zeros (n = 1001..1400, mpmath.zetazero), unfolded by theta(t)/pi + 1 and
rescaled to the spacing s = 1/rho(X) of row_exact.py: a real environment with zeta's own local statistics."""
import mpmath as mp, json, math
X = 6000000185827.0
RHO = math.log(X / (4 * math.pi)) / (4 * math.pi); S = 1 / RHO
mp.mp.dps = 20
n0, N = 1001, 400
g = [mp.zetazero(n).imag for n in range(n0, n0 + N)]
u = [float(mp.siegeltheta(t) / mp.pi + 1) for t in g]
x = [(ui - u[0] + 0.5) * S for ui in u]
gaps = [u[i + 1] - u[i] for i in range(N - 1)]
print('gamma range', float(g[0]), float(g[-1]), 'unfolded span', u[-1] - u[0], 'min gap', min(gaps), 'max gap', max(gaps))
json.dump({'n0': n0, 'N': N, 'gamma': [float(t) for t in g], 'x': x, 'S': S}, open('Riemann/Deep/ExploreRelax/logs/zeta_unfolded.json', 'w'))
