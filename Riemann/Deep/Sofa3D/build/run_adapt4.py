import sys, json, time, math, argparse, numpy as np
sys.path.insert(0, '.')
import sofa2 as S, adapt4
ap = argparse.ArgumentParser()
ap.add_argument('--path', default='w5_night2'); ap.add_argument('--nc0', type=int, default=50); ap.add_argument('--L', type=int, default=3)
ap.add_argument('--tau_rel', type=float, default=0.03); ap.add_argument('--tau_abs', type=float, default=1e-9)
ap.add_argument('--tt', type=float, default=0.005); ap.add_argument('--tr', type=float, default=0.05); ap.add_argument('--off', type=float, default=0.38196601125)
ap.add_argument('--gmerge', type=float, default=0.02); ap.add_argument('--save', default='')
ap.add_argument('--chunks', type=int, default=1); ap.add_argument('--minvol', type=float, default=0.0)
a = ap.parse_args()
P = S.load_path(a.path if a.path.endswith('.json') else '../gist/%s_path.json' % a.path)
sm = S.smooth(P, a.tt, a.tr)
C, th = S.checks_linear(sm)
t0 = time.time()
ap2 = int(a.chunks); edges = np.linspace(0, a.nc0, ap2 + 1).astype(int)
tot = []; subs = []; cls = []
for ci in range(ap2):
    for ck in range(ap2):
        t_, s_, c_, st = adapt4.certify(C, a.nc0, a.L, a.tau_rel, a.tau_abs, off=a.off, gmerge=a.gmerge, log=lambda m: None,
                                        jr=(int(edges[ci]), int(edges[ci + 1])), kr=(int(edges[ck]), int(edges[ck + 1])))
        vv = 4 * s_[:, 2] * s_[:, 3] * (s_[:, 7] - s_[:, 4]); keep = vv >= a.minvol
        tot.append(t_); subs.append(s_[keep]); cls.append(c_[keep])
        print('chunk', ci, ck, 'cert %.6f kept %d of %d' % (t_, int(keep.sum()), len(s_)), flush=True)
total = math.fsum(tot); sub = np.vstack(subs); cells = np.vstack(cls)
v = math.fsum((4 * sub[:, 2] * sub[:, 3] * (sub[:, 7] - sub[:, 4])).tolist())
print(json.dumps(dict(total=total, subvol=v, nsub=len(sub), wall=round(time.time() - t0, 1))))
if a.save: np.savez(a.save, sub=sub, cells=cells, meta=np.array([a.nc0, a.L, a.off, a.tt, a.tr]))
