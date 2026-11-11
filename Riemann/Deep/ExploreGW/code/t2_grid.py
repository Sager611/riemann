
# T2 grid: Hardy Z of the even primitive characters mod 13 on [t0, T] (Arb via python-flint)
import sys, time, numpy as np
from flint import dirichlet_char, acb, ctx
ctx.prec = 64
T = float(sys.argv[1]); h = float(sys.argv[2])
labels = [12, 3, 9, 4, 10]           # quad; order-3 pair; order-6 pair
ch = {n: dirichlet_char(13, n) for n in labels}
# conjugate check: Z_9(t) vs Z_3(-t), Z_10(t) vs Z_4(-t)
for t in (17.3, 123.4):
    print("conj check t=%g: Z9=%.12f Z3(-t)=%.12f | Z10=%.12f Z4(-t)=%.12f" % (t,
      float(ch[9].hardy_z(acb(t)).real.mid()), float(ch[3].hardy_z(acb(-t)).real.mid()),
      float(ch[10].hardy_z(acb(t)).real.mid()), float(ch[4].hardy_z(acb(-t)).real.mid())), flush=True)
ts = np.arange(h, T + h/2, h)
Z = np.zeros((len(labels), len(ts)))
t0 = time.time()
for i, t in enumerate(ts):
    for j, n in enumerate(labels):
        Z[j, i] = float(ch[n].hardy_z(acb(float(t))).real.mid())
    if i % 4000 == 0:
        np.savez("t2_grid_partial.npz", ts=ts[:i+1], Z=Z[:, :i+1], labels=labels)
        print("i=%d t=%.2f elapsed %.1f s" % (i, t, time.time()-t0), flush=True)
np.savez("t2_grid.npz", ts=ts, Z=Z, labels=labels)
print("done", len(ts), "points, %.1f s" % (time.time()-t0))
