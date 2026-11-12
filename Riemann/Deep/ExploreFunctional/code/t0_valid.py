import sys, numpy as np, math, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form
t0 = time.time()
for c in [1.2, 2.0]:
    F = Form(c)
    for sec in ['even', 'odd']:
        for n in [10, 16, 24]:
            Qm, p = F.Q(n, sec)
            e = np.linalg.eigvalsh(Qm)
            ea = np.linalg.eigvalsh(p['arch'])
            print(c, sec, n, 'Q:', ['%.6e' % v for v in e[:3]], 'arch:', '%.6e' % ea[0], flush=True)
print('sec', round(time.time()-t0, 1))
