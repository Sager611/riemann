# nb2fals f4: export per-N values (N in [2,640] csv; [2,5500] csv.gz) for C2, kappa, J - Sigma_N, R, d_N^2.
import numpy as np, gzip
O = '/Users/adriansager/playground/riemann/.scratch/nb2fals/'; D = '/Users/adriansager/playground/riemann/Riemann/Deep/NB2Fals/data/'
import os; os.makedirs(D, exist_ok=True)
P = np.load(O + 'perN5500.npz')['perN']
hdr = 'N,C2=L^2d^2-CL,kappa=S_N-Sigma_N,J-Sigma_N,R=L^2I_N-CL-J,d_N^2'
rows = [f'{int(r[0])},{r[2]:.10f},{r[3]:.10f},{r[4]:.10f},{r[5]:.10f},{r[6]:.16e}' for r in P]
open(D + 'perN_640.csv', 'w').write(hdr + '\n' + '\n'.join(rows[:639]) + '\n')
with gzip.open(D + 'perN_5500.csv.gz', 'wt') as f: f.write(hdr + '\n' + '\n'.join(rows) + '\n')
print(len(rows), rows[638][:60])
