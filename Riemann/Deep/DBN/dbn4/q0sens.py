
import sys, math, json
sys.path.insert(0,'.')
exec(open('sens.py').read().split('U0=replay')[0])
for q0 in (1.000002, 0.99, 0.9752, 0.95):
    U=replay(S,K57,rate_fn=make_rate(ELL0),q0=q0)
    print('q0=%.4f  U=%.6f'%(q0,U))
