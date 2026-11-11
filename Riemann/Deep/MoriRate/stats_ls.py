
import json, math, sys, os
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs')
def g2(y):
    y=abs(y); return math.log((2+math.sqrt(4-y*y))/y) if y<2 else 0.0
for fn in ['ls_c6_n15.log','ls_c10.log','ls_c14.log']:
    rows=[json.loads(l) for l in open(os.path.join(HERE, fn)) if l.startswith('{') and '"k"' in l]
    for sec in ('even','odd'):
        R=[r for r in rows if r['sector']==sec]
        mid=[r['logratio']-4*g2(r['y']) for r in R if 0.3<=r['y']<=1.8]
        lo=[r for r in R if r['y']<0.3]
        hi=[r['logratio'] for r in R if r['y']>=2.2]
        print(fn, sec, 'n_mid=%d mean_off=%.3f min=%.3f max=%.3f'%(len(mid), sum(mid)/len(mid), min(mid), max(mid)),
              'first(y=%.3f): meas %.3f pred %.3f'%(R[0]['y'], R[0]['logratio'], 4*g2(R[0]['y'])),
              'y>=2.2: n=%d mean=%.3f max=%.3f'%(len(hi), sum(hi)/len(hi), max(hi)))
