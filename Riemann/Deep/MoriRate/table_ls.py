
import json, math, sys
def g2(y):
    y=abs(y); return math.log((2+math.sqrt(4-y*y))/y) if y<2 else 0.0
for fn in sys.argv[1:]:
    rows={}
    head=None
    for l in open(fn):
        l=l.strip()
        if not l.startswith('{'): continue
        d=json.loads(l)
        if 'N' in d: head=d
        if 'E1' in d: print(fn.split('/')[-1], d)
        if 'k' in d: rows.setdefault(d['k'],{})[d['sector']]=d
    print(head)
    print(' k  gamma    y     meas_e  meas_o  4g2    diff_e diff_o')
    for k in sorted(rows):
        r=rows[k]; e=r.get('even',{}); o=r.get('odd',{})
        y=e.get('y',o.get('y'))
        p=4*g2(y)
        me=e.get('logratio',float('nan')); mo=o.get('logratio',float('nan'))
        if k<=12 or k%3==0:
            print('%3d %7.3f %.3f %7.3f %7.3f %7.3f %6.3f %6.3f'%(k, e.get('gamma',0), y, me, mo, p, me-p, mo-p))
