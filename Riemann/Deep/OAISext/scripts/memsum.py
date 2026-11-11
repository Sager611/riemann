import re
mx={}
for l in open('/Users/adriansager/playground/riemann/Riemann/Deep/OAISext/logs/samples.log'):
    m=re.match(r'(\S+) (\d+) lean (\d+) mem (\d+)M(?: faults (\d*))?',l.strip())
    if not m: continue
    c=m.group(2); mem=int(m.group(4)); f=m.group(5)
    d=mx.setdefault(c,{'n':0,'max':0,'over450':0,'faults':[]})
    if mem>0: d['n']+=1; d['max']=max(d['max'],mem); d['over450']+= mem>450
    if f: d['faults'].append(int(f))
for c,d in mx.items():
    fs=sorted(d['faults'])
    print(c,'samples',d['n'],'max mem',d['max'],'M; >450M',d['over450'],'; faults at sample: n',len(fs),'median',fs[len(fs)//2] if fs else None,'max',fs[-1] if fs else None)
