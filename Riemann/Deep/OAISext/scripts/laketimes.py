import re,statistics
base='/Users/adriansager/playground/riemann/Riemann/Deep/OAISext/logs/'
def secs(s):
    m=re.match(r'(?:(\d+)m)?\s*([\d.]+)(ms|s)',s)
    if not m: return None
    v=float(m.group(2)); v=v/1000 if m.group(3)=='ms' else v
    return v+60*int(m.group(1) or 0)
for f in ['chunk_001.log','chunk_002.log']:
    ts=[]
    for l in open(base+f,encoding='utf-8'):
        m=re.search(r'Built (\S+) \(([^)]*)\)',l)
        if m:
            v=secs(m.group(2))
            if v is not None: ts.append(v)
    print(f,'n',len(ts),'mean %.1f s'%statistics.mean(ts),'median %.1f s'%statistics.median(ts),'min %.1f'%min(ts),'max %.1f'%max(ts),'sum %.0f s'%sum(ts))
