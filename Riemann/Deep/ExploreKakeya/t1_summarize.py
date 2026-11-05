
import re
rows={}
for l in open('t1_out_raw.txt'):
    m=re.match(r'(\S+) s=(\S+) A=(\S+)',l)
    if m: rows.setdefault(m.group(1),[]).append((float(m.group(2)),float(m.group(3))))
sig=[s for s,_ in rows['H0']]
print('sigma  '+' '.join('%6.2f'%s for s in sig[::2]))
for H,v in rows.items():
    print('%-6s '%H+' '.join('%6.3f'%a for _,a in v[::2]))
for H,v in rows.items():
    s,a=max(v,key=lambda x:x[1]); d=dict(v); print(H,'sup A on [0.60,0.84] = %.4f at s=%.2f'%(a,s), ' A(0.70)=%.4f A(0.75)=%.4f'%(d[0.70],d[0.75]))
print('GM 15/(3+5s) at .7,.75,.8:', [round(15/(3+5*s),4) for s in (0.7,0.75,0.8)], ' 30/13=%.4f'%(30/13))
