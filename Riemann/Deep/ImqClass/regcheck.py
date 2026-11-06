import sys
t4={}
for l in open('runs/reg2000.t4'):
    d,M=l.split(); t4[int(d)]=int(M)
lst={}
for l in open('runs/reg2000.list'):
    d,h=l.split(); lst[int(d)]=int(h)
missing=[d for d in lst if d not in t4]
print('list fields', len(lst), 'in t4', len(lst)-len(missing), 'MISSING', len(missing), missing[:5])
bad=[d for d in lst if d in t4 and t4[d]>lst[d]]
print('M>h violations', len(bad), bad[:5])
extra=[d for d in t4 if d not in lst]
print('t4 not in list (need h>2000):', len(extra))
open('runs/reg2000.extra','w').write(''.join(str(d)+'\n' for d in extra))
