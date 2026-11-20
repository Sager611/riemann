
import os,re,sys
base = sys.argv[1]
out = []
for f in sorted(os.listdir(base)):
    if not f.endswith('.lean'): continue
    L = open(os.path.join(base,f)).read().split('\n')
    imps = [l.split()[1].split('.')[-1] for l in L if l.startswith('import ')]
    out.append('### %s (%d lines) imports: %s' % (f, len(L), ', '.join(imps)))
    i = 0
    while i < len(L):
        l = L[i]
        m = re.match(r'^(private |protected )?(theorem|lemma|def|structure|abbrev|noncomputable def|instance|class)\s', l)
        if m:
            blk = [l]
            j = i
            while not re.search(r':=\s*(by)?\s*$|:= by|where\s*$|:=\s*$', L[j]) and j < i+25 and j+1 < len(L):
                j += 1
                blk.append(L[j])
            out.append('%d: ' % (i+1) + '\n    '.join(x.strip() for x in blk)[:1500])
            i = j+1
            continue
        i += 1
open(sys.argv[2],'w').write('\n'.join(out))
print(len(out))
