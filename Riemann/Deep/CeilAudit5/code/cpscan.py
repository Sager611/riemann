import re
L = open('.scratch/ceilaudit5/dl/CP.md').read().split(chr(10))
starts = [(15023,'H'),(15924,'I'),(16596,'J'),(16860,'K'),(17093,'L'),(17411,'M'),(17767,'N'),(17985,'O'),(18310,'P'),(18492,'Q'),(18575,'END')]
pat = re.compile(r'(gain over|The gain|pays? the return|Two changes|The improvement|new (?:source|input|estimate)|only change|changes? (?:pay|from)|improves|The new |is paid by|surplus)')
for (a,name),(b,_) in zip(starts, starts[1:]):
    seg = L[a-1:b-1]
    hits = [(a+i, l.strip()[:400]) for i,l in enumerate(seg) if pat.search(l)]
    print('==', name, len(seg), 'lines')
    for h in hits[:3]: print('  %d: %s' % h)
