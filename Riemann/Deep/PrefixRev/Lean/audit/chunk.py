import sys
for name in ["23_18_4", "23_16_6", "23_12_10"]:
    L, code = open(".scratch/prefixlean/word_%s.txt" % name).read().split()
    L = int(L); code = int(code)
    digits = []
    for _ in range(L):
        digits.append(code % 3); code //= 3
    assert code == 0
    chunks = []
    for i in range(0, L, 20):
        ch = digits[i:i+20]
        c = sum(d * 3**j for j, d in enumerate(ch))
        chunks.append((len(ch), c))
    s = ", ".join("(%d, %d)" % t for t in chunks)
    open(".scratch/prefixlean/chunks_%s.txt" % name, "w").write(s)
    print(name, len(chunks))
