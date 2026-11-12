
import re
out=[]
for f in ["LemmaZ2","LemmaZW","WiredB"]:
    s=open(f"Riemann/Deep/DPCBeyond/Lean/{f}.lean").read()
    lines=s.split("\n"); i=0
    out.append(f"-- ===== {f}.lean")
    while i<len(lines):
        L=lines[i]
        if re.match(r'^(theorem|lemma|def|noncomputable def|structure) ',L):
            blk=[L]
            if L.startswith("structure"):
                j=i+1
                while j<len(lines) and lines[j].strip()!="":
                    blk.append(lines[j]); j+=1
                out.append("\n".join(blk)); i=j; continue
            j=i
            while ":= by" not in lines[j] and not re.search(r":=\s*$",lines[j]) and j<len(lines)-1:
                j+=1; blk.append(lines[j])
            if L.startswith("def") or L.startswith("noncomputable def"):
                k=j+1
                while k<len(lines) and lines[k].strip()!="":
                    blk.append(lines[k]); k+=1
            txt="\n".join(blk); txt=re.sub(r"\s*:= by\s*$","",txt)
            out.append(txt); i=j+1; continue
        i+=1
open(".scratch/dpcblean/statements_new.txt","w").write("\n".join(out)+"\n")
p=open("Riemann/Deep/DPCBeyond/PROOF.md").read()
a=p.index("## Lemma Z"); b=p.index("## Step 6e′")
c=p.index("(a) Lemma Z with"); e=p.index("(b) y-integral")
open(".scratch/dpcblean/proof_parts.txt","w").write(p[a:b]+"\n...\n"+p[c:e])
print(len("\n".join(out)), len(p[a:b])+len(p[c:e]))

