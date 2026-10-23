#!/usr/bin/env python3
"""Rebuild instance.json from the arXiv 2610.10108 e-print.
Usage: curl -s -L -o e.tgz https://arxiv.org/e-print/2610.10108 && mkdir src && tar xzf e.tgz -C src
       python3 parse_eprint.py src
"""
import json, os, re, sys

d = sys.argv[1]
src = open(os.path.join(d, "a1-missing-proof.tex")).read()
causal = open(os.path.join(d, "4-causal.tex")).read()
strip = lambda b: re.sub(r"\\begin\{aligned\}|\\end\{aligned\}", "", b)
nums = lambda b: [int(x) for x in re.findall(r"\d+", b)]
R = nums(strip(re.search(r"R=\\left\\{(.*?)\\right\\}", src, re.S).group(1)))
P = nums(strip(re.search(r"\\pi=\\left\((.*?)\\right\)", src, re.S).group(1)))
C = nums(re.search(r"C&=\\\{([^}]*)\\\}", causal).group(1))
J = nums(re.search(r"J&=\\\{([^}]*)\\\}", causal).group(1))
assert len(set(R)) == 175 and sorted(P) == list(range(182)) and len(C) == len(J) == 9
inst = {"source": "arXiv 2610.10108 e-print (a1-missing-proof.tex, 4-causal.tex)",
        "field": "F9 = F3[theta]/(theta^2+1); element index k = a + 3b for a + b*theta; eta_r = (r mod 3) + floor(r/3)*theta, so eta_r has index r",
        "T": "v_0=(0,0,1), v_{1+r}=(0,1,eta_r), v_{10+9r+s}=(1,eta_r,eta_s)",
        "labels": "P_i -> i, L_j -> 91+j; e_{10i+r} = {i, 91+j_{i,r}} with j_{i,0}<...<j_{i,9}",
        "R": R, "pi": P, "C": C, "J": J}
json.dump(inst, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "instance.json"), "w"), indent=0)
print("wrote instance.json")

