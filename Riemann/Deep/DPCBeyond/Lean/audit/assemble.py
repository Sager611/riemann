
import re,os,subprocess
R='Riemann/Deep/DPCBeyond/'
first=open('.scratch/dpcblean/first.txt').read().strip()
log=open('.scratch/dpcblean/verify.log').read()
pk={}
for m in re.finditer(r'^(\d\d:\d\d:\d\d) (\w+) rc=(\d+) (\d+)s footprint=(\d+)MiB memcap: peak (\d+) MB',log,re.M):
    pk[m.group(2)]=(m.group(3),m.group(4),m.group(5),m.group(6))
desc={
'ParamsB':'(a) one_lt_of, step0B_nonempty, step0B_signs, step8B_empty, step0B_iff, exists_q_s, exists_u_w, order_of_choices, conrey_threshold, conrey_sharp',
'Rpow':'rp (x^e) and lgex, lg_lt, lg_le, lg_nonneg, lg_nonpos, rp_of_pos, rp_of_nonpos, rp_pos, rp_nonneg, rp_one, rp_zero, rp_add, rp_mul, rp_rp, rp_neg, rp_le_rp, rp_lt_rp, rp_le_rp_exp, rp_le_max',
'Step7B':'(e) step7B_meas, step7B',
'HairpinB':'(d) DensityK, Step6ClassesB, class_count, class_bound',
'SumB':'(d) step6B_of_classes, PhiB_small',
'DLamB':'(f) UFrom, dpos_of_from, nToInf_of_from, dlam_pos_of_from',
'TheoremUB':'TheoremUB, DataB, DataB.base, HolderB, MomentM, Step5cB, InputsB, DPCInputsB, cover_of, moments_of, split_c, theoremUB',
'FinalB':'(f) D_pos_of_kappa, theoremUB_conrey, D_pos_seven_elevenths',
'LemmaZ':'(b) ex_sum, sum_mul_le, prod_ex_le, prod_ex_le_sup, rp_sup_le_sum, mul_prod_eq, tuple_one, tuple_sum, tuple_bound, packet_case1, packet_case2, packet_nofree, pos_part_shift, pos_part_split, term_bound, mono_bound, lemmaZ_sum',
'LemmaM':'(c) rp_ex, M2_bound, geom_mul, geom_le, M3_series, M3_exists_w',
'LemmaM4':'(c) M4_count, M4_moment, key_ineq, inv_sq_sum, kernel_sum',
'LemmaM3':'(c) off_bound, momentM_of, exists_r_w, M3Inputs, momentM_of_inputs',
'FineB':'InputsFineB, DPCInputsFineB, inputsB_of_fine, theoremUB_fine, D_pos_seven_elevenths_fine',
'CheckB':'#print axioms (62 theorems)'}
order=['ParamsB','Rpow','Step7B','HairpinB','SumB','DLamB','TheoremUB','FinalB','LemmaZ','LemmaM','LemmaM4','LemmaM3','FineB','CheckB']
rows=[]; tot=0
for m in order:
    n=sum(1 for _ in open(R+'Lean/'+m+'.lean')); tot+=n
    rc,t,fp,mc=pk.get(m,('?','?','?','?'))
    rows.append(f'| {m} | {desc[m]} | {n} | {fp} MiB | {mc} MB | {t} s |')
deps=['FinSum','Model','Liminf','Params','Gaps','Tail','Phi','Step7','TheoremU','DLam','Hairpin','Fine']
dl=', '.join(f'{d} {pk[d][2]}' for d in deps if d in pk)
allfp=[int(pk[m][2]) for m in order+deps if m in pk]
allrc=[pk[m][0] for m in order+deps if m in pk]
ax=open('.scratch/dpcblean/log_CheckB.txt').read()
nax=len(re.findall(r"depends on axioms: \[propext, Classical.choice, Quot.sound\]",ax))
nbad=len(re.findall(r"depends on axioms",ax))-nax
mainst=open('.scratch/dpcblean/main_statements.txt').read().rstrip()
prog=open(R+'LEAN.md').read()
prog=prog[prog.index('## Progress'):].rstrip()
parts={k:open(f'.scratch/dpcblean/part{k}.md').read().rstrip() for k in 'BCEFT'}
hdr=f'''{first}

# DPCBeyond (Theorem U′) in Lean 4 + Mathlib (dpcblean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Riemann/Deep/DPCBeyond/Lean/, {len(order)} files, {tot} lines, namespace DPCBLean; not imported by Riemann.lean. Imports dpclean's Theorem U modules unchanged (Riemann/Deep/DPCUncond/Lean: Tail's ExpLog, Phi, Hairpin's hairpin_sum, Step7's step7_bad, Gaps' gapSet/offSet/Reflect/card_off_le, TheoremU's Data/Step3/MainTerm/Step1Pos/Covering, DLam's RvMUp/ScaleHalf/scaleToInf_of/scaleHalf_of/gapSet_half_sub, Fine's deepSet/classSet) and through them pplean's zero model and liminf D (Riemann/Deep/ExplorePointProc/Lean). No Mathlib zeta, Fourier, LSeries, complex analysis, measure theory, Real.log, Real.rpow or π (290 MB cap): integrals and norms are abstract real functions of T, exp/log abstract (ExpLog), x^e = rp ex lg x e, every analytic input a named hypothesis.

| file | content | lines | footprint (time -l) | memcap reading | time |
|---|---|---|---|---|---|
''' + "\n".join(rows) + f'''

Imported modules in the same final run (footprint MiB): {dl}. Final run: all {len(allrc)} compiles rc = {'0' if all(r=='0' for r in allrc) else 'NOT ALL 0'}, peak footprint {min(allfp)}–{max(allfp)} MiB, no watchdog kill.

## Main statements (verbatim; proofs omitted)

    -- pplean: Z.γ n (n ≥ 1, with multiplicity), Z.N T, Z.sc T = log T/(2π); dpclean: gapSet Z φ T = 𝒩(T), Data, Step3,
    -- MainTerm, Step1Pos, Covering, Reflect, offSet, classSet, deepSet, ExpLog, Phi (Riemann/Deep/DPCUncond/LEAN.md)
''' + "\n".join("    "+l if l.strip() else "" for l in mainst.split("\n")) + '''
'''
build=f'''## Build

Each module separately, one at a time, lean directly (no lake process), oleans to .scratch/dpcblean/olean (build output, deleted at the end):

    rm .lake && bash .scratch/memcap.sh 290 3600 lake exe cache get        # 14:59–15:01 UTC: peak 243 MB, rc = 0
    cp Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt .scratch/dpcblean/leanpath.txt     # lake's LEAN_PATH
    # per module m (dir d): pplean's FinSum Model Liminf, dpclean's Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine, then
    # ParamsB Rpow Step7B HairpinB SumB DLamB TheoremUB FinalB LemmaZ LemmaM LemmaM4 LemmaM3 FineB CheckB (d = Riemann/Deep/DPCBeyond/Lean):
    LEAN_PATH=.scratch/dpcblean/olean:$(cat .scratch/dpcblean/leanpath.txt) bash .scratch/memcap.sh 290 1800 \\
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/dpcblean/olean/$d/$m.olean $d/$m.lean

Scripts and logs: .scratch/dpcblean/verify.sh (clean rebuild of all {len(allrc)} modules), mk.sh (one module), verify.log, verify1.log (first full run), peaks.log (every compile of the session), log_<module>.txt, review_prompt.txt, review_out.md, statements.txt. rg over Lean/*.lean: no sorry, admit, axiom, native_decide, decide, skipKernelTC. CheckB.lean: {nax} #print axioms lines, all [propext, Classical.choice, Quot.sound]{'' if nbad==0 else f'; {nbad} OTHER'}:

    'DPCBLean.theoremUB' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for step0B_nonempty, step0B_signs, step8B_empty, step0B_iff, exists_q_s, exists_u_w, order_of_choices, conrey_threshold,
     conrey_sharp, rp_add, rp_mul, rp_rp, rp_le_rp, rp_lt_rp, rp_le_rp_exp, rp_le_max, step7B_meas, step7B, class_count, class_bound,
     step6B_of_classes, PhiB_small, dpos_of_from, nToInf_of_from, dlam_pos_of_from, cover_of, moments_of, split_c, D_pos_of_kappa,
     theoremUB_conrey, D_pos_seven_elevenths, ex_sum, sum_mul_le, prod_ex_le, prod_ex_le_sup, rp_sup_le_sum, tuple_sum, tuple_bound,
     packet_case1, packet_case2, packet_nofree, pos_part_shift, pos_part_split, term_bound, mono_bound, lemmaZ_sum, M2_bound, geom_le,
     M3_series, M3_exists_w, M4_count, M4_moment, inv_sq_sum, kernel_sum, off_bound, momentM_of, exists_r_w, momentM_of_inputs,
     inputsB_of_fine, theoremUB_fine, D_pos_seven_elevenths_fine)
'''
fz=open('.scratch/dpcblean/fz.txt').read().strip().replace('NLINES',str(tot)).replace('MAXFP',str(max(allfp)))
out=hdr+"\n"+parts['B']+"\n\n"+build+"\n"+parts['E']+"\n\n"+parts['C']+"\n\n"+parts['T']+"\n\n"+parts['F']+"\n\n"+prog+"\n\n## FORMALIZATIONS.md entry\n\n"+fz+"\n"
open(R+'LEAN.md','w').write(out)
print("written", out.count("\n"), "lines;", "rc all 0:", all(r=='0' for r in allrc), "maxfp", max(allfp), "nax", nax, "nbad", nbad, "missing", [m for m in order+deps if m not in pk])

