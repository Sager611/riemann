# orderslean STATE

Start: 2026-11-12 16:30 UTC

- reading inputs
- 16:50 TeX v1+v2 read (rules identical; v2 l.702 still open; v2 l.640-660 no c=3 patch in its shape). LEAN.md fidelity lines written.
- Plan: Defs (spec+checker iff) -> Ext (Steps, perfect_extend) -> Bulk -> Patch (Lin checker soundness) -> Search (layers, sound) -> Fail -> D_* per order (generated) -> Main -> Small -> Check.
- 17:01 Defs, Ext, Bulk compiled (Bulk 253 MB footprint: watch). Patch.lean: checker + identB/segB/pool lemmas compiling; next pool_nodup, patchOK, patch_sound.
- 17:10 Patch.lean (patch_sound, Lemma 3 generic) + Induct.lean (works_of) compiled. D_pxmd (order +x-/) compiled: Theorem 2 for N>=9, 28 s, 210 MB. Generating/compiling other 10 orders. Next: Fail (layers search), Main corollary, Small (exact sets), Check.
- 17:22 All 11 D_* compiled (max 260 MB). Search.lean (sound search), Fail.lean (13 bounds), Small1-6 (exact small sets; split, max 234 MB) compiled. Next Main (Thm1 + Corollary), Check, LEAN.md, reviewer.
- 17:34 Main + Check compiled; clean rebuild 26/26 rc=0, max 258 MB, 433 s. 82 axiom lines standard. Reviewer: FAITHFUL. Writing LEAN.md, then cleanup.
- 17:35 done: LEAN.md written (FORMALIZED), cleanup done (downloads, tests, oleans, __pycache__). .lake untouched.
