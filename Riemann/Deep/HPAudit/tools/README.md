Audit scripts (python3, stdlib only). scan.py OUT ROOT...: red-flag scan, comments/strings stripped -> OUT.hits.txt; show.py OUT.hits.txt cat1,~cat2 N: list hits (~ = grouped).
closure.py MOD1,MOD2 SRCDIR...: local import closure; decls.py ROOTS NAMES: namespace-aware declaration collisions; insts.py / insts2.py ROOT: global instances, those on C/R or bare variables.
cmpoai.py / cmpoai_diff.py: our QRH port vs the tracker's openai-baseline (env TRK_OAI = tracker clone .../Submissions/openai_baseline/OAI).
