# RCAudit checks

- l12.py: Lemma 1 bound and exact sharp constant of Lemma 2 (uv run --no-project --with numpy python3 l12.py out.json 10000 all|big); sum12.py -> l12_summary.out (raw JSON 4.7 MB not kept).
- l4g.gp + mkl4.py: Lemma 4 second moment via PARI 2.19 (gp -q --default parisizemax=150000000 l4run.gp); raw l4g_out.txt; sum4.py -> l4_summary.out.
- epstein.py: exact X=1, u=0 second moment by Chowla-Selberg, fundamental D to 1e9; ep_big*.json; sumep.py -> epstein_summary.out.
- lit_scan1.out, lit_scan2.out: arXiv Oct-Nov 2026 title scans and Zenodo queries (.scratch/lst.py).
