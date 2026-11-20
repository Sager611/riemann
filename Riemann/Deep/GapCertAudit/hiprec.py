"""gapcertaudit: rerun certJ.py with EXTRA working bits (monkeypatch of gcert.set_prec), to separate Arb's arithmetic radius from
the Lemma T tail in the D = 60 test.  Run from a copy of ExploreShape.  usage: hiprec.py amp.json extra_bits"""
import sys, os
sys.path.insert(0, os.getcwd())
import certJ as CJ
extra = int(sys.argv[2])
orig = CJ.GC.set_prec
CJ.GC.set_prec = lambda bits: orig(bits + extra)
sys.argv = [sys.argv[0], sys.argv[1]]
CJ.main()
