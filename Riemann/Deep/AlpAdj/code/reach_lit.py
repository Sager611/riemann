#!/usr/bin/env python3
# BFS over alpdep's dependency edges from given declaration names; report reachable literature Props
# (field types of LitLoc/LitAut/LitCoh/LitCM, Hypotheses', HypothesesH') and Hida/R=T-looking names.
# usage: python3 reach_lit.py <name1> [<name2> ...]   (memory: numpy arrays + one set of visited ids)
import sys, re, numpy as np
OUT = '/Users/adriansager/playground/riemann/.scratch/alpdep'
eoff = np.load(OUT + '/eoff.npy', mmap_mode='r'); elen = np.load(OUT + '/elen.npy', mmap_mode='r')
store = np.memmap(OUT + '/estore.bin', dtype=np.int32, mode='r')
dist = np.load(OUT + '/dist.npy', mmap_mode='r')
targets = set(sys.argv[1:])
LIT = set('''BernsteinZelevinskyFiniteness BernsteinUnitaryProductIrreducible GeometricLemma ZelevinskyPolynomialRing
TotallyDisjointProductIrreducible MoyPrasadJacquetIso MinguezThetaGL TemperedProductOfSquareIntegrable EssSquareIntegrableSegments
SquareIntegrableSegment UnitaryDualGL GenericUnlinkedGL LanglandsMultisegments LanglandsExhaustionGL StandardCosocleGL LocalLanglandsExists
LocalConverseGL DepthZeroCuspidalLL DepthZeroInduction DepthZeroTypeInertia BrauerDegreeZero ChevalleyUnits GreenRegularCharacters
CuspRapidDecay HarishChandraFinite AutomorphicFormsFiniteGL JacquetShalika CuspidalSpectrumGL DiscreteSpectrumResidual BaseChangeIsobaric
GelfandKazhdanDual GaloisRepOfPolarizedPair CaraianiLocalGlobalDual CaraianiTempered ClozelPurity ClozelCohomological GlobalArtinReciprocity
RayClassFieldSt HeckeLAdicCharacterAbove AlgebraicCharacterWD DeRhamPotSemistable PoitouTateNormalised LocalFiniteness LocalTateDuality
LocalEulerChar GlobalFiniteness GlobalEulerChar UnramifiedOrthogonality GreenbergWiles LocalFinitenessK LocalTateDualityK LocalEulerCharK
BlochKatoLocal HilbertBlumenthal GrothendieckDeformation AbelianDeRham AbelianDualCrystal PeriodComparison GaussManinPeriods
SerreTateGoodReductionSt EllipticMainTheoremCM HilbertClassFieldJ SiegelKlingenDeligneRibet PartialZetaFunctionalEquation
HeckeTateFunctionalEquation ClozelCohC ClozelConjC BorelInjC GaloisRepExistsTotallyRealGL2 StrongMultOneReal GelfandKazhdanReal
AdmissibleReal QuadBC0 NewformRep DescWeak CycC Desc2 DefQuatExists CuspSemisimpleReal LiteratureRankNHypothesis CohP'''.split())
ids = {}
names = {}
with open(OUT + '/names.tsv', encoding='utf-8') as f:
    for line in f:
        p = line.split('\t', 2)
        if p[1] in targets: ids[p[1]] = int(p[0])
print('targets', ids)
N = len(eoff); seen = np.zeros(N, dtype=np.uint8); stack = list(ids.values())
for i in stack: seen[i] = 1
cnt = 0; unproc = 0
while stack:
    u = stack.pop(); cnt += 1
    o = int(eoff[u]); l_ = int(elen[u])
    if o < 0: unproc += 1; continue
    for v in store[o:o + l_].tolist():
        if 0 <= v < N and not seen[v]: seen[v] = 1; stack.append(v)
vis = set(np.nonzero(seen)[0].tolist())
print('reached', len(vis), 'nodes; processed', cnt, 'unprocessed leaves', unproc)
hits = []; hida = []
with open(OUT + '/names.tsv', encoding='utf-8') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        i = int(p[0])
        if i in vis:
            last = p[1].rsplit('.', 1)[-1]
            if last in LIT: hits.append((last, p[1], p[2]))
            elif re.search(r'(Hida|Taylor|Wiles|Patch|Kisin|Thorne|Allen|Newton|Geraghty|Calegari|Skinner|Eisenstein|Bellaiche|Chenevier|RisT|RequalsT|RT)', p[1]) and p[2] in ('D', 'S'):
                hida.append(p[1])
print('literature Props reached:', len(hits))
for h in sorted(hits): print('  ', h[0], h[1], h[2])
print('Hida/patching-named defs reached:', len(hida))
for h in sorted(hida)[:60]: print('  ', h)
