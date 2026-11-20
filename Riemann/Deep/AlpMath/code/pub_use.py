import re, collections
P = "/Users/adriansager/playground/riemann/.scratch/alpmath/pubuse.txt"
names = "HidaGal|SelfTw|IsNormChar|WeakBC|AdjointLift|SymCube|CyclicBaseChangeC|QuadInduction|HilbertEisOne|HilbertOne|KPlusExists|KAcyclicHomology|KMilnor|KHurewicz|KStability|IsSel|coefField|IsLat|FerreroWashington|KisinCrFull|KFinGen|KBorelRank|KFiniteField|KLocalisation|KDwyerFriedlander|DescentUnitary|TypesNZero|SignCsd|CyclicBaseChange|JacquetLanglands|DefQuatExists|FormsFiniteReal|CuspSemisimpleReal|DeRhamCharAlg|ClozelCohC|ClozelCohR|ClozelConjC|ClozelConjR|BorelInjC|BorelInjR|BorelSerreFL|GaloisRepReal|StrongMultOneReal|GelfandKazhdanReal|AdmissibleReal|IsQuadCharOf|QuadBaseChange|ThorneAdjoint|NewformRep|RepNewform|IsGalRepOfR|IsBC".split("|")
files = collections.defaultdict(set); proved = collections.defaultdict(list); assumed = collections.defaultdict(list)
pubdir = "Development/AutomorphicArithmetic/AutomorphicVarieties/PublishedTheorems/"
for line in open(P, encoding="utf-8", errors="replace"):
    m = re.match(r"^\./([^:]+):(\d+):(.*)$", line.rstrip("\n"))
    if not m: continue
    f, ln, txt = m.group(1), m.group(2), m.group(3)
    for n in names:
        if re.search(r"(?<![A-Za-z0-9_])(Published\.)?" + n + r"(?![A-Za-z0-9_])", txt):
            files[n].add(f)
            if re.match(r"\s*(theorem|lemma|def|abbrev)\s+\S+\s*:\s*[^=]*" + n, txt) and not f.startswith(pubdir):
                proved[n].append(f + ":" + ln)
            if re.search(r"\(\w+ : [^()]*" + n + r"\b[^()]*\)", txt):
                assumed[n].append(f + ":" + ln)
for n in names:
    fs = [x for x in files[n] if not x.startswith(pubdir)]
    print(n, "files:", len(fs), "proved-at:", proved[n][:2], "assumed-at:", len(assumed[n]), assumed[n][:1])
