import ROOT
import math
import itertools

ROOT.EnableImplicitMT()
from LoadTree import loadTree

# ============================
# CONFIGURATION
# ============================

directory = '/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

category = 'ggHcat'
year = '_2024'        

signal_files = {
    "ggHcat":  ["11"],             # ggH
    "VBFcat":  ["10"],             # VBF
    "VHcat":   ["12", "13", "14"], # VH
    "VLcat":   ["12", "13", "14"], # VL
    "Zinvcat": ["14"],             # VL
    "TTHcat":  ["15"],             # TTH
    "TTLcat":  ["15"],             # TTL
}

mytree = ROOT.TChain('events')
mytree = loadTree(mytree, directory, category, year )

MASS_CENTER = 125.0
MASS_WINDOW = 25.0

MVAdiscr_min = 0.0
MVAdiscr_max = 1.0

# minimum background events per bin
#MIN_BKG = 100   # corresponds to 10% stat uncertainty
MIN_BKG = 20  ## from Zgamma

# max categories to try
MAX_CAT = 3

# number of scan points for MVAdiscr
SCAN_POINTS = 50

#SIGNAL_WEIGHT_BRANCH = "w*lumiIntegrated"
#SIGNAL_WEIGHT_BRANCH = "w"
WEIGHT_BRANCH = "w_allSF"

# ============================
# SIGNIFICANCE FUNCTION
# ============================

def compute_significance(S, B):

    total = 0.0

    for s, b in zip(S, B):

        if b > 0:
            total += (s*s)/b

    return math.sqrt(total)

# ============================
# LOAD DATAFRAMES
# ============================

df_sig = ROOT.RDataFrame(mytree).Filter(f"mc==10 or mc==11 or mc==12 or mc==13 or mc==14 or mc==15")
df_bkg = ROOT.RDataFrame(mytree).Filter("mc<0")

sig = df_sig.Filter(
    f"abs(HiggsCandCorrMass - {MASS_CENTER}) <= {MASS_WINDOW}"
)

bkg = df_bkg.Filter(
    f"abs(HiggsCandCorrMass - {MASS_CENTER}) <= {MASS_WINDOW}"
)

print("Total signal weighted:",
      sig.Sum(WEIGHT_BRANCH).GetValue())

print("Total signal raw:",
      sig.Count().GetValue())

print("Total background:",
      bkg.Count().GetValue())

print("Total background weighted:",
      bkg.Sum(WEIGHT_BRANCH).GetValue())

# ============================
# PRECOMPUTE CUMULATIVE COUNTS
# ============================

grid = [
    MVAdiscr_min + i*(MVAdiscr_max-MVAdiscr_min)/SCAN_POINTS
    for i in range(SCAN_POINTS+1)
]

sig_counts = []
bkg_counts = []

for edge in grid:

    s = sig.Filter(f"discrMVA0 >= {edge}").Sum(WEIGHT_BRANCH).GetValue()
    b = bkg.Filter(f"discrMVA0 >= {edge}").Sum(WEIGHT_BRANCH).GetValue()

    sig_counts.append(s)
    bkg_counts.append(b)

# helper function
def count_between(i, j):

    s = sig_counts[i] - sig_counts[j]
    b = bkg_counts[i] - bkg_counts[j]

    return s, b

# ============================
# OPTIMIZATION LOOP
# ============================

best_Z = 0
best_edges = None

print("\nStarting optimization...\n")

for Ncat in range(1, MAX_CAT+1):

    print(f"Trying {Ncat} categories...")

    for indices in itertools.combinations(range(1, SCAN_POINTS), Ncat-1):

        edges_idx = (0,) + indices + (SCAN_POINTS,)

        S = []
        B = []

        valid = True

        for i in range(len(edges_idx)-1):

            s, b = count_between(edges_idx[i], edges_idx[i+1])

            if b < MIN_BKG:
                valid = False
                break

            S.append(s)
            B.append(b)

        if not valid:
            continue

        Z = compute_significance(S, B)

        if Z > best_Z:

            best_Z = Z
            best_edges = [grid[i] for i in edges_idx]

# ============================
# PRINT RESULT
# ============================

print("\n====================================")
print("BEST RESULT")
print("====================================")

print(f"Best significance: {best_Z:.4f}")
print(f"Number of categories: {len(best_edges)-1}")

print("\nEdges:")

for e in best_edges:
    print(f"{e:.4f}")

print("\nDetailed bins:\n")
print(category)

for i in range(len(best_edges)-1):

    low = best_edges[i]
    high = best_edges[i+1]

    s = sig.Filter(
        f"discrMVA0 >= {low} && discrMVA0 < {high}"
    ).Sum(WEIGHT_BRANCH).GetValue()

    b = bkg.Filter(
        f"discrMVA0 >= {low} && discrMVA0 < {high}"
    ).Sum(WEIGHT_BRANCH).GetValue()

    print(f"[{low:.4f}, {high:.4f}]  S={s}  B={b}  S/sqrt(B)={s/math.sqrt(b):.3f}")
