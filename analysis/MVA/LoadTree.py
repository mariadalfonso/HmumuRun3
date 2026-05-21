import ROOT
import os

def safe_add_tree(chain, filepath, treename="events"):
    """Add file to TChain only if it exists and contains the desired TTree."""
    if not os.path.exists(filepath):
        print(f"⚠️  File not found: {filepath}")
        return False

    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        print(f"⚠️  Could not open file: {filepath}")
        return False

    tree = f.Get(treename)
    if not tree:
        print(f"⚠️  File has no TTree '{treename}': {filepath}")
        f.Close()
        return False

    # All good
    f.Close()
    chain.Add(filepath)
    return True

# Common groups
hmumu   = [f"{i}" for i in range(10, 16)]       # 10–15
hzgamma = [f"{i}" for i in range(20, 25)]       # 20–24
vv      = [f"{i}" for i in range(201, 207)] + [f"{i}" for i in range(207, 213)] + [f"{i}" for i in range(213, 217)]  # 201–206, 211–214
tt2l    = ["140", "141", "142"]
ttV     = [f"{i}" for i in range(221, 238)] + ["107","105","106"] #ttV, tt1l , singletop

def loadTree(mytree, directory_, category, year ):

    directory = directory_ + category + "/"

    # Year-specific samples
    data_samples = {
        "_12022": ["-11", "-13", "-14"],
        "_22022": ["-15", "-16", "-17"],
        "_12023": ["-23", "-24"],
        "_22023": ["-31", "-32"],
        "_2024" : [str(i) for i in range(-41, -55, -1)],  # -41 … -54
    }

    # Add Hmumu
    for tag in hmumu:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    ####
    #### DY
    ####

    if category in ["VHcat"]:

        if year in {"_12022", "_22022", "_12023", "_22023"}:
#            safe_add_tree(mytree, f"{directory}snapshot_mc_111{year}_{category}.root")  # DY 0j - 2j
#            safe_add_tree(mytree, f"{directory}snapshot_mc_112{year}_{category}.root")  # DY 0j - 2j
#            safe_add_tree(mytree, f"{directory}snapshot_mc_113{year}_{category}.root")  # DY 0j - 2j

            safe_add_tree(mytree, f"{directory}snapshot_mc_114{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_115{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_116{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_117{year}_{category}.root")  # PT binned
        elif year == "_2024":
#            safe_add_tree(mytree, f"{directory}snapshot_mc_119{year}_{category}.root")  # DY 0j - 2j
#            safe_add_tree(mytree, f"{directory}snapshot_mc_120{year}_{category}.root")  # DY 0j - 2j
#            safe_add_tree(mytree, f"{directory}snapshot_mc_121{year}_{category}.root")  # DY 0j - 2j

            safe_add_tree(mytree, f"{directory}snapshot_mc_122{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_123{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_124{year}_{category}.root")  # PT binned
            safe_add_tree(mytree, f"{directory}snapshot_mc_125{year}_{category}.root")  # PT binned

    if category in ["VBFcat", "ggHcat"]:

        # standard
        if year in {"_12022", "_22022", "_12023", "_22023"}:
            safe_add_tree(mytree, f"{directory}snapshot_mc_100{year}_{category}.root")  # DY
        elif year == "_2024":
            safe_add_tree(mytree, f"{directory}snapshot_mc_103{year}_{category}.root")  # DY
            safe_add_tree(mytree, f"{directory}snapshot_mc_104{year}_{category}.root")  # DY

        # DY EWK
        safe_add_tree(mytree, f"{directory}snapshot_mc_101{year}_{category}.root")
        safe_add_tree(mytree, f"{directory}snapshot_mc_99{year}_{category}.root")
        safe_add_tree(mytree, f"{directory}snapshot_mc_98{year}_{category}.root")

        # extra DY QCD
        safe_add_tree(mytree, f"{directory}snapshot_mc_109{year}_{category}.root")  # DY 105-160
        safe_add_tree(mytree, f"{directory}snapshot_mc_108{year}_{category}.root")  # powheg 50-120
        safe_add_tree(mytree, f"{directory}snapshot_mc_110{year}_{category}.root")  # powheg 120-200

    ####
    #### TOP
    ####

    if category in ["Zinvcat","TTHcat","TTLcat","VHcat"]:
        for tag in tt2l:  # TT2L
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    if category in ["TTLcat"]:
        for tag in ttV:  # 1L ttv
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    ####
    #### diboson
    ####

    if category in ["VLcat"]:
        # Add VV
        for tag in vv:
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")


    return mytree
