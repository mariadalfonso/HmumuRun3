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

    # Common groups
    hmumu   = [f"{i}" for i in range(10, 16)]       # 10–15
    hzgamma = [f"{i}" for i in range(20, 25)]       # 20–24
    vv      = [f"{i}" for i in range(201, 206)] + [f"{i}" for i in range(211, 215)]  # 201–205, 211–214
    ttV     = [f"{i}" for i in range(221, 226)]

    # Add Hmumu
    for tag in hmumu:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # extra DY
    safe_add_tree(mytree, f"{directory}snapshot_mc_109{year}_{category}.root")  # DY 105-160
    safe_add_tree(mytree, f"{directory}snapshot_mc_108{year}_{category}.root")  # powheg 50-120
    safe_add_tree(mytree, f"{directory}snapshot_mc_110{year}_{category}.root")  # powheg 120-200
    if category in ["VHcat"]:
        safe_add_tree(mytree, f"{directory}snapshot_mc_111{year}_{category}.root")  # DY 0j - 2j
        safe_add_tree(mytree, f"{directory}snapshot_mc_112{year}_{category}.root")  # DY 0j - 2j
        safe_add_tree(mytree, f"{directory}snapshot_mc_113{year}_{category}.root")  # DY 0j - 2j
        safe_add_tree(mytree, f"{directory}snapshot_mc_114{year}_{category}.root")  # PT binned
        safe_add_tree(mytree, f"{directory}snapshot_mc_115{year}_{category}.root")  # PT binned
        safe_add_tree(mytree, f"{directory}snapshot_mc_116{year}_{category}.root")  # PT binned
        safe_add_tree(mytree, f"{directory}snapshot_mc_117{year}_{category}.root")  # PT binned

    # Special cases
    if year in {"_12022", "_22022", "_12023", "_22023"}:
        safe_add_tree(mytree, f"{directory}snapshot_mc_100{year}_{category}.root")  # DY
    elif year == "_2024":
        safe_add_tree(mytree, f"{directory}snapshot_mc_103{year}_{category}.root")  # DY

    if category in ["VLcat","TTLcat"]:
        # Add VV
        for tag in vv:
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Always add
    if category in ["TTHcat","TTLcat"]:
        for tag in ["102"]:  # TT2L
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    return mytree

    # Add year-specific samples
    for tag in data_samples.get(year):
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add Hmumu
    for tag in hmumu:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add H→Zγ
    for tag in hzgamma:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Special cases
    if year in {"_12022", "_22022", "_12023", "_22023"}:
        safe_add_tree(mytree, f"{directory}snapshot_mc_25{year}_{category}.root") #Zg
        safe_add_tree(mytree, f"{directory}snapshot_mc_100{year}_{category}.root")  # DY
    elif year == "_2024":
        safe_add_tree(mytree, f"{directory}snapshot_mc_26{year}_{category}.root") #Zg
        safe_add_tree(mytree, f"{directory}snapshot_mc_103{year}_{category}.root")  # DY
        safe_add_tree(mytree, f"{directory}snapshot_mc_104{year}_{category}.root")  # DY

    safe_add_tree(mytree, f"{directory}snapshot_mc_101{year}_{category}.root")  # DY-EWK

    # Always add
    for tag in ["102"]:  # TT2L
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add VV
    for tag in vv:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    if category in ["VLcat", "TTLcat", "TTHcat", "VHcat", "ggHcat"]:
        other = ["105", "106", "107"]  # TW, TT1L
        for tag in other:
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

        # Add TTV
        for tag in ttV:
            if year in {"_12022", "_22022", "_12023", "_22023"}:
                safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")
            if year == "_2024":
                safe_add_tree(mytree, f"{directory}snapshot_mc_222{year}_{category}.root")
                safe_add_tree(mytree, f"{directory}snapshot_mc_223{year}_{category}.root")
                safe_add_tree(mytree, f"{directory}snapshot_mc_224{year}_{category}.root")

    return mytree
