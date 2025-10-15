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

def loadTree(mytree, directory, category, year ):

    # Year-specific samples
    data_samples = {
        "_12022": ["-11", "-12", "-13", "-14"],
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

    # Add year-specific samples
    for tag in data_samples.get(year):
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add Hmumu
    for tag in hmumu:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add H→Zγ
    for tag in hzgamma:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Special cases
    if year in {"_12022", "_22022", "_12023", "_22023"}:
        safe_add_tree(mytree, f"{directory}snapshot_mc_25{year}{category}.root") #Zg
        safe_add_tree(mytree, f"{directory}snapshot_mc_100{year}{category}.root")  # DY
    elif year == "_2024":
        safe_add_tree(mytree, f"{directory}snapshot_mc_26{year}{category}.root") #Zg
        safe_add_tree(mytree, f"{directory}snapshot_mc_103{year}{category}.root")  # DY
        safe_add_tree(mytree, f"{directory}snapshot_mc_104{year}{category}.root")  # DY

    safe_add_tree(mytree, f"{directory}snapshot_mc_101{year}{category}.root")  # DY-EWK

    # Always add
    for tag in ["102"]:  # TT2L
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add VV
    for tag in vv:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

    if category in ["_VLcat", "_TTLcat", "_TTHcat", "_VHcat", "_ggHcat"]:
        other = ["105", "106", "107"]  # TW, TT1L
        for tag in other:
            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")

        # Add TTV
        for tag in ttV:
            if year in {"_12022", "_22022", "_12023", "_22023"}:
                safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}{category}.root")
            if year == "_2024":
                safe_add_tree(mytree, f"{directory}snapshot_mc_222{year}{category}.root")
                safe_add_tree(mytree, f"{directory}snapshot_mc_223{year}{category}.root")
                safe_add_tree(mytree, f"{directory}snapshot_mc_224{year}{category}.root")

    return mytree
