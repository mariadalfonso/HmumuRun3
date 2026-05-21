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
hzgamma = [f"{i}" for i in range(20, 27)]       # 20-26
hww     = [f"{i}" for i in range(30, 38)]       # 30–37
vv      = [f"{i}" for i in range(201, 207)] + [f"{i}" for i in range(213, 217)]  # 201–206, 213–216
ttV     = [f"{i}" for i in range(221, 238)]
tt2l    = ["140"]
top     = ["105", "106", "107"]
dyewk   = ["101"]
dy_2223 = ["100"]
dy_24     = ["103","104"]
dy_pt2223 = ["114","115","116","117"]
dy_pt24 = ["122","123","124","125"]
dy_minllo = ["126","127"]
dy_j    = ["111","112","113"]

def loadTree(mytree, directory_, category, year ):

    directory = directory_ + category + "/"

    # Year-specific samples
    data_samples = {
        "_12022": ["-11", "-13", "-14"],
        "_22022": ["-15", "-16", "-17"],
        "_12023": ["-23", "-24"],
        "_22023": ["-31", "-32"],
        "_2024" : [str(i) for i in range(-41, -55, -1)],  # -41 … -54
        "_12024" : [str(i) for i in range(-41, -47, -1)],  # -41 … -46
        "_22024" : [str(i) for i in range(-47, -55, -1)],  # -47 … -54
        "_2025" : [str(i) for i in range(-61, -73, -1)],  # -61 … -72
    }

    print(data_samples[year])
    if year == '_12024': yearX = '_2024'
    elif year == '_22024': yearX = '_2024'
    elif year == '_22025': yearX = '_2025'
    else: yearX = year
    # Add year-specific samples
    for tag in data_samples.get(year):
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{yearX}_{category}.root")

    if year == '_12024': year = '_2024'
    if year == '_22024': year = '_2024'
    if year == '_2025': year = '_2024'

    # Add Hmumu
    for tag in hmumu:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add H→Zγ
    for tag in hzgamma:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add H→WW
#    if category in ["TTLcat"] and year in {"_2024"}:
#        for tag in hww:
#            safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")  # DY 0j

    # Add W+jets
#    if category in ["VHcat"]:
#        for tag in hww:
#            safe_add_tree(mytree, f"{directory}snapshot_mc_129{year}_{category}.root")

    # Special cases
    if year in {"_12022", "_22022", "_12023", "_22023"}:
        if category in ["VHcat"]: # later more stat , but normalization need to be fixed
            for tag in dy_pt2223: safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")
        else:
            for tag in dy_2223: safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")  # DY
    elif year in {"_2024"}:
        if category in ["VHcat"]:
            for tag in dy_pt24: safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")
        else:
            for tag in dy_24: safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    for tag in dy_minllo:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")  # DY-EWK

    # Add DYEWK
    for tag in dyewk:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")  # DY-EWK

    # Add TT2L
    for tag in tt2l:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add VV
    for tag in vv:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add TW, TT1L
    for tag in top:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    # Add TTV
    for tag in ttV:
        safe_add_tree(mytree, f"{directory}snapshot_mc_{tag}{year}_{category}.root")

    return mytree
