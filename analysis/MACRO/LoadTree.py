import ROOT
import os

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

    # Add year-specific samples
    for tag in data_samples.get(year):
        mytree.Add(f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add Hmumu
    for tag in hmumu:
        mytree.Add(f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add H→Zγ
    for tag in hzgamma:
        mytree.Add(f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Special cases
    if year in {"_12022", "_22022", "_12023", "_22023"}:
        mytree.Add(f"{directory}snapshot_mc_25{year}{category}.root")
        mytree.Add(f"{directory}snapshot_mc_100{year}{category}.root")  # DY
    elif year == "_2024":
        mytree.Add(f"{directory}snapshot_mc_26{year}{category}.root")
        mytree.Add(f"{directory}snapshot_mc_103{year}{category}.root")  # DY
        mytree.Add(f"{directory}snapshot_mc_104{year}{category}.root")  # DY

    # Always add
    common = ["102", "105", "106"]  # TT2L, TW2L, TW2L
    for tag in common:
        mytree.Add(f"{directory}snapshot_mc_{tag}{year}{category}.root")

    # Add VV
    for tag in vv:
        mytree.Add(f"{directory}snapshot_mc_{tag}{year}{category}.root")

    return mytree
