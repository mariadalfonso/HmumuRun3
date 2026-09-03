import ROOT
import glob
from subprocess import check_output
import sys
import json
import os

def loadUserCode():
    print('loadUserCode()')
    ROOT.gSystem.AddDynamicPath("./.")
    ROOT.gROOT.ProcessLine(".include ./.")
    ROOT.gInterpreter.AddIncludePath("./.")
    ROOT.gInterpreter.ProcessLine('#include "./config/functions.h"')

def loadtmvahelper():
    print('loadtmvahelper()')
    ROOT.gInterpreter.ProcessLine('#include "./config/tmva_helper_xml.h"')
    ROOT.gInterpreter.ProcessLine('#include "./config/tmva_helper_xgb.h"')

# lumis with golden json
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#ReReco_ERAS_A_B_C_D_E
lumis={
    '12022':7.99, # C-D
    '22022':26.68, # E, F, G
    '12023':17.96, # C
    '22023':9.68, # D
    '2024':109.82, # C-I
    '2025':110.59, # C-G
    '2026':25.31 # C, B
}

btagPNetBM={
    '12022':0.245,
    '22022':0.2605,
    '12023':0.1917,
    '22023':0.1919,
    '2024':0.245, # found only UparT
    '2025':0.245, # copy the 2024 value
    '2026':0.245, # copy the 2024 value
}

btagPNetBL={
    '12022':0.047,
    '22022':0.0499,
    '12023':0.0358,
    '22023':0.0359,
    '2024':0.047, # found only UparT https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer24/#general-remarks
    '2025':0.047, # copy the 2024 value
    '2026':0.047, # copy the 2024 value
}

btagParTL={
    '12022':0.0849,
    '22022':0.0897,
    '12023':0.0681,
    '22023':0.0683,
    '2024':0.0246, # found only UparT
    '2025':0.0246, # copy the 2024 value
    '2026':0.0246, # copy the 2024 value
}

btagParTM={
    '12022':0.4319,
    '22022':0.451,
    '12023':0.3487,
    '22023':0.3494,
    '2024':0.1272, # found only UparT
    '2025':0.1272, # copy the 2024 value
    '2026':0.1272, # copy the 2024 value
}


xsecRun2={
    'ggH':48580,
    'VBFH':3781.7,
}

# Run3 xsection from https://xsecdb-xsdb-official.app.cern.ch/xsdb/
# HIGGS https://arxiv.org/pdf/2402.09955
# https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCTopWG
# top https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO#Updated_reference_cross_sections
# top BR https://pdg.lbl.gov/2017/reviews/rpp2017-rev-top-quark.pdf
# SM https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV


# from Marini 31 july
#kfactor = 2094 *3 / (5378 + 1017 + 385.5 ) = 0.9264803481 (MATRIX MLL>50)

kfactor={
    'Z':0.9265,
}

def SwitchSample(thisdict,argument):

    return thisdict.get(argument, "BKGdefault, xsecDefault")

def computeWeigths(rdf,xsec):

    genEventSumWeight = rdf.Sum("genEventSumw").GetValue()
    genEventSumNoWeight = rdf.Sum("genEventCount").GetValue()

    print('genEventSumWeight',genEventSumWeight)
    print('genEventSumNoWeight',genEventSumNoWeight)

    weight = xsec / genEventSumWeight

    lumiEq = (genEventSumNoWeight / xsec)
    print("lumi equivalent fb %s" %lumiEq)

    # neglecting the lumi
    #    weight = (1./genEventSumWeight)

    print('weight',weight)

    return weight


def amendJsonFile(file_path):
    import json
    import gzip

#    with gzip.open(file_path, "rt") as f:
#        data = json.load(f)
    with open(file_path, "rt") as f:
        data = json.load(f)

    # Your new correction
    new_corr = {
        "name": "RandomSmearing",
        "version": 1,
        "inputs": [
            {"name": "evtNr", "type": "int"},
            {"name": "lumiNr", "type": "int"},
            {"name": "phi", "type": "real"}
        ],
        "output": {
            "name": "random_number",
            "type": "real"
        },
        "data": {
            "nodetype": "hashprng",
            "inputs": ["evtNr", "lumiNr", "phi"],
            "distribution": "stdflat"
        }
    }

    data["corrections"].append(new_corr)

    # Convert back to JSON string
    json_str = json.dumps(data)

    escaped = json_str.replace('"', '\\"')

    ROOT.gROOT.ProcessLine(
        f'auto cset = correction::CorrectionSet::from_string("{escaped}");'
    )

def loadCorrectionSet(year):
    print('loadCorrectionSet()')
    import correctionlib
    correctionlib.register_pyroot_binding()
    ROOT.gInterpreter.Declare('#include "./config/sfCorrLib.h"')

    ROOT.gInterpreter.ProcessLine('auto corr_sf = MyCorrections(%d);' % (year))

    subDirName = ""
    '''
    if year == 12022: subDirName = "Run3-22CDSep23-Summer22-NanoAODv12"
    if year == 22022: subDirName = "Run3-22EFGSep23-Summer22EE-NanoAODv12"
    if year == 12023: subDirName = "Run3-23CSep23-Summer23-NanoAODv12"
    if year == 22023: subDirName = "Run3-23DSep23-Summer23BPix-NanoAODv12"
    if year == 2024: subDirName = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15"
    if year == 2025: subDirName = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15" # temporary
    '''
    if year == 12022: subDirName = "2022_Summer22"
    if year == 22022: subDirName = "2022_Summer22EE"
    if year == 12023: subDirName = "2023_Summer23"
    if year == 22023: subDirName = "2023_Summer23BPix"
    if year == 2024: subDirName = "2024"
    if year == 2025: subDirName = "2025"
    if year == 2026: subDirName = "2025"

    if(year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024 or year==2025 or year==2026):
        print('loadMuonScale()')
#        amendJsonFile("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/"+subDirName+"/latest/muon_scalesmearing.json.gz")
        amendJsonFile("config/POG/MUO/JSON_VXBS/"+subDirName+"/schemaV2.json")

#        ROOT.gROOT.ProcessLine(
#            f'auto cset = correction::CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/'+subDirName+'/latest/muon_scalesmearing.json.gz");'
#        )

        ROOT.gInterpreter.ProcessLine('#include "./config/MuonScaRe.cc"')
        ROOT.gInterpreter.Declare('#include "./config/functionsMuCorr.h"')

    print('loadHiggsNNLOPS()')
    ROOT.gROOT.ProcessLine(
        f'auto csetH = correction::CorrectionSet::from_file("./config/THeory/NNLOPS_reweight_13p6TeV.json");'
    )

    print('loadDYturbo()')
    ROOT.gROOT.ProcessLine(
        f'auto csetZ1 = correction::CorrectionSet::from_file("./config/THeory/qt_norm_reweight_ZR.json");'
    )
    ROOT.gROOT.ProcessLine(
        f'auto csetZ2 = correction::CorrectionSet::from_file("./config/THeory/qt_norm_reweight_SB.json");'
    )

    ROOT.gInterpreter.Declare('#include "./config/functionsObjCor.h"')

def loadJSON(fIn):

    if not os.path.isfile(fIn):
        print("JSON file %s does not exist" % fIn)
        return

    if not hasattr(ROOT, "jsonMap"):
        print("jsonMap not found in ROOT dict")
        return

    info = json.load(open(fIn))
    print("JSON file %s loaded" % fIn)
    for k,v in info.items():

        vec = ROOT.std.vector["std::pair<unsigned int, unsigned int>"]()
        for combo in v:
            pair = ROOT.std.pair["unsigned int", "unsigned int"](*[int(c) for c in combo])
            vec.push_back(pair)
            ROOT.jsonMap[int(k)] = vec

def readDataQuality(year):
    print("HELLO readDataQuality", year)
#    if year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024: dirJson = "/cvmfs/cms-griddata.cern.ch/cat/metadata/DC/Collisions22/latest/"

    dirJson = "./config"
    print('dirJson = ',dirJson)

    json_map = {
        "2018":  "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
        "2017":  "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
        "12016": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "22016": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "12022": "Cert_Collisions2022_355100_362760_Golden.json",
        "22022": "Cert_Collisions2022_355100_362760_Golden.json",
        "12023": "Cert_Collisions2023_366442_370790_Golden.json",
        "22023": "Cert_Collisions2023_366442_370790_Golden.json",
        "2024":  "Cert_Collisions2024_378981_386951_Golden.json",
        "2025":  "Cert_Collisions2025_391658_398903_Golden.json",
        "2026":  "Cert_Collisions2026_401624_403937_golden.json",
    }

    fname = json_map.get(str(year))
    if fname:
        loadJSON(f"{dirJson}/cert/{fname}")
    else:
        print(f"No JSON mapping found for year={year}")
