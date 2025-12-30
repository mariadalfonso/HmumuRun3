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
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, # C
    '22023':9.451, # D
    '2024':108.95, # C-I
    '2025':115.54 # C-G
}

btagPNetBM={
    '12022':0.245,
    '22022':0.2605,
    '12023':0.1917,
    '22023':0.1919,
    '2024':0.245, # found only UparT
    '2025':0.245, # copy the 2024 value
}

btagPNetBL={
    '12022':0.047,
    '22022':0.0499,
    '12023':0.0358,
    '22023':0.0359,
    '2024':0.047, # found only UparT https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer24/#general-remarks
    '2025':0.047, # copy the 2024 value
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
xsecRun3={
    'ggH':52230, # 0.4 × σ(13 TeV) + 0.6 × σ(14 TeV)
    'VBFH':4078,
    'Wm':567.7,
    'Wp':888.9,
    'ZH':943.9,
    'TTH':570,
    #'ggZH':XXX,
    ###
#    'Z':6688.0*1000,
    'Z':6244.8*1000, # TO USE TURBO XSECION 6244.8
    'EWKZ':1, # TEMP: placeholder for now
    'TT2l2n':923.6*1000*0.105, # needed the 2L BR (NNLO)
    'TTln':923.6*1000*0.438,
    'TW2l2n':87.9*0.5*1000*0.105,
    'WZto2L2Q':7.568*1000, # NLO
    'ZZto2L2Nu':1.031*1000, # NLO
    'ZZto2L2Q':6.788*1000, # NLO
    'WZto3LNu':4.924*1000, # NLO
    'ZZto4L':1.39*1000, # NLO
    'WWW':0.2328*1000,
    'WWZ':0.1851*1000,
    'WZZ':0.06206*1000,
    'ZZZ':0.01591*1000,
    'TTW':0.2505*1000,
    'TTLNu-EWK':0.01769*1000,
    'TTLL_MLL-4to50':0.03949*1000,
    'TTLL_MLL-50':0.08646*1000,
    'TTZ-ZtoQQ':0.6603*1000,
}


# from Marini 31 july
#kfactor = 2094 *3 / (5378 + 1017 + 385.5 ) = 0.9264803481 (MATRIX MLL>50)

kfactor={
    'Z':0.9265,
}

def findDIR(directory,useXROOTD=False):

    print('HELLO')
    print(directory)
    counter = 0
    rootFiles = ROOT.vector('string')()
    maxFiles = 1000000000

    if(useXROOTD == True and "/data/submit/cms" in directory):
        xrd = "root://submit50.mit.edu/"
        xrdpath = directory.replace("/data/submit/cms","")
        f = check_output(['xrdfs', f'{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
        stringFiles = f.split()
        for e in range(len(stringFiles)):
            filePath = xrd + stringFiles[e]
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            if ".txt" in filePath: continue
            counter+=1
            if(counter > maxFiles): break
            rootFiles.push_back(filePath)
    else:

        filename = ("{}".format(directory)+"/*.root")
        for filenames in glob.glob(filename):
            counter+=1
            if(counter > maxFiles): break
            rootFiles.push_back(filenames)

    print(len(rootFiles))
    return rootFiles

def BuildDict(year):

    dirName="/ceph/submit/data/group/cms/store/Hmumu/v12/"
    if (str(year) == '2024' or str(year) == '2025'): dirName="/ceph/submit/data/group/cms/store/Hmumu/v15/"

    campaign_map = {
        "12022": "/NANOAODSIM/130X_mcRun3_2022_realistic_v*/*",
        "22022": "/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v*/*",
        "12023": "/NANOAODSIM/130X_mcRun3_2023_realistic_v*/*",
        "22023": "/NANOAODSIM/130X_mcRun3_2023_realistic_postBPix_v*-v*/*",
        "2024":  "/NANOAODSIM/150X_mcRun3_2024_realistic_v*-v*/*",
    }
    campaign = campaign_map.get(year, "")

    # Helper for building paths
    def path(suffix):
        return f"{dirName}{year}/{suffix}{campaign}"

    thisdict = {
        10: (findDIR(path("/VBF*Hto2Mu_*M-125_TuneCP5*_13p6TeV_powheg-pythia8")),xsecRun3['VBFH']*0.00022),
        11: (findDIR(path("/GluGlu*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggH']*0.00022),
        12: (findDIR(path("/WminusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wm']*0.00022),
        13: (findDIR(path("/WplusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wp']*0.00022),
        14: (findDIR(path("/ZH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['ZH']*0.00022),
        15: (findDIR(path("/TTH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.00022),
        ##
        20: (findDIR(path("/VBF*HtoZG*to2L*M-125_TuneCP5*_13p6TeV_powheg-pythia8")),xsecRun3['VBFH']*0.0015),
        21: (findDIR(path("/GluGlu*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggH']*0.0015),
        22: (findDIR(path("/WminusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wm']*0.0015),
        23: (findDIR(path("/WplusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wp']*0.0015),
        24: (findDIR(path("/ZH*HtoZGto2LG_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['ZH']*0.0015),
        25: (findDIR(path("/ttHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.0015),
        26: (findDIR(path("/TTH-HtoZGto2LG_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.0015),
        ## DY-QCD
        103: (findDIR(path("/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']*(1./3)),
        104: (findDIR(path("/DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']*(1./3)),
        100: (findDIR(path("/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']),
        108: (findDIR(path("/DYto2Mu_*MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8")),2219*1000), # NLO from xsecDB
        110: (findDIR(path("/DYto2Mu_*MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8")),21.65*1000), # NLO from xsecDB
        109: (findDIR(path("/DYto2Mu-2Jets*MLL-105*o160_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),49.21*1000), # NLO from xsecDB
        111: (findDIR(path("/DYto2L-2Jets_MLL-50_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),5378*1000), # LO from xsecDB
        112: (findDIR(path("/DYto2L-2Jets_MLL-50_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),1017*1000), # LO from xsecDB
        113: (findDIR(path("/DYto2L-2Jets_MLL-50_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),375*1000), # LO from xsecDB
        114: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),58*1000), # taken from Hinv these are temporary
        115: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),6.6*1000), # taken from Hinv these are temporary
        116: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.38*1000), # taken from Hinv these are temporary
        117: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.067*1000), # taken from Hinv these are temporary
        ## DY-EWK
        101: (findDIR(path("/EWK*2L2J_*TuneCH3_13p6TeV_madgraph-herwig7")),xsecRun3['EWKZ']), # when ready add the 105-120 Hrewig7 powheg dipole
        ## TTBAR
        102: (findDIR(path("/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']), # this is the default
        107: (findDIR(path("/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTln']),
        105: (findDIR(path("/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        106: (findDIR(path("/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        118: (findDIR(path("/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['TW2l2n']),
        ## DI/TRI- BOSON
        201: (findDIR(path("/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto2L2Q']),
        202: (findDIR(path("/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Nu']),
        203: (findDIR(path("/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Q']),
        204: (findDIR(path("/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto3LNu']),
        205: (findDIR(path("/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto4L']),
        211: (findDIR(path("/WWW*_TuneCP5_13p6TeV_amcatnlo*-pythia8")),xsecRun3['WWW']),
        212: (findDIR(path("/WWZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WWZ']),
        213: (findDIR(path("/WZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WZZ']),
        214: (findDIR(path("/ZZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['ZZZ']),
        221: (findDIR(path("/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['TTW']),
        222: (findDIR(path("/TTLNu-EWK_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLNu-EWK']),
        223: (findDIR(path("/TTLL_*MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-4to50']),
        224: (findDIR(path("/TTLL_*MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-50']),
        225: (findDIR(path("/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX*-pythia8")),xsecRun3['TTZ-ZtoQQ']),
    }
    # TODO: add the signal aMCNLO sample, for the training we can use the powheg


   # Data-taking runs, grouped by year
    data_runs = {
        "12022": {
            -11: "Run2022C/SingleMuon/*/*/",
            -12: "Run2022C/DoubleMuon/*/*/",
            -13: "Run2022C/Muon/*/*/",
            -14: "Run2022D/Muon/*/*/",
            -113: "Run2022C/ParkingDoubleMuonLowMass*/*/*/",
            -114: "Run2022D/ParkingDoubleMuonLowMass*/*/*/",
        },
        "22022": {
            -15: "Run2022E/Muon/*/*/",
            -16: "Run2022F/Muon/*/*/",
            -17: "Run2022G/Muon/*/*/",
            -115: "Run2022E/ParkingDoubleMuonLowMass*/*/*/",
            -116: "Run2022F/ParkingDoubleMuonLowMass*/*/*/",
            -117: "Run2022G/ParkingDoubleMuonLowMass*/*/*/",
        },
        "12023": {
            -23: "Run2023C/Muon0/*/*/",
            -24: "Run2023C/Muon1/*/*/",
            -123: "Run2023C/ParkingDoubleMuonLowMass*/*/*/",
        },
        "22023": {
            -31: "Run2023D/Muon0/*/*/",
            -32: "Run2023D/Muon1/*/*/",
            -131: "Run2023D/ParkingDoubleMuonLowMass*/*/*/",
        },
        "2024": {
            -41: "Run2024C/Muon0/*/*/",
            -42: "Run2024C/Muon1/*/*/",
            -43: "Run2024D/Muon0/*/*/",
            -44: "Run2024D/Muon1/*/*/",
            -45: "Run2024E/Muon0/*/*/",
            -46: "Run2024E/Muon1/*/*/",
            -47: "Run2024F/Muon0/*/*/",
            -48: "Run2024F/Muon1/*/*/",
            -49: "Run2024G/Muon0/*/*/",
            -50: "Run2024G/Muon1/*/*/",
            -51: "Run2024H/Muon0/*/*/",
            -52: "Run2024H/Muon1/*/*/",
            -53: "Run2024I/Muon0/*/*/",
            -54: "Run2024I/Muon1/*/*/",
            -141: "Run2024C/ParkingDoubleMuonLowMass*/*/*/",
            -142: "Run2024D/ParkingDoubleMuonLowMass*/*/*/",
            -143: "Run2024E/ParkingDoubleMuonLowMass*/*/*/",
            -144: "Run2024F/ParkingDoubleMuonLowMass*/*/*/",
            -145: "Run2024G/ParkingDoubleMuonLowMass*/*/*/",
            -146: "Run2024H/ParkingDoubleMuonLowMass*/*/*/",
            -147: "Run2024I/ParkingDoubleMuonLowMass*/*/*/",
        },
        "2025": {
            -61: "Run2025B/Muon0/*/*/*/*/*",
            -62: "Run2025B/Muon1/*/*/*/*/*",
            -63: "Run2025C/Muon0/*/*/*/*/*",
            -64: "Run2025C/Muon1/*/*/*/*/*",
            -65: "Run2025D/Muon0/*/*/*/*/*",
            -66: "Run2025D/Muon1/*/*/*/*/*",
            -67: "Run2025E/Muon0/*/*/*/*/*",
            -68: "Run2025E/Muon1/*/*/*/*/*",
            -69: "Run2025F/Muon0/*/*/*/*/*",
            -70: "Run2025F/Muon1/*/*/*/*/*",
            -71: "Run2025G/Muon0/*/*/*/*/*",
            -72: "Run2025G/Muon1/*/*/*/*/*",
            -161: "Run2025B/ParkingDoubleMuonLowMass*/*/*/*/*/*", # this is Prompt
            -162: "Run2025C/ParkingDoubleMuonLowMass*/*/*/*/*/*",
            -163: "Run2025D/ParkingDoubleMuonLowMass*/*/*/*/*/*",
            -164: "Run2025E/ParkingDoubleMuonLowMass*/*/*/*/*/*",
            -165: "Run2025F/ParkingDoubleMuonLowMass*/*/*/*/*/*",
            -166: "Run2025G/ParkingDoubleMuonLowMass*/*/*/*/*/*",
         }
    }

    # Add year-specific runs if available
    if year in data_runs:
        thisdict.update({k: (findDIR(f"{dirName}{year}/{v}"), -1) for k, v in data_runs[year].items()})
    return thisdict

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


def loadCorrectionSet(year):
    print('loadCorrectionSet()')
    import correctionlib
    correctionlib.register_pyroot_binding()
    ROOT.gInterpreter.Declare('#include "./config/sfCorrLib.h"')

    ROOT.gInterpreter.ProcessLine('auto corr_sf = MyCorrections(%d);' % (year))

    print('loadMuonScale()')

    subDirName = ""
    if year == 12022: subDirName = "Run3-22CDSep23-Summer22-NanoAODv12"
    if year == 22022: subDirName = "Run3-22EFGSep23-Summer22EE-NanoAODv12"
    if year == 12023: subDirName = "Run3-23CSep23-Summer23-NanoAODv12"
    if year == 22023: subDirName = "Run3-23DSep23-Summer23BPix-NanoAODv12"
    if year == 2024: subDirName = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15"

    if(year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024):
        print('loadMuonScale()')
        ROOT.gROOT.ProcessLine(
            f'auto cset = correction::CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/'+subDirName+'/latest/muon_scalesmearing.json.gz");'
        )
        ROOT.gInterpreter.ProcessLine('#include "./config/MuonScaRe.cc"')
        ROOT.gInterpreter.Declare('#include "./config/functionsMuCorr.h"')

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
#    else: dirJson = "./config"

    dirJson = "./config"

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
        "2025":  "Cert_Collisions2025_391658_398860_Golden.json",
    }

    fname = json_map.get(str(year))
    if fname:
        loadJSON(f"{dirJson}/cert/{fname}")
    else:
        print(f"No JSON mapping found for year={year}")
