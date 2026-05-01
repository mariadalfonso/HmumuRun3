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
    '2024':109.95, # C-I
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

btagParTL={
    '12022':0.0849,
    '22022':0.0897,
    '12023':0.0681,
    '22023':0.0683,
    '2024':0.0246, # found only UparT
    '2025':0.0246, # copy the 2024 value
}

btagParTM={
    '12022':0.4319,
    '22022':0.451,
    '12023':0.3487,
    '22023':0.3494,
    '2024':0.1272, # found only UparT
    '2025':0.1272, # copy the 2024 value
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
#    'ggZH':XXX,
#    'bbH':XXX,
    ###
#    'Z':6688.0*1000,
    'Z':6244.8*1000, # TO USE TURBO XSECION 6244.8
    'EWKZ':1.71*1000, # TEMP: placeholder for now from Andrea
    'W':67710.0*1000,
    'TT2l2n':923.6*1000*0.105, # needed the 2L BR (NNLO)
    'TTln':923.6*1000*0.438,
    'TW2l2n':87.9*0.5*1000*0.105,
    'WWto2L2Nu':11.79*1000, # NLO
    'WZto2L2Q':7.568*1000, # NLO
    'ZZto2L2Nu':1.031*1000, # NLO # Didar 1.175
    'ZZto2L2Q':6.788*1000, # NLO # Didar 0.449
    'WZto3LNu':4.924*1000, # NLO # Didara 5.297
    'ZZto4L':1.39*1000, # NLO  # Didar 1.65
    'WWW':0.2328*1000,
    'WWZ':0.1851*1000,
    'WZZ':0.06206*1000,
    'ZZZ':0.01591*1000,
    ###
    'TTW':0.2505*1000,
    'TTLNu-EWK':0.01769*1000,
    'TTLL_MLL-4to50':0.03949*1000,
    'TTLL_MLL-50':0.08646*1000,
    'TTZ-ZtoQQ':0.6603*1000,
    'TZQ':0.07968*1000,
    'TWZ3l':0.00334*1000,
    'TWZ4l':0.00167*1000,
    'TTTT':0.009652*1000,
    'TTWZ':0.002715*1000,
    'TTWW':0.008191*1000,
    'TTZZ':0.001579*1000,
    'TTW-WtoQQ':0.4678*1000,
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
    dirNameScratch="/scratch/submit/cms/mariadlf/Hmumu/v12/"
    if (str(year) == '2024' or str(year) == '2025'):
        dirName="/ceph/submit/data/group/cms/store/Hmumu/v15/"
        dirNameScratch="/scratch/submit/cms/mariadlf/Hmumu/v15/"

    campaign_map = {
        "12022": "/NANOAODSIM/130X_mcRun3_2022_realistic_v*/*",
        "22022": "/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v*/*",
        "12023": "/NANOAODSIM/130X_mcRun3_2023_realistic_v*/*",
        "22023": "/NANOAODSIM/130X_mcRun3_2023_realistic_postBPix_v*-v*/*",
        "2024":  "/NANOAODSIM/150X_mcRun3_2024_realistic_v*-v*/*",
        "2025":  "/NANOAODSIM/150X_mcRun3_2025*-v*/*",
    }
    campaign = campaign_map.get(year, "")

    # Helper for building paths
    def path(suffix):
        return f"{dirName}{year}/{suffix}{campaign}"

    def pathEx(suffix,specialCampaign):
        if year in "2024":
            return f"{dirName}{year}/{suffix}/NANOAODSIM/{specialCampaign}/*"
        else:
            return f"{dirName}{year}/{suffix}{campaign}"

    def pathScratch(suffix):
        return f"{dirNameScratch}{year}/{suffix}{campaign}"

    thisdict = {
        10: (findDIR(pathScratch("/VBF*Hto2Mu_*M-125_TuneCP5*_13p6TeV_powheg-pythia8")),xsecRun3['VBFH']*0.00022),
        11: (findDIR(pathScratch("/GluGlu*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggH']*0.00022),   # reweighting for NNLO PS should be done
        12: (findDIR(pathScratch("/WminusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wm']*0.00022),
        13: (findDIR(pathScratch("/WplusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wp']*0.00022),
        14: (findDIR(pathScratch("/ZH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['ZH']*0.00022),
        15: (findDIR(pathScratch("/TTH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.00022),
#        16: (findDIR(path("/GluGluZH-Zto2L-Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggZH']*0.00022),
#        17: (findDIR(path("/BBH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['bbH']*0.00022),
        ##
        20: (findDIR(path("/VBF*HtoZG*to2L*M-125_TuneCP5*_13p6TeV_powheg-pythia8")),xsecRun3['VBFH']*0.0015),
        21: (findDIR(path("/GluGlu*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggH']*0.0015),
        22: (findDIR(path("/WminusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wm']*0.0015),
        23: (findDIR(path("/WplusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wp']*0.0015),
        24: (findDIR(path("/ZH*HtoZGto2LG_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['ZH']*0.0015),
        25: (findDIR(path("/ttHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.0015),
        26: (findDIR(path("/TTH-HtoZGto2LG_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.0015),
        ##
        30: (findDIR(path("/VBF*Hto2Wto2L2Nu*M-125_TuneCP5*_13p6TeV_powheg-jhugen-pythia8")),xsecRun3['VBFH']*0.2137),
        31: (findDIR(path("/GluGlu*Hto2Wto2L2Nu*M-125_TuneCP5_13p6TeV_powheg-jhugen-pythia8")),xsecRun3['ggH']*0.2137),
        32: (findDIR(path("/WminusH*Wto2Q-Hto2Wto2L2Nu*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['Wm']*0.2137*(1-0.325)),
        33: (findDIR(path("/WminusH*WtoLNu-Hto2Wto2L2Nu*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['Wm']*0.2137*0.325),
        34: (findDIR(path("/WplusH*Wto2Q_Hto2Wto2L2Nu*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['Wp']*0.2137*(1-0.325)),
        35: (findDIR(path("/WplusH*WtoLNu-Hto2Wto2L2Nu*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['Wp']*0.2137*0.325),
        36: (findDIR(path("/ZH*Zto2Q-Hto2Wto2L2Nu_*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['ZH']*0.2137*0.70),
        37: (findDIR(path("/ZH*Zto2L-Hto2Wto2L2Nu_*M-125_TuneCP5_13p6TeV_powheg*-jhugen-pythia8")),xsecRun3['ZH']*0.2137*0.10),
        # this ZH add for all missing xsection
        #/GluGluZH-Zto2L-Hto2Wto2L2Nu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-jhugen-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
        #/GluGluZH-Zto2Q-Hto2Wto2L2Nu_Par-M-125_TuneCP5_13p6TeV_powhegMINLO-jhugen-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
        ## DY-QCD
        103: (findDIR(pathScratch("/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']*(1./3)), # I think this has a genCut at 130
#        103: (findDIR(path("/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']*(1./3)),
        104: (findDIR(path("/DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']*(1./3)),
        100: (findDIR(pathScratch("/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['Z']),
        108: (findDIR(path("/DYto2Mu_*MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8")),2219*1000), # NLO from xsecDB
        110: (findDIR(path("/DYto2Mu_*MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8")),21.65*1000), # NLO from xsecDB
        109: (findDIR(path("/DYto2Mu-2Jets*MLL-105*o160_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),49.21*1000*0.95), # NLO from xsecDB *0.05  reduction empirical
        111: (findDIR(path("/DYto2L-2Jets_MLL-50_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),5378*1000), # LO from xsecDB
        112: (findDIR(path("/DYto2L-2Jets_MLL-50_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),1017*1000), # LO from xsecDB
        113: (findDIR(path("/DYto2L-2Jets_MLL-50_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),385.5*1000), # LO from xsecDB
        119: (findDIR(path("/DYto2Mu-2Jets_Bin-0J-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),5378*1000*1./3), # LO from xsecDB
        120: (findDIR(path("/DYto2Mu-2Jets_Bin-1J-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),1017.5*1000*1./3), # LO from xsecDB
        121: (findDIR(path("/DYto2Mu-2Jets_Bin-2J-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),385*1000*1./3), # LO from xsecDB # on the xsecDB there is no
        114: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),110.9*1000*(0.90)), # NLO # empirically scaled
        115: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),11.36*1000*(0.08)), # NLO # empirically scaled
        116: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.6278*1000*(0.08)), # NLO # empirically scaled
        117: (findDIR(path("/DYto2L-2Jets_MLL-50_PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.08555*100*(0.08)), # NLO # empirically scaled
        122: (findDIR(path("/DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),106.6*1000*0.95), # NLO
        123: (findDIR(path("/DYto2L-2Jets_Bin-MLL-50-PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),10.83*1000*0.10), # NLO
        124: (findDIR(path("/DYto2L-2Jets_Bin-MLL-50-PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.5943*1000*0.10), # NLO
        125: (findDIR(path("/DYto2L-2Jets_Bin-MLL-50-PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),0.08108*1000*(0.10)), # NLO
        126: (findDIR(path("/DYto2Mu_*MLL-50to130_TuneCP5_13p6TeV_powhegMINNLO-pythia8-photos")),xsecRun3['Z']*(1./3)*0.95), # temporary reduction empirical
        127: (findDIR(path("/DYto2Mu_*MLL-130to200_TuneCP5_13p6TeV_powhegMINNLO-pythia8-photos")),0.60*21.65*1000), # temporaty presa quella di 120to200 empirically scaled
        128: (findDIR(path("/DYto2Mu_*MLL-200to400_TuneCP5_13p6TeV_powhegMINNLO-pythia8-photos")),3.058*1000), # temporaty presa quella di powhegV2
        129: (findDIR(path("/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['W']),
#        130: (findDIR(path("/WtoLNu-2Jets-EWK_TuneCP5_13p6TeV_madgraph-pythia8")) # no xsection so far
        ## DY-EWK
        101: (findDIR(path("/EWK*2L2J_*TuneCH3_13p6TeV_madgraph-herwig7")),xsecRun3['EWKZ']),
        99:  (findDIR(path("/EWK-2Mu2J*M2Mu-105to160*M2J-120_TuneCP5_13p6TeV_madgraph-pythia8")),0.06443*1000), # pythia8 i.e. dipole ?
        98:  (findDIR(path("/EWK*2Mu2J_*MLL-105to160_TuneCH3_13p6TeV_madgraph-herwig7")),xsecRun3['EWKZ']/120.), # added 120. empirical scaling from the DYQCD
        ## TTBAR
        102: (findDIR(pathScratch("/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']), #<-- can add the AMC@NLO
#        102: (findDIR(path("/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']), #<-- can add the AMC@NLO
        107: (findDIR(path("/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTln']),
        105: (findDIR(path("/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        106: (findDIR(path("/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        118: (findDIR(path("/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['TT2l2n']),
        ## DI/TRI- BOSON
        201: (findDIR(path("/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto2L2Q']),
        202: (findDIR(path("/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Nu']),
        203: (findDIR(path("/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Q']),
        204: (findDIR(path("/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto3LNu']),
        205: (findDIR(path("/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto4L']),
        206: (findDIR(path("/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WWto2L2Nu']),
        211: (findDIR(path("/WWW*_TuneCP5_13p6TeV_amcatnlo*-pythia8")),xsecRun3['WWW']),
        212: (findDIR(path("/WWZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WWZ']),
        213: (findDIR(path("/WZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WZZ']),
        214: (findDIR(path("/ZZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['ZZZ']),
        ## tX
        221: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8","mg35x_150X_mcRun3_2024_realistic_v2-v2")),xsecRun3['TTW']),
        222: (findDIR(path("/TTLNu-EWK_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLNu-EWK']),
        223: (findDIR(path("/TTLL_*MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-4to50']),
        224: (findDIR(path("/TTLL_*MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-50']),
        225: (findDIR(path("/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX*-pythia8")),xsecRun3['TTZ-ZtoQQ']),
        226: (findDIR(pathEx("/TZQB-Zto2L-4FS_*MLL-30_TuneCP5_13p6TeV_amcatnlo-pythia8","Madgraph_2_6_5_150X_mcRun3_2024_realistic_v2-v2")),xsecRun3['TZQ']),
        227: (findDIR(path("TWZ*Tto2Q_WtoLNu_Zto2L_DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ3l']),
        228: (findDIR(path("TWZ*TtoLNu_Wto2Q_Zto2L_DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ3l']),
        229: (findDIR(path("TWZ*TtoLNu_WtoLNu_Zto2L_DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ4l']),
        230: (findDIR(path("TTTT_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTTT']),
        231: (findDIR(path("TTWZ_TuneCP5_13p6TeV_madgraph-pythia8")),xsecRun3['TTWZ']),
        232: (findDIR(path("TTWW_TuneCP5_13p6TeV_madgraph*-pythia8")),xsecRun3['TTWW']),
        233: (findDIR(path("TTZZ_TuneCP5_13p6TeV_madgraph*-pythia8")),xsecRun3['TTZZ']),
        234: (findDIR(path("TTW-WtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),xsecRun3['TTW-WtoQQ']),
    }

    #submitted
    #TWZ-Tto2Q-WtoLNu-Zto2L-DR1_TuneCP5_13p6TeV_amcatnlo-pythia8
    #TWZ-TtoLNu-Wto2Q-Zto2L-DR1_TuneCP5_13p6TeV_amcatnlo-pythia8
    #TWZ-TtoLNu-WtoLNu-Zto2L-DR1_TuneCP5_13p6TeV_amcatnlo-pythia8

    # TODO: add the signal aMCNLO sample, for the training we can use the powheg

    folder = 'MINIv6NANOv15-v1'

   # Data-taking runs, grouped by year
    data_runs = {
        "12022": {
            -11: ("Run2022C/SingleMuon/*/*/",dirNameScratch),
            -12: ("Run2022C/DoubleMuon/*/*/",dirNameScratch),
            -13: ("Run2022C/Muon/*/*/",dirNameScratch),
            -14: ("Run2022D/Muon/*/*/",dirNameScratch),
            -113: ("Run2022C/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -114: ("Run2022D/ParkingDoubleMuonLowMass*/*/*/",dirName),
        },
        "22022": {
            -15: ("Run2022E/Muon/*/*/",dirNameScratch),
            -16: ("Run2022F/Muon/*/*/",dirNameScratch),
            -17: ("Run2022G/Muon/*/*/",dirNameScratch),
            -115: ("Run2022E/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -116: ("Run2022F/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -117: ("Run2022G/ParkingDoubleMuonLowMass*/*/*/",dirName),
        },
        "12023": {
            -23: ("Run2023C/Muon0/*/*/",dirNameScratch),
            -24: ("Run2023C/Muon1/*/*/",dirNameScratch),
            -123: ("Run2023C/ParkingDoubleMuonLowMass*/*/*/",dirName),
        },
        "22023": {
            -31: ("Run2023D/Muon0/*/*/",dirNameScratch),
            -32: ("Run2023D/Muon1/*/*/",dirNameScratch),
            -131: ("Run2023D/ParkingDoubleMuonLowMass*/*/*/",dirName),
        },
        "2024": {
            -41: ("Run2024C/Muon0/*/*/",dirNameScratch),
            -42: ("Run2024C/Muon1/*/*/",dirNameScratch),
            -43: ("Run2024D/Muon0/*/*/",dirNameScratch),
            -44: ("Run2024D/Muon1/*/*/",dirNameScratch),
            -45: ("Run2024E/Muon0/*/*/",dirNameScratch),
            -46: ("Run2024E/Muon1/*/*/",dirNameScratch),
            -47: ("Run2024F/Muon0/*/*/",dirNameScratch),
            -48: ("Run2024F/Muon1/*/*/",dirNameScratch),
            -49: ("Run2024G/Muon0/*/*/",dirNameScratch),
            -50: ("Run2024G/Muon1/*/*/",dirNameScratch),
            -51: ("Run2024H/Muon0/*/*/",dirNameScratch),
            -52: ("Run2024H/Muon1/*/*/",dirNameScratch),
            -53: ("Run2024I/Muon0/*/*/",dirNameScratch),
            -54: ("Run2024I/Muon1/*/*/",dirNameScratch),
            -141: ("Run2024C/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -142: ("Run2024D/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -143: ("Run2024E/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -144: ("Run2024F/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -145: ("Run2024G/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -146: ("Run2024H/ParkingDoubleMuonLowMass*/*/*/",dirName),
            -147: ("Run2024I/ParkingDoubleMuonLowMass*/*/*/",dirName),
        },
        "2025": {
            -61: ("Run2025B/Muon0/*/*/*/*/*",dirNameScratch),
            -62: ("Run2025B/Muon1/*/*/*/*/*",dirNameScratch),
            -63: ("Run2025C/Muon0/*/*/*/*/*",dirNameScratch),
            -64: ("Run2025C/Muon1/*/*/*/*/*",dirNameScratch),
            -65: ("Run2025D/Muon0/*/*/*/*/*",dirNameScratch),
            -66: ("Run2025D/Muon1/*/*/*/*/*",dirNameScratch),
            -67: ("Run2025E/Muon0/*/*/*/*/*",dirNameScratch),
            -68: ("Run2025E/Muon1/*/*/*/*/*",dirNameScratch),
            #
            -69: ("Run2025F/Muon0/*/*/*/*/*",dirNameScratch),
            -70: ("Run2025F/Muon1/*/*/*/*/*",dirNameScratch),
            -71: ("Run2025G/Muon0/*/*/*/*/*",dirNameScratch),
            -72: ("Run2025G/Muon1/*/*/*/*/*",dirNameScratch),
            -161: ("Run2025B/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -162: ("Run2025C/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -163: ("Run2025D/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -164: ("Run2025E/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -165: ("Run2025F/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -166: ("Run2025G/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
        }
    }

    # Add year-specific runs if available
    if year in data_runs:
        for run_code, (pattern, baseDir) in data_runs[year].items():
            # Build full search path
            full_path = f"{baseDir}{year}/{pattern}"

            # Search directory
            found = findDIR(full_path)

            if not found:
                print(f"⚠️ WARNING: No files found for run {run_code}: {full_path}")

            thisdict[run_code] = (found, -1)

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

    if(year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024 or year==2025):
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
        f'auto csetH = correction::CorrectionSet::from_file("./config/NNLOPS_reweight_13p6TeV.json");'
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
    }

    fname = json_map.get(str(year))
    if fname:
        loadJSON(f"{dirJson}/cert/{fname}")
    else:
        print(f"No JSON mapping found for year={year}")
