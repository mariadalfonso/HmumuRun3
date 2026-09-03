"""
datasets.py — dataset selection and definitions for the HmumuRun3 analysis.

Single module for everything "which samples and where":

  Selection  (which IDs to run for a year/mode):
    - getMCList(year, mode)   -> list of MC sample IDs
    - getDataList(year, mode) -> list of data (negative) sample IDs

  Definition (what each ID is: files + cross-section):
    - xsecRun3   : cross-section dictionary
    - findDIR    : glob/xrootd helper resolving a directory to a ROOT file list
    - BuildDict(year) : maps each sample ID -> (file list, xsec)

Execution (the event loop) lives in Hmm.loopOnDataset.
"""

import ROOT
import glob
from subprocess import check_output
import sys


# ===========================================================================
#  Cross-sections
# ===========================================================================
xsecRun3={
    'ggH':52230, # 0.4 × σ(13 TeV) + 0.6 × σ(14 TeV)
    'VBFH':4078,
    'Wm':567.7,
    'Wp':888.9,
    'ZH':943.9,
    'TTH':570,
    'ggZH':217.09, # from xsecDB
    'bbH':568, # 0.568 /pb
    ###
#    'Z':6688.0*1000,
    'Z':6244.8*1000, # TO USE TURBO XSECION 6244.8
    'EWKZ':1.71*1000, # TEMP: placeholder for now from Andrea
    'W':67710.0*1000,
    'EWKW':91.48*1000,
    'TT2l2n':923.6*1000*0.105, # needed the 2L BR (NNLO)
    'TTln':923.6*1000*0.438,
    'TW2l2n':87.9*0.5*1000*0.105,
    #
    'WWto2L2Nu':11.79*1000, # NLO
    'WZto2L2Q':7.568*1000, # NLO
    'ZZto2L2Nu':1.031*1000, # NLO # Didar 1.175
    'ZZto2L2Q':6.788*1000, # NLO # Didar 0.449
    'WZto3LNu':4.924*1000, # NLO # Didara 5.297
    'ZZto4L':1.39*1000, # NLO  # Didar 1.65
    #
    'WZto2L2Q-1Jets':6.724*1000,
    'ZZto2L2Q-1Jets':3.86*1000,
    'VVto2L2Nu_MLL-50-1Jets':10.48*1000,
    'VVto2L2Nu_MLL-4to50-1Jets':3.304*1000,
    'WZto3LNu-1Jets':5.315*1000,
    'ZZto4L-1Jets':1.499*1000,
    #
    'WWW':0.2328*1000,
    'WWZ':0.1851*1000,
    'WZZ':0.06206*1000,
    'ZZZ':0.01591*1000,
    ###
    'TTW':0.2505*1000,
    'TTLNu-EWK':0.01769*1000,
    'TTW-WtoQQ':0.4678*1000,
    'TTLL_MLL-4to50':0.03949*1000,
    'TTLL_MLL-50':0.08646*1000,
    'TTZ-ZtoQQ':0.6603*1000,
    'TTNuNu':0.1638*1000,
    'TZQ':0.07968*1000,
    'TWZ3l':0.00334*1000,
    'TWZ4l':0.00167*1000,
    'TTTT':0.009652*1000,
    'TTWZ':0.002715*1000,
    'TTWW':0.008191*1000,
    'TTZZ':0.001579*1000,
    'TTTWm':0.0006029*1000,
    'TTTWp':0.0006061*1000,
}


# ===========================================================================
#  File discovery
# ===========================================================================
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


# ===========================================================================
#  Dataset definition: ID -> (files, xsec)
# ===========================================================================
def BuildDict(year):

    dirName="/ceph/submit/data/group/cms/store/Hmumu/v12/"
    dirNameScratch="/scratch/submit/cms/mariadlf/Hmumu/v12/"
    if (str(year) == '2024' or str(year) == '2025' or str(year) == '2026'):
        dirName="/ceph/submit/data/group/cms/store/Hmumu/v15/"
        dirNameScratch="/scratch/submit/cms/mariadlf/Hmumu/v15/"

    campaign_map = {
        "12022": "130X_mcRun3_2022_realistic_v*/*",
        "22022": "130X_mcRun3_2022_realistic_postEE_v*/*",
        "12023": "130X_mcRun3_2023_realistic_v*/*",
        "22023": "130X_mcRun3_2023_realistic_postBPix_v*-v*/*",
        "2024":  "150X_mcRun3_2024_realistic_v*-v*/*",
        "2025":  "150X_mcRun3_2025*-v*/*", # noMC
        "2026":  "150X_mcRun3_2025*-v*/*", # noMC
    }
    campaign = campaign_map.get(year, "")

    # Helper for building paths
    def path(suffix):
        return f"{dirName}{year}/{suffix}/NANOAODSIM/{campaign}"

    def pathEx(suffix,specialCampaignPrefix, year_):
        if year in year_:
            return f"{dirName}{year}/{suffix}/NANOAODSIM/{specialCampaignPrefix}{campaign}"
        else:
            return f"{dirName}{year}/{suffix}/NANOAODSIM/{campaign}"

    def pathScratch(suffix):
        return f"{dirNameScratch}{year}/{suffix}/NANOAODSIM/{campaign}"

    thisdict = {
        10: (findDIR(pathScratch("/VBF*Hto2Mu_*M-125_TuneCP5*_13p6TeV_powheg-pythia8")),xsecRun3['VBFH']*0.00022),
        11: (findDIR(pathScratch("/GluGlu*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggH']*0.00022),   # reweighting for NNLO PS should be done
        12: (findDIR(pathScratch("/WminusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wm']*0.00022),
        13: (findDIR(pathScratch("/WplusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['Wp']*0.00022),
        14: (findDIR(pathScratch("/ZH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8")),xsecRun3['ZH']*0.00022),
        15: (findDIR(pathScratch("/TTH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTH']*0.00022),
#        16: (findDIR(path("/GluGluZH-Zto2L-Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ggZH']*0.00022),
        17: (findDIR(path("/BBH-Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['bbH']*0.00022),
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
        # 102 is now free
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
        130: (findDIR(path("/WtoLNu-2Jets-EWK_TuneCP5_13p6TeV_madgraph-pythia8")),xsecRun3['EWKW']),
        ## DY-EWK
        101: (findDIR(path("/EWK*2L2J_*TuneCH3_13p6TeV_madgraph-herwig7")),xsecRun3['EWKZ']),
        99:  (findDIR(path("/EWK-2Mu2J*M2Mu-105to160*M2J-120_TuneCP5_13p6TeV_madgraph-pythia8")),0.06443*1000), # pythia8 i.e. dipole ?
        98:  (findDIR(path("/EWK*2Mu2J_*MLL-105to160_TuneCH3_13p6TeV_madgraph-herwig7")),0.0527*1000), # from Filippo, otherwise if xsecRun3['EWKZ']/120. is empirical scaling from the DYQCD
        ## TTBAR
        140: (findDIR(pathScratch("/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        141: (findDIR(path("/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['TT2l2n']),
        142: (findDIR(path("/TTto2L2Nu-3Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8")),xsecRun3['TT2l2n']),
        143: (findDIR(path("/TTto2L2Nu_TuneCP5CR2_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        144: (findDIR(path("/TTto2L2Nu_TuneCP5CR1_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        145: (findDIR(path("/TTto2L2Nu_TuneCP5Down_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        146: (findDIR(path("/TTto2L2Nu_TuneCP5Up_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        147: (findDIR(path("/TTto2L2Nu_*ERD*On*_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        148: (findDIR(path("/TTto2L2Nu_*Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        149: (findDIR(path("/TTto2L2Nu_*Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TT2l2n']),
        107: (findDIR(path("/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TTln']),  # this also has the same syst as the
        105: (findDIR(path("/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        106: (findDIR(path("/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['TW2l2n']),
        ## DI/TRI- BOSON
        201: (findDIR(path("/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto2L2Q']),
        202: (findDIR(path("/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Nu']),
        203: (findDIR(path("/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto2L2Q']),
        204: (findDIR(path("/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WZto3LNu']),
        205: (findDIR(path("/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['ZZto4L']),
        206: (findDIR(path("/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8")),xsecRun3['WWto2L2Nu']),
        #
        207: (findDIR(path("/WZto2L2Q-1Jets-4FS_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['WZto2L2Q-1Jets']),
        208: (findDIR(path("/ZZto2L2Q-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['ZZto2L2Q-1Jets']),
        209: (findDIR(path("/VVto2L2Nu-1Jets-4FS_*MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['VVto2L2Nu_MLL-50-1Jets']),
        210: (findDIR(path("/VVto2L2Nu-1Jets-4FS_*MLL-4to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['VVto2L2Nu_MLL-4to50-1Jets']),
        211: (findDIR(path("/WZto3LNu-1Jets-4FS_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['WZto3LNu-1Jets']),
        212: (findDIR(path("/ZZto4L-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8")),xsecRun3['ZZto4L-1Jets']),
        ## VVV
        213: (findDIR(path("/WZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WZZ']),
        214: (findDIR(path("/ZZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['ZZZ']),
        215: (findDIR(path("/WWW*_TuneCP5_13p6TeV_amcatnlo*-pythia8")),xsecRun3['WWW']),
        216: (findDIR(path("/WWZ*_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['WWZ']),
        ## ttV
        221: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8","mg35x_","2024")),xsecRun3['TTW']), # this we have syst
        222: (findDIR(path("/TTLNu-EWK_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLNu-EWK']),
        234: (findDIR(path("TTW-WtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),xsecRun3['TTW-WtoQQ']),
        #
        223: (findDIR(path("/TTLL_*MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-4to50']),
        224: (findDIR(path("/TTLL_*MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTLL_MLL-50']),
        225: (findDIR(path("/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX*-pythia8")),xsecRun3['TTZ-ZtoQQ']),                  # nevents 22EE 1 334 574
        238: (findDIR(path("TTZ-ZtoQQ-TTto2L2Nu-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),xsecRun3['TTZ-ZtoQQ']*0.105), # nevents 22EE 4 865 313
        239: (findDIR(path("TTZ-ZtoQQ-TTtoLNu2Q-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),xsecRun3['TTZ-ZtoQQ']*0.444), # nevents 22EE 5 387 360
        237: (findDIR(path("TTNuNu_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTNuNu']),
        ## tZX + TT+(TT,VV,TW)
        226: (findDIR(pathEx("/TZQB-Zto2L-4FS_*MLL-30_TuneCP5_13p6TeV_amcatnlo-pythia8","Madgraph_2_6_5_","2024")),xsecRun3['TZQ']),
        227: (findDIR(path("TWZ*Tto2Q*WtoLNu*Zto2L*DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ3l']),
        228: (findDIR(path("TWZ*TtoLNu*Wto2Q*Zto2L*DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ3l']),
        229: (findDIR(path("TWZ*TtoLNu*WtoLNu*Zto2L*DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TWZ4l']),
        #
        230: (findDIR(path("TTTT_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTTT']),
        231: (findDIR(path("TTWZ_TuneCP5_13p6TeV_madgraph-pythia8")),xsecRun3['TTWZ']),
        232: (findDIR(path("TTWW_TuneCP5_13p6TeV_madgraph*-pythia8")),xsecRun3['TTWW']),
        233: (findDIR(path("TTZZ_TuneCP5_13p6TeV_madgraph*-pythia8")),xsecRun3['TTZZ']),
        235: (findDIR(path("TTTWminus-DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTTWm']),
        236: (findDIR(path("TTTWplus-DR1_TuneCP5_13p6TeV_amcatnlo-pythia8")),xsecRun3['TTTWp']),
        #
        242: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5CR1_13p6TeV_amcatnloFXFX-pythia8","mg35x_",year)),xsecRun3['TTW']), # TTW we have syst
        243: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5CR2_13p6TeV_amcatnloFXFX-pythia8","mg35x_",year)),xsecRun3['TTW']), # TTW we have syst
        244: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5Up_13p6TeV_amcatnloFXFX-pythia8","mg35x_",year)),xsecRun3['TTW']), # TTW we have syst
        245: (findDIR(pathEx("/TTLNu-1Jets_TuneCP5Down_13p6TeV_amcatnloFXFX-pythia8","mg35x_",year)),xsecRun3['TTW']), # TTW we have syst
        #
        246: (findDIR(path("TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),4.22*1000),
        247: (findDIR(path("TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),0.41*1000),
        248: (findDIR(path("TTG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),0.13*1000),
        249: (findDIR(path("TTG-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8")),4.629*1000),
    }

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
            -69: ("Run2025F/Muon0/*/*/*/*/*",dirNameScratch),
            -70: ("Run2025F/Muon1/*/*/*/*/*",dirNameScratch),
            -71: ("Run2025G/Muon0/*/*/*/*/*",dirNameScratch),
            -72: ("Run2025G/Muon1/*/*/*/*/*",dirNameScratch),
            #
            -161: ("Run2025B/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -162: ("Run2025C/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -163: ("Run2025D/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -164: ("Run2025E/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -165: ("Run2025F/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
            -166: ("Run2025G/ParkingDoubleMuonLowMass*/*/*/*/*/*",dirName),
        },
        "2026": {
            -81: ("Run2026B/Muon0/*/*/*/*/*",dirNameScratch),
            -82: ("Run2026B/Muon1/*/*/*/*/*",dirNameScratch),
            -83: ("Run2026B/Muon2/*/*/*/*/*",dirNameScratch),
            -84: ("Run2026B/Muon3/*/*/*/*/*",dirNameScratch),
            -85: ("Run2026D/Muon0/*/*/*/*/*",dirNameScratch),
            -86: ("Run2026D/Muon1/*/*/*/*/*",dirNameScratch),
            -87: ("Run2026D/Muon2/*/*/*/*/*",dirNameScratch),
            -88: ("Run2026D/Muon3/*/*/*/*/*",dirNameScratch),
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


# ===========================================================================
#  Sample selection: which IDs to run for a given (year, mode)
#  (moved from samples.py — verbatim from the original loopOnDataset)
# ===========================================================================


def getMCList(year, mode):
    """MC sample IDs to process for this (year, mode). Verbatim from the
    original loopOnDataset construction."""

    mc = []
    mc.extend([10,11,12,13,14,15,17])
    if year in ["12022", "22022", "12023", "22023"]: mc.extend([20,21,22,23,24,25]) #Zgamma
    if year=="2024": mc.extend([20,21,22,23,24,26]) #Zgamma
    #if year=="2024": mc.extend([30,31,32,33,34,35,36,37])

    if mode == "isVhad":
        if year in ["2024"]: mc.extend([122,123,124,125]) # extra DY pt binned for ML training (new)
        if year in ["12022", "22022", "12023", "22023"]: mc.extend([114,115,116,117]) # extra DY pt binned for ML training (check Zeynep xsection ..)
    else:
        if year=="2024": mc.extend([103,104])   #DY
        else: mc.extend([100])    #DY madgpragh

    mc.extend([101])    #DY EWK
    if mode == "isVBF": mc.extend([99,98])    #DY EWK

    mc.extend([140])    #TT2l (was 102)

    mc.extend([201,202,203,204,205,206]) #VV powheg
    mc.extend([213,214,215,216])     #VVV
    mc.extend([221,222,223,224,226,227,228,229,230,231,232,233,234,235,236,237])    #ttV

    if year=="2024": mc.extend([225]) #ttZqq (for 22-23 use the TT lep binned)
    else: mc.extend([238,239])
    if year=="2024": mc.extend([249]) #ttGamma
    else: mc.extend([246,247,248])

    mc.extend([107,105,106]) # tt1l, tW

    '''
    # below for training
    if mode == "isTThad" or mode == "isTTlep" or mode == "isZinv" or mode == "isVhad": mc.extend([141,142]) # extra TTbar (2024 not there yet)
    if mode == "isVlep": mc.extend([207,208,209,210,211,212]) #VV amcNLO

#    if year in ["2024"]:
#        if mode == "isVhad": mc.extend([119,120,121])  # extra DY jet binned for ML training
#    if year in ["12022", "22022", "12023", "22023"]:
#        if mode == "isVhad": mc.extend([111,112,113])  # extra DY jet binned for ML training

    if mode == "isTTlep" or mode == "isVlep" or mode == "isZinv" or mode == "isTThad":
        if year in ["12022", "22022", "12023", "22023"]: mc.extend([242,243,244,245])  # extra TTW syst var
        mc.extend([143,144,145,146,147,148,149])  # extra TT2l syst var

    if mode == "isVBF" or mode == "isGGH" or mode == "isVhad" or mode == "isTThad": mc.extend([109]) # extra DY jet mass binned
#    if mode == "isVBF" or mode == "isGGH": mc.extend([108,110]) # extra DY jet mass binned (no point those will break MVA)

    if mode == "isGGH" or mode == "isTThad":
        mc.extend([126,127]) # MINNLO
        if year in ["12022", "2024"]: mc.extend([128])
    '''

    return mc


def getDataList(year, mode):
    """Data (negative) sample IDs for this year.

    NOTE: in the original loopOnDataset the data loop was disabled
    (`data = []`, with `data_map.get(year, [])` commented out). Here the real
    map is returned; enable/disable at the call site in loopOnDataset.
    """

    data_map = {
        "12022": [-11, -12, -13, -14],
        "22022": [-15, -16, -17],
        "12023": [-23, -24],
        "22023": [-31, -32],
        "2024":  list(range(-41, -55, -1)),  # generates -41 to -54
        "2025":  list(range(-61, -73, -1)),  # generates -61 to -70
        "2026":  list(range(-81, -89, -1)),  # generates -81 to -88
    }

    data_map_jpsi = {
#        "12022": [-11, -12, -13, -14],
        "22022": [-116,-117], # missing (E -115)
        "12023": [-123],
        "22023": [-131],
        "2024":  list(range(-141, -148, -1)),
        "2025":  list(range(-161, -167, -1)),  # generates -61 to -70
    }

    return data_map.get(year, [])
