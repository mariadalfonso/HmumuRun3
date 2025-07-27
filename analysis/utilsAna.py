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

lumis={
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, #C
    '22023':9.451, #D
}

xsecRun2={
    'ggH':48580,
    'VBFH':3781.7,
}

# Run3 xsection from https://xsecdb-xsdb-official.app.cern.ch/xsdb/
#https://arxiv.org/pdf/2402.09955
xsecRun3={
    'ggH':52230, # 0.4 × σ(13 TeV) + 0.6 × σ(14 TeV)
    'VBFH':4078,
    'Wm':567.7,
    'Wp':888.9,
    'ZH':943.9,
    'TTH':570,
    #'ggZH':XXX,
    ###
    'Z':6688.0*1000,
    'EWKZ':1, # TEMP: placeholder for now
    'TT2l2n':762*1000,
    'WZto2L2Q':7.568*1000,
    'ZZto2L2Nu':1.031*1000,
    'ZZto2L2Q':6.788*1000,
    'WZto3LNu':4.924*1000
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
    campaign=""
    if (str(year) == '12022'): campaign="/NANOAODSIM/130X_mcRun3_2022_realistic_v*/*"
    if (str(year) == '22022'): campaign="/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v*/*"
    if (str(year) == '12023'): campaign="/NANOAODSIM/130X_mcRun3_2023_realistic_v*/*"
    if (str(year) == '22023'): campaign="/NANOAODSIM/130X_mcRun3_2023_realistic_postBPix_v*-v*/*"
    
    thisdict = {
        10: (findDIR(dirName+"/"+str(year)+"/VBFHto2Mu_M-125_TuneCP5*_13p6TeV_powheg-pythia8"+campaign),xsecRun3['VBFH']*0.00022),
        11: (findDIR(dirName+"/"+str(year)+"/GluGluHto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['VBFH']*0.00022),
        12: (findDIR(dirName+"/"+str(year)+"/WminusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wm']*0.00022),
        13: (findDIR(dirName+"/"+str(year)+"/WplusH_Hto2Mu_WtoAll_M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wp']*0.00022),
        14: (findDIR(dirName+"/"+str(year)+"/ZH_Hto2Mu_ZtoAll_M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['ZH']*0.00022),
##        15: (findDIR(dirName+"/"+str(year)+"/..../"+campaign),xsecRun3['ggZH']),
        15: (findDIR(dirName+"/"+str(year)+"/TTH_Hto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TTH']*0.00022),
        #
        100: (findDIR(dirName+"/"+str(year)+"/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+campaign),xsecRun3['Z']),
        101: (findDIR(dirName+"/"+str(year)+"/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7"+campaign),xsecRun3['EWKZ']),
        102: (findDIR(dirName+"/"+str(year)+"/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TT2l2n']),
        201: (findDIR(dirName+"/"+str(year)+"/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['WZto2L2Q']),
        202: (findDIR(dirName+"/"+str(year)+"/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ZZto2L2Nu']),
        203: (findDIR(dirName+"/"+str(year)+"/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ZZto2L2Q']),
        204: (findDIR(dirName+"/"+str(year)+"/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['WZto3LNu']),
    }
    # TODO: add the signal aMCNLO sample, for the training we can use the powheg

    if(year == '12022'):
        dict_ = {
            -11: (findDIR(dirName+str(year)+"/Run2022C/SingleMuon/*/*/"),-1),
            -12: (findDIR(dirName+str(year)+"/Run2022C/DoubleMuon/*/*/"),-1),
            -13: (findDIR(dirName+str(year)+"/Run2022C/Muon/*/*/"),-1),
            -14: (findDIR(dirName+str(year)+"/Run2022D/Muon/*/*/"),-1),
        }
        thisdict.update(dict_)

    if(year == '22022'):
        dict_ = {
            -15: (findDIR(dirName+str(year)+"/Run2022E/Muon/*/*/"),-1),
            -16: (findDIR(dirName+str(year)+"/Run2022F/Muon/*/*/"),-1),
            -17: (findDIR(dirName+str(year)+"/Run2022G/Muon/*/*/"),-1),            
        }
        thisdict.update(dict_)

    if(year == '12023'):
        dict_ = {
            -21: (findDIR(dirName+str(year)+"/Run2023B/Muon0/*/*/"),-1),
            -22: (findDIR(dirName+str(year)+"/Run2023B/Muon1/*/*/"),-1),
            -23: (findDIR(dirName+str(year)+"/Run2023C/Muon0/*/*/"),-1),
            -24: (findDIR(dirName+str(year)+"/Run2023C/Muon1/*/*/"),-1),                        
        }
        thisdict.update(dict_)

    # to finalize
    if(year == '22023'):
        dict_ = {
            -31: (findDIR(dirName+str(year)+"/Run2023D/Muon0/*/*/"),-1),
            -32: (findDIR(dirName+str(year)+"/Run2023D/Muon1/*/*/"),-1),                                    
        }
        thisdict.update(dict_)

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
    # to add the JEC scale and uncertainties

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
    print('HELLO readDataQuality', year)
    dirJson = "./config"
    if(str(year) == '2018'):
        loadJSON("{}/cert/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt".format(dirJson))
    if(str(year) == '2017'):
        loadJSON("{}/cert/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt".format(dirJson))
    if(str(year) == '22016' or year == '12016'):
        loadJSON("{}/cert/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt".format(dirJson))
    if(str(year) == '12022'):
        loadJSON("{}/cert/Cert_Collisions2022_355100_357900_eraBCD_Golden.json".format(dirJson))
    if(str(year) == '22022'):
        loadJSON("{}/cert/Cert_Collisions2022_359022_362760_eraEFG_Golden.json".format(dirJson))
    if(str(year) == '12023'):
        loadJSON("{}/cert/Cert_Collisions2023_366403_369802_eraBC_Golden.json".format(dirJson))
    if(str(year) == '22023'):
        loadJSON("{}/cert/Cert_Collisions2023_369803_370790_eraD_Golden.json".format(dirJson))

