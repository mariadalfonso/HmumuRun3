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

# lumis with golden json
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#ReReco_ERAS_A_B_C_D_E
lumis={
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, #C
    '22023':9.451, #D
    '2024':107.3, #C-I
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
    'WZto3LNu':4.924*1000,
    'ZZto4L':1.39*1000,
    'WWW':0.2328*1000,
    'WWZ':0.1851*1000,
    'WZZ':0.06206*1000,
    'ZZZ':0.01591*1000,
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
    if (str(year) == '2024'): dirName="/ceph/submit/data/group/cms/store/Hmumu/v15/"

    campaign=""
    if (str(year) == '12022'): campaign="/NANOAODSIM/130X_mcRun3_2022_realistic_v*/*"
    if (str(year) == '22022'): campaign="/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v*/*"
    if (str(year) == '12023'): campaign="/NANOAODSIM/130X_mcRun3_2023_realistic_v*/*"
    if (str(year) == '22023'): campaign="/NANOAODSIM/130X_mcRun3_2023_realistic_postBPix_v*-v*/*"
    if (str(year) == '2024'): campaign="/NANOAODSIM/150X_mcRun3_2024_realistic_v*-v*/*"

    thisdict = {
        10: (findDIR(dirName+"/"+str(year)+"/VBF*Hto2Mu_*M-125_TuneCP5*_13p6TeV_powheg-pythia8"+campaign),xsecRun3['VBFH']*0.00022),
        11: (findDIR(dirName+"/"+str(year)+"/GluGlu*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ggH']*0.00022),
        12: (findDIR(dirName+"/"+str(year)+"/WminusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wm']*0.00022),
        13: (findDIR(dirName+"/"+str(year)+"/WplusH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wp']*0.00022),
        14: (findDIR(dirName+"/"+str(year)+"/ZH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['ZH']*0.00022),
        15: (findDIR(dirName+"/"+str(year)+"/TTH*Hto2Mu_*M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TTH']*0.00022),
        ##
        20: (findDIR(dirName+"/"+str(year)+"/VBF*HtoZG*to2L*M-125_TuneCP5*_13p6TeV_powheg-pythia8"+campaign),xsecRun3['VBFH']*0.0015),
        21: (findDIR(dirName+"/"+str(year)+"/GluGlu*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ggH']*0.0015),
        22: (findDIR(dirName+"/"+str(year)+"/WminusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wm']*0.0015),
        23: (findDIR(dirName+"/"+str(year)+"/WplusH*HtoZG*to2L*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['Wp']*0.0015),
        24: (findDIR(dirName+"/"+str(year)+"/ZH*HtoZGto2LG_*M-125_TuneCP5_13p6TeV_powheg*-pythia8"+campaign),xsecRun3['ZH']*0.0015),
        25: (findDIR(dirName+"/"+str(year)+"/ttHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TTH']*0.0015),
        26: (findDIR(dirName+"/"+str(year)+"/TTH-HtoZGto2LG_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TTH']*0.0015),
        ##
        103: (findDIR(dirName+"/"+str(year)+"/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+campaign),xsecRun3['Z']*(1./3)),
        104: (findDIR(dirName+"/"+str(year)+"/DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+campaign),xsecRun3['Z']*(1./3)),
        100: (findDIR(dirName+"/"+str(year)+"/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+campaign),xsecRun3['Z']),
        101: (findDIR(dirName+"/"+str(year)+"/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7"+campaign),xsecRun3['EWKZ']),
        102: (findDIR(dirName+"/"+str(year)+"/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['TT2l2n']),
        #
        201: (findDIR(dirName+"/"+str(year)+"/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['WZto2L2Q']),
        202: (findDIR(dirName+"/"+str(year)+"/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ZZto2L2Nu']),
        203: (findDIR(dirName+"/"+str(year)+"/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ZZto2L2Q']),
        204: (findDIR(dirName+"/"+str(year)+"/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['WZto3LNu']),
        205: (findDIR(dirName+"/"+str(year)+"/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8"+campaign),xsecRun3['ZZto4L']),
        211: (findDIR(dirName+"/"+str(year)+"/WWW*_TuneCP5_13p6TeV_amcatnlo*-pythia8"+campaign),xsecRun3['WWW']),
        212: (findDIR(dirName+"/"+str(year)+"/WWZ*_TuneCP5_13p6TeV_amcatnlo-pythia8"+campaign),xsecRun3['WWZ']),
        213: (findDIR(dirName+"/"+str(year)+"/WZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8"+campaign),xsecRun3['WZZ']),
        214: (findDIR(dirName+"/"+str(year)+"/ZZZ*_TuneCP5_13p6TeV_amcatnlo-pythia8"+campaign),xsecRun3['ZZZ']),
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
#            -21: (findDIR(dirName+str(year)+"/Run2023B/Muon0/*/*/"),-1),
#            -22: (findDIR(dirName+str(year)+"/Run2023B/Muon1/*/*/"),-1),
            -23: (findDIR(dirName+str(year)+"/Run2023C/Muon0/*/*/"),-1),
            -24: (findDIR(dirName+str(year)+"/Run2023C/Muon1/*/*/"),-1),                        
        }
        thisdict.update(dict_)

    if(year == '22023'):
        dict_ = {
            -31: (findDIR(dirName+str(year)+"/Run2023D/Muon0/*/*/"),-1),
            -32: (findDIR(dirName+str(year)+"/Run2023D/Muon1/*/*/"),-1),                                    
        }
        thisdict.update(dict_)

    if(year == '2024'):
        dict_ = {
            -41: (findDIR(dirName+str(year)+"/Run2024C/Muon0/*/*/"),-1),
            -42: (findDIR(dirName+str(year)+"/Run2024C/Muon1/*/*/"),-1),
            -43: (findDIR(dirName+str(year)+"/Run2024D/Muon0/*/*/"),-1),
            -44: (findDIR(dirName+str(year)+"/Run2024D/Muon1/*/*/"),-1),
            -45: (findDIR(dirName+str(year)+"/Run2024E/Muon0/*/*/"),-1),
            -46: (findDIR(dirName+str(year)+"/Run2024E/Muon1/*/*/"),-1),
            -47: (findDIR(dirName+str(year)+"/Run2024F/Muon0/*/*/"),-1),
            -48: (findDIR(dirName+str(year)+"/Run2024F/Muon1/*/*/"),-1),
            -49: (findDIR(dirName+str(year)+"/Run2024G/Muon0/*/*/"),-1),
            -50: (findDIR(dirName+str(year)+"/Run2024G/Muon1/*/*/"),-1),
            -51: (findDIR(dirName+str(year)+"/Run2024H/Muon0/*/*/"),-1),
            -52: (findDIR(dirName+str(year)+"/Run2024H/Muon1/*/*/"),-1),
            -53: (findDIR(dirName+str(year)+"/Run2024I/Muon0/*/*/"),-1),
            -54: (findDIR(dirName+str(year)+"/Run2024I/Muon1/*/*/"),-1),
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
    ROOT.gInterpreter.Declare('''
        #ifndef MYFUN
        #define MYFUN
        Vec_f computeJECcorrection(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_rawFactor, Vec_f jet_eta, Vec_f jet_phi, Vec_f jet_area, float rho, float run, bool isData, string year, string mc){
        Vec_f new_jet(jet_pt.size(), 1.0);
        Vec_f raw_jet(jet_pt.size(), 1.0);
        for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) {
             raw_jet[idx] = jet_pt[idx] * (1.0 - jet_rawFactor[idx]);
             new_jet[idx] = raw_jet[idx] * corrSFs.eval_jetCORR(jet_area[idx], jet_eta[idx], jet_phi[idx], raw_jet[idx], rho, isData, run, year, mc );
        }
        return new_jet;
        }

        Vec_f computeJECuncertainties(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_eta){
        Vec_f new_jet_delta(jet_pt.size(), 1.0);
        int type = 0;
        for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) new_jet_delta[idx] = corrSFs.eval_jesUnc(jet_eta[idx], jet_pt[idx], type );
        return new_jet_delta;
        }

        Vec_b cleaningJetVetoMapMask(const Vec_f& jet_eta, const Vec_f& jet_phi, const string year) {
        Vec_b jet_vetoMap_mask(jet_eta.size(), true);
        for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {
        double jetVetoMap = corr_sf.eval_jetVeto(jet_eta[idx], jet_phi[idx]);
        if(jetVetoMap > 0) jet_vetoMap_mask[idx] = false;
        }
        return jet_vetoMap_mask;
        }

        Vec_b cleaningJetIDMask(Vec_f jet_eta, Vec_f jet_chHEF, Vec_f jet_neHEF, Vec_f jet_chEmEF, Vec_f jet_neEmEF, Vec_f jet_muEF, Vec_f jet_chMultiplicity, Vec_f jet_neMultiplicity, string year) {
        Vec_b jetID_mask(jet_eta.size(), true);
        for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {
        double jetID = corr_sf.eval_jetID(jet_eta[idx], jet_chHEF[idx], jet_neHEF[idx], jet_chEmEF[idx], jet_neEmEF[idx], jet_muEF[idx], jet_chMultiplicity[idx], jet_neMultiplicity[idx]);
        if (jetID < 0) jetID_mask[idx] = false;
        }
        return jetID_mask;
        }
        #endif
        '''
    )


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
    ##
    if((str(year) == '12022') or (str(year) == '22022')):
        loadJSON("{}/cert/Cert_Collisions2022_355100_362760_Golden.json".format(dirJson))
    if((str(year) == '12023') or (str(year) == '22023')):
        loadJSON("{}/cert/Cert_Collisions2023_366442_370790_Golden.json".format(dirJson))
    if(str(year) == '2024'):
        loadJSON("{}/cert/Cert_Collisions2024_378981_386951_Golden.json".format(dirJson))

