#!env python
import sys, os
import re
#from array import array
import math
from optparse import OptionParser,OptionGroup

parser= OptionParser()

parser.add_option("","--inputFileSig",type='string',help="Input ROOT file. [%default]", default="WS/Signal_ggHcat_2024_workspace.root")
parser.add_option("","--inputFileBKG",type='string',help="Input ROOT file bkg model. [%default]", default="WS/Bkg_ggHcat_2024_workspace.root")

parser.add_option("-c","--whichCat",type='string',help="Which category (Wcat, Zcat Zinvcatm, VBFcat)", default="Wcat")
parser.add_option("-b","--whichBIN",type='string',help="Which bin of BDT (_bin1 or empty)", default="")
parser.add_option("-y","--whichYear",type='string',help="Which year", default="2024")
parser.add_option("-o","--output",type='string',help="Output ROOT file. [%default]", default="workspace_STAT_Rho_2018.root")
parser.add_option("-d","--datCardName",type='string',help="Output txt file. [%default]", default="datacard_STAT_2024.txt")

opts, args = parser.parse_args()

sys.argv=[]
import ROOT
ROOT.gROOT.SetBatch()

ROOT.gSystem.Load("/code/HiggsAnalysis/CombinedLimit/build/lib/libHiggsAnalysisCombinedLimit.so")

doSyst=False
MultiPdf=True

############ CONFIGURABLES ###########

if MultiPdf:
    BkgPdf={
        'TTLcat': 'multipdf',
        'TTHcat': 'multipdf',
        'VLcat': 'multipdf',
        'VHcat': 'multipdf',
        'Zinvcat': 'multipdf',
        'VBFcat': 'multipdf',
        'ggHcat': 'multipdf',
    }
else:
    #OLD
    BkgPdf={
        'Vcat': 'exp1',
        'Wcat': 'exp1',
        'Zcat': 'exp1',
        'VBFcat': 'bxg',
        'Zinvcat': 'exp1',
        'VBFcatlow': 'bxg',
        'GFcat': 'bxg',
    }

SigPdf={
    'TTLcat': 'crystal_ball',
    'TTHcat': 'crystal_ball',
    'VLcat': 'crystal_ball',
    'VHcat': 'crystal_ball',
    'Zinvcat': 'crystal_ball',
    'VBFcat': 'crystal_ball',
    'ggHcat': 'crystal_ball',
}

# these are the various signal datasets: TO UPDATES
ENUM={
    'ggH': 0,
    'qqH': -1,
    'VH': -2,
    'ttH': -3,
}

#OLD
QCDscale={
    'ggH': '0.961/1.0039',
    'qqH': '0.997/1.004',
    'WH': '0.993/1.006',
    'ZH': '0.995/1.005',
    'ZinvH': '0.995/1.005',
    'WHl': '0.993/1.006',
    'ZHl': '0.995/1.005',
    'TTH': '0.995/1.005', ## TMP
}

#OLD
pdf_Higgs={
    'ggH': '0.968/1.032',
    'VBFH': '0.979/1.021',
    'WH': '0.98/1.020',
    'ZH': '0.981/1.019', #assume the majority is not ggZH
    'ZinvH': '0.981/1.019', #assume the majority is not ggZH
    'WHl': '0.98/1.020',
    'ZHl': '0.981/1.019', #assume the majority is not ggZH
    'TTH': '0.981/1.019', ## TMP
}

#OLD
lumi={
    '_2016': '1.007',
    '_2017': '1.008',
    '_2018': '1.011',
    '_12022': '1.',
    '_12022': '1.',
    '_22023': '1.',
    '_22023': '1.',
    '_2024': '1.',
}

def addSystematics(systname, systtype, value, whichProc, category, mcAll, datacard):

    datacard.write(systname+" \t"+systtype)
    for cat in category:
        print("cat",cat)
        for proc in mcAll:
            print("proc",proc)
            if (proc==whichProc): datacard.write("\t"+value)
            else:
                datacard.write("\t-")
    datacard.write("\n")

############ cat and processes ###########

category_suffix = {
    "_ggHcat":  "ggHcat",
    "_VBFcat":  "VBFcat",
    "_VHcat":   "VHcat",
    "_VLcat":   "VLcat",
    "_Zinvcat": "Zinvcat",
    "_TTHcat":  "TTHcat",
    "_TTLcat":  "TTLcat",
}

suffix = category_suffix["_"+opts.whichCat]

# Define defaults so they always exist
sigAll = []
mcAll  = []
category = []


#    n_ggH = w.var(f"{myVar}_ggH_N_SM")
#    n_qqH = w.var(f"{myVar}_qqH_N_SM")
#    n_VH  = w.var(f"{myVar}_VH_N_SM")
#    n_ttH = w.var(f"{myVar}_ttH_N_SM")

if opts.whichCat=='ggHcat':
    sigAll = ['ggH','qqH','VH','ttH']
    mcAll = ['ggH','qqH','VH','ttH','bkgGF']
    category = ['ggHcat']

if opts.whichCat=='VBFcat':
    sigAll = ['ggH','qqH','VH','ttH']
    mcAll = ['ggH','qqH','VH','ttH','bkgVBF']
    category = ['VBFcat']

if opts.whichCat=='Zinvcat':
    sigAll = ['qqH','VH','ttH']
    mcAll = ['qqH','VH','ttH','bkgV']
    category = ['Zinvcat']

if opts.whichCat=='VLcat':
    sigAll = ['VH','ttH']
    mcAll = ['VH','ttH','bkgV']
    category = ['VLcat']

if opts.whichCat=='VHcat':
    sigAll = ['VH','ttH']
    mcAll = ['VH','ttH','bkgV']
    category = ['VHcat']

if opts.whichCat=='TTLcat':
    sigAll = ['VH','ttH']
    mcAll = ['VH','ttH','bkgtth']
    category = ['TTLcat']

if opts.whichCat=='TTHcat':
    sigAll = ['VH','ttH']
    mcAll = ['ttH','bkgtth']
    category = ['TTHcat']

################### OPEN OUTPUT ############
w = ROOT.RooWorkspace("w","w")
w.Print()

############ DATACARD ###########
#datName = "cms_datacard_ws"
#datName += ".txt"
datName = opts.datCardName

numChannel = len(mcAll)-1
print("-> Opening datacard",datName)
datacard=open(datName,"w")
datacard.write("-------------------------------------\n")
datacard.write("imax "+str(len(category))+" number of channels\n")
datacard.write("jmax "+ str(numChannel)+" number of background minus 1\n")
datacard.write("kmax * number of nuisance parameters\n")
datacard.write("-------------------------------------\n")

########################## IMPORT DATA #############

xlowRange = 110
xhighRange = 150

w.factory(f"mh{suffix}[110,150]") # RooRealVar
mh=w.var(f"mh{suffix}")

arglist_obs = ROOT.RooArgList(mh)
argset_obs = ROOT.RooArgSet(mh)

################# Import SIGNAL/BKG CONTRIBUTIONS ##############

fSigIn = ROOT.TFile.Open(opts.inputFileSig,"READ")
if fSigIn == None: 
    print("ERROR: file",opts.inputFileSig,"doesn't exist")
    exit(1)    

fBkgIn = ROOT.TFile.Open(opts.inputFileBKG,"READ")
if fBkgIn == None: 
    print("ERROR: file",opts.inputFileBKG,"doesn't exist")
    exit(1)        

#Backgrounds are given a positive number, while 0 and negative numbers are used for signal processes. Different process identifiers must be used for different processes.

for cat in category:
    tag = cat+"_"+opts.whichBIN+"_"+opts.whichYear
    for proc in mcAll:
        print('doing', proc)
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow' or proc=='bkgtth':
            if opts.inputFileBKG != "":
                wInput=fBkgIn.Get("w")
                name = BkgPdf[cat]+"_"+tag+"_"+'bkg'
                nameNorm = name+"_norm"
        else:
            if opts.inputFileSig != "":
                wInput=fSigIn.Get("w")
                name = SigPdf[cat]+"_"+tag+"_"+proc
                nameNorm = name+"_norm"
                
        print("proc=",proc," cat=",cat," name=",name)
        func = wInput.pdf(name)
        if func == None: raise IOError("Unable to get func" + name)
        getattr(w,'import')(func)

        datacard.write("shapes")
        datacard.write("\t"+proc )
        datacard.write("\t"+tag)
        datacard.write("\t"+opts.output )
        datacard.write("\tw:"+name )
        datacard.write("\n")

    wInput=fBkgIn.Get("w")
    wInput.Print()
    hist_data = wInput.data("datahist_"+tag)
    hist_data.SetName("observed_data_"+tag)
    getattr(w,'import')(hist_data)
    datacard.write("shapes")
    datacard.write("\tdata_obs" )
    datacard.write("\t"+tag )
    datacard.write("\t"+opts.output )
    datacard.write("\tw:"+"observed_data_"+tag)
    datacard.write("\n")

#### OBSERVATION
datacard.write("-------------------------------------\n")
datacard.write("bin")
for cat in category:
    datacard.write("\t"+tag)
datacard.write("\n")
datacard.write("observation")
for cat in category:
    datacard.write("\t-1")
datacard.write("\n")

#### RATE
datacard.write("-------------------------------------\n")
datacard.write("bin\t\t")
#for cat in range(0,opts.ncat):
for cat in category:    
    for proc in mcAll:
        datacard.write("\t"+tag)
datacard.write("\n")
datacard.write("process\t\t")
for cat in category:    
    for proc in mcAll:
        datacard.write("\t"+proc)
datacard.write("\n")
datacard.write("process\t\t")
for cat in category:
    for idx,proc in enumerate(mcAll): ## TRICK, put the only signal first
        print(idx,proc,(-1)*idx)
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgtth' or proc=='bkgVBF': datacard.write("\t1")
        else:
#            newIDX=(-1)*idx
            newIDX=ENUM[proc]
            datacard.write("\t%d"%newIDX)
datacard.write("\n")

############ UL_r ########

datacard.write("rate\t\t")
for cat in category:
    tag = cat+"_"+opts.whichBIN+"_"+opts.whichYear
    for proc in mcAll:
        #        datacard.write("\t1")
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgtth' or proc=='bkgVBF':
            datacard.write("\t1") # I handle this via rateParma
        else:
            wInput=fSigIn.Get("w")
            name = SigPdf[cat]+"_"+tag+"_"+proc
            var = wInput.obj(name+"_norm")
            datacard.write("\t%6f"%var.getVal())
datacard.write("\n")

############ SYST ########
datacard.write("-------------------------------------\n")


## rateParam: take the norm from the workspace
for cat in category:
    tag = cat+"_"+opts.whichBIN+"_"+opts.whichYear
    for proc in mcAll:
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgtth' or proc=='bkgVBF':
            wInput=fBkgIn.Get("w")
            name = BkgPdf[cat]+"_"+tag+"_"+'bkg'
            nameNorm = name+"_norm"

            datacard.write(nameNorm)
            datacard.write("\trateParam")
            datacard.write("\t")
            datacard.write(tag)
            datacard.write("\t")
            datacard.write(proc)
            datacard.write("\t")
            var = wInput.var(nameNorm)
            datacard.write("\t%.6f" % var.getVal())
datacard.write("\n")

####

datacard.write("lumi_13p6TeV \tlnN ")

for cat in category:
    for proc in mcAll:
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgtth': datacard.write("\t-")
        else:
            if opts.whichCat=='GFcat' or opts.whichCat=='VBFcatlow': datacard.write("\t%.3f"%(1.025) )
            else: datacard.write("\t%.3f"%(1.016) )
    datacard.write("\n")


    datacard.write("BR_mm   \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgtth': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.012) )
    datacard.write("\n")

if doSyst:

    datacard.write("CMS_pileup  \tlnN ")
    for cat in category:
        for proc in mcAll:
            datacard.write("\t-")            
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgtth': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.01) )
    datacard.write("\n")

#    ### Add autoMCStats
#    datacard.write("\n")
#    datacard.write("* autoMCStats 0\n")

datacard.write("-------------------------------------\n")

if MultiPdf:
    tag = cat+"_"+opts.whichBIN+"_"+opts.whichYear
    pdfindexSTR='pdfindex_'+tag
    datacard.write("\n")
    datacard.write("%s discrete\n"%pdfindexSTR)

############ DONE ########
w.writeToFile(opts.output)
print("->Done")
