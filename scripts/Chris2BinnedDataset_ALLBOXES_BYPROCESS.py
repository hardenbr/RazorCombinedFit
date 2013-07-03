from optparse import OptionParser
import os
from array import array
from pdfShit import *

import ROOT as rt
import RootTools
from RazorCombinedFit.Framework import Config

boxMap = {'MuEle':0,'MuMu':1,'EleEle':2,'Mu':3,'Ele':4,'Had':5, 'BJet':6}
cross_sections = {'SingleTop_s':4.21,'SingleTop_t':64.6,'SingleTop_tw':10.6,\
                               'TTj':157.5,'Zll':3048,'Znn':2*3048,'Wln':31314,\
                               'WW':43,'WZ':18.2,'ZZ':5.9,'Vgamma':173
                               }
lumi = 1.0

def isInFitRegion(x, y, box, config, doMultijet):
    myworkspace = rt.RooWorkspace(box)
    myvariables = config.getVariablesRange(box,"variables",myworkspace)
    myargs = myworkspace.allVars()
    
    def testInRegion(region, x, y):
        inMR = x >= myargs['MR'].getMin(region) and x < myargs['MR'].getMax(region)
        inRsq = y >= myargs['Rsq'].getMin(region) and y < myargs['Rsq'].getMax(region)
        return inMR and inRsq

    regions = ['fR1','fR2','fR3','fR4']
    if doMultijet:
        regions.append('fR5')
    
    result = False
    for r in regions:
        if testInRegion(r, x, y):
            result = True
            break

    del myvariables
    del myworkspace
    return result
    
def cutFitRegion(histo, box, config, doMultijet):

    iBinx = histo.GetNbinsX()
    iBiny = histo.GetNbinsY()
    xAxis = histo.GetXaxis()
    yAxis = histo.GetYaxis()
    for ix in range(1, iBinx+1):
        for iy in range(1, iBiny+1):
            if isInFitRegion(xAxis.GetBinCenter(ix),yAxis.GetBinCenter(iy),box,config,doMultijet): histo.SetBinContent(ix,iy,0.)
    return histo

def getMeanSigma(n0, nP, nM):
    maxVal = max(n0, nP, nM)
    minVal = min(n0, nP, nM)
    return (maxVal+minVal)/2., (maxVal-minVal)/2.

def writeTree2DataSet(data, outputFile, outputBox, rMin, mRmin):
    
    output = rt.TFile.Open(outputFile+"_MR"+str(mRmin)+"_R"+str(rMin)+'_'+outputBox,'UPDATE')
    for mydata in data:
        mydata.Write()
    print 'Writing',output.GetName()
    output.Close()

def convertTree2Dataset(tree, outputFile, config, minH, maxH, btag, nToys, varBin, doSMS, doXsec, doMultijet, write = True):
    """This defines the format of the RooDataSet"""
    boxes = config.getBoxes()

    wHisto = []
    wHisto_JESup = []
    wHisto_JESdown = []
    #wHisto_xsecup = []
    #wHisto_xsecdown = []
    wHisto_pdfCEN = []
    wHisto_pdfSYS = []
    wHisto_btagup = []
    wHisto_btagdown = []    

    plotByProcess = []
    plotByProcessJESUP = []
    plotByProcessJESDOWN = []
    #plotByProcessXSECUP = []
    #plotByProcessXSECDOWN = []
    plotByProcessPDFCEN = []
    plotByProcessPDFSYS = []
    
    plotByProcessBTAGUP = []
    plotByProcessBTAGDOWN = []

    yieldByProcess = []

    # define the loosest bin ranges
    loose_bin = 'MuMu'
    if loose_bin not in boxes:
        loose_bin = 'Had'
    workspace = rt.RooWorkspace(loose_bin)
    variables = config.getVariablesRange(loose_bin,"variables",workspace)
    args = workspace.allVars()
                    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()
    
    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()
    rMin = rt.TMath.Sqrt(rsqMin)
    rMax = rt.TMath.Sqrt(rsqMax)

    binedgexLIST = []
    binedgeyLIST = []
    #either use a binning scheme defined here or take from the config
    if not config.hasBinning(loose_bin):
        # if the bin is fixed, do 50 GeV in mR
        # and 0.1 in R^2
        binwMR = 50.
        binwR2 = 0.1
        #use a fixed bin for mR
        if varBin != 1: maxVal = mRmax
        else: maxVal = 700.
        mRedge = mRmin
        while mRedge < maxVal: 
            binedgexLIST.append(mRedge)
            mRedge = mRedge + binwMR
        binedgexLIST.append(maxVal)
        if varBin == 1:
            if mRmax> 800: binedgexLIST.append(800)
            if mRmax> 900: binedgexLIST.append(900)
            if mRmax> 1000: binedgexLIST.append(1000)
            if mRmax> 1200: binedgexLIST.append(1200)
            if mRmax> 1600: binedgexLIST.append(1600)
            if mRmax> 2000: binedgexLIST.append(2000)
            if mRmax> 2800: binedgexLIST.append(2800)
            binedgexLIST.append(mRmax)

            #use a fixed bin for R^2
            if varBin != 1:
                R2edge = rsqMin
                while R2edge <rsqMax: 
                    binedgexLIST.append(R2edge)
                    R2edge = R2edge + binwR2
                binedgeyLIST.append(rsqMax)
            else: 
                #use fixed binning 
                binedgeyLIST = [rsqMin,0.18,0.2,0.3,0.4,0.5]
    else:
        #take from the config
        binning = config.getBinning(loose_bin)
        #MR and Rsq in that order
        binedgexLIST.extend(binning[0])
        binedgeyLIST.extend(binning[1])

    nbinx =  len(binedgexLIST)-1
    nbiny = len(binedgeyLIST)-1    

    binedgex = array('d',binedgexLIST)
    binedgey = array('d',binedgeyLIST)

    if doSMS == True:
        nProcess   = 1
        weight     = "W_EFF"
        #weightUP   = "W_EFF"
        #weightDOWN = "W_EFF"
    else:
        nProcess = 10
        if doXsec == 1:    weight = "W_UP"
        elif doXsec == -1: weight = "W_DOWN"
        else: weight = "W"

    # get weight
    for box in boxes:
        yieldByProcessThisBox = []        
        for i in range(0,nProcess): 
            histo = rt.TH2D("wHisto_TMP","wHisto_TMP", nbinx, binedgex, nbiny, binedgey) 
            #tree.Project("wHisto_TMP", "RSQ:MR", 'LEP_W*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BOX_NUM == %i) && (BTAG_NUM>=%i))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,boxMap[box], btag))
            tree.Project("wHisto_TMP", "RSQ:MR", 'LEP_W*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BOX_NUM == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,boxMap[box], btag,btag))
            yieldByProcessThisBox.append(histo.Integral())
            del histo
        yieldByProcess.append([box,yieldByProcessThisBox])

    # plots by process TO CHECK
    for i in range(0,nProcess): 
        
        # nominal
        histo = rt.TH2D("wHisto_PROC%i" %i,"wHisto_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
        tree.Project("wHisto_PROC%i" %i, "RSQ:MR", 'LEP_W*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
        plotByProcess.append(histo.Clone())

        # JES UP
        histo = rt.TH2D("wHisto_JESup_PROC%i" %i,"wHisto_JESup_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
        tree.Project("wHisto_JESup_PROC%i" %i, "RSQ_JES_UP:MR_JES_UP", 'LEP_W*%s*(MR_JES_UP >= %f && MR_JES_UP <= %f && RSQ_JES_UP >= %f && RSQ_JES_UP <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
        plotByProcessJESUP.append(histo.Clone())

        # JES DOWN
        histo = rt.TH2D("wHisto_JESdown_PROC%i" %i,"wHisto_JESdown_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
        tree.Project("wHisto_JESdown_PROC%i" %i, "RSQ_JES_DOWN:MR_JES_DOWN", 'LEP_W*%s*(MR_JES_DOWN >= %f && MR_JES_DOWN <= %f && RSQ_JES_DOWN >= %f && RSQ_JES_DOWN <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
        plotByProcessJESDOWN.append(histo.Clone())
        
    
        if doMultijet:
            #BTag up
            histo = rt.TH2D("wHisto_btagup_PROC%i" %i,"wHisto_btagup_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
            tree.Project("wHisto_btagup_PROC%i" %i,  "RSQ:MR",  '(LEP_W+LEP_W_SYS)*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
            plotByProcessBTAGUP.append(histo.Clone())            
            
            #BTag down
            histo = rt.TH2D("wHisto_btagdown_PROC%i" %i,"wHisto_btagdown_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
            tree.Project("wHisto_btagdown_PROC%i" %i,  "RSQ:MR",  '(LEP_W-LEP_W_SYS)*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weight,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
            plotByProcessBTAGDOWN.append(histo.Clone())            

        # XSEC UP
        #histo = rt.TH2D("wHisto_xsecup_PROC%i" %i,"wHisto_xsecup_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
        #tree.Project("wHisto_xsecup_PROC%i" %i,  "RSQ:MR",  'LEP_W*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weightUP,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
        #plotByProcessXSECUP.append(histo.Clone())

        # XSEC DOWN
        #histo = rt.TH2D("wHisto_xsecdown_PROC%i" %i,"wHisto_xsecdown_PROC%i" %i, nbinx, binedgex, nbiny, binedgey)
        #tree.Project("wHisto_xsecdown_PROC%i" %i,  "RSQ:MR",  'LEP_W*%s*(MR >= %f && MR <= %f && RSQ >= %f && RSQ <= %f && (PROC == %i) && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f))' % (weightDOWN,mRmin,mRmax,rsqMin,rsqMax,i,btag,btag))
        #plotByProcessXSECDOWN.append(histo.Clone())

        # PDF CEN and SIGMA
        histo_pdfCEN, histo_pdfSYS = makePDFPlotCONDARRAY(tree, histo, nbinx, binedgex, nbiny, binedgey, "PROC == %i && (BTAG_TrackCount[0]>=%f || BTAG_TrackCount[1]>=%f)" %(i,btag,btag))
        histo_pdfCEN.SetName("wHisto_pdfCEN_PROC%i" %i)
        histo_pdfSYS.SetName("wHisto_pdfSYS_PROC%i" %i)
        plotByProcessPDFCEN.append(histo_pdfCEN.Clone())
        plotByProcessPDFSYS.append(histo_pdfSYS.Clone())

    #this is the nominal histogram
    for box in boxes:               
        histo = rt.TH2D("wHisto_%s" %box,"wHisto_%s" %box, nbinx, binedgex, nbiny, binedgey)
        for yieldByProcessThisBox in yieldByProcess:
            if yieldByProcessThisBox[0] != box: continue
            for iProc in range(0,len(yieldByProcessThisBox[1])):
                thishisto = plotByProcess[iProc]
                thisYield = yieldByProcessThisBox[1][iProc]
                if thishisto.Integral() != 0. : histo.Add(thishisto,thisYield/thishisto.Integral())
        wHisto.append(histo.Clone())

        # JES correctiobns UP
        histo_JESup = rt.TH2D("wHisto_JESup_%s" %box,"wHisto_JESup_%s" %box, nbinx, binedgex, nbiny, binedgey)
        for yieldByProcessThisBox in yieldByProcess:
            if yieldByProcessThisBox[0] != box: continue
            for iProc in range(0,len(yieldByProcessThisBox[1])):
                thishisto = plotByProcessJESUP[iProc]
                thisYield = yieldByProcessThisBox[1][iProc]
                if thishisto.Integral() != 0. : histo_JESup.Add(thishisto,thisYield/thishisto.Integral())
        wHisto_JESup.append(histo_JESup.Clone())

        # JES correctiobns DOWN
        histo_JESdown = rt.TH2D("wHisto_JESdown_%s" %box,"wHisto_JESdown_%s" %box, nbinx, binedgex, nbiny, binedgey)
        for yieldByProcessThisBox in yieldByProcess:
            if yieldByProcessThisBox[0] != box: continue
            for iProc in range(0,len(yieldByProcessThisBox[1])):
                thishisto = plotByProcessJESDOWN[iProc]
                thisYield = yieldByProcessThisBox[1][iProc]
                if thishisto.Integral() != 0. : histo_JESdown.Add(thishisto,thisYield/thishisto.Integral())
        wHisto_JESdown.append(histo_JESdown.Clone())

        if doMultijet:
            
            #Btag UP
            histo_btagup = rt.TH2D("wHisto_btagup_%s" %box,"wHisto_btagup_%s" %box, nbinx, binedgex, nbiny, binedgey)
            for yieldByProcessThisBox in yieldByProcess:
                if yieldByProcessThisBox[0] != box: continue
                for iProc in range(0,len(yieldByProcessThisBox[1])):
                    thishisto = plotByProcessBTAGUP[iProc]
                    thisYield = yieldByProcessThisBox[1][iProc]
                    if thishisto.Integral() != 0. : histo_btagup.Add(thishisto,thisYield/thishisto.Integral())                
            wHisto_btagup.append(histo_btagup.Clone())

            #Btag down            
            histo_btagdown = rt.TH2D("wHisto_btagdown_%s" %box,"wHisto_btagdown_%s" %box, nbinx, binedgex, nbiny, binedgey)
            for yieldByProcessThisBox in yieldByProcess:
                if yieldByProcessThisBox[0] != box: continue
                for iProc in range(0,len(yieldByProcessThisBox[1])):
                    thishisto = plotByProcessBTAGDOWN[iProc]
                    thisYield = yieldByProcessThisBox[1][iProc]
                    if thishisto.Integral() != 0. : histo_btagdown.Add(thishisto,thisYield/thishisto.Integral())                
            wHisto_btagdown.append(histo_btagdown.Clone())
            
            

        # xsec UP
        #histo_xsecup = rt.TH2D("wHisto_xsecup_%s" %box,"wHisto_xsecup_%s" %box, nbinx, binedgex, nbiny, binedgey)
        #for yieldByProcessThisBox in yieldByProcess:
        #    if yieldByProcessThisBox[0] != box: continue
        #    for iProc in range(0,len(yieldByProcessThisBox[1])):
        #        thishisto = plotByProcessXSECUP[iProc]
        #        thisYield = yieldByProcessThisBox[1][iProc]
        #        if thishisto.Integral() != 0. : histo_xsecup.Add(thishisto,thisYield/thishisto.Integral())                
        #wHisto_xsecup.append(histo_xsecup.Clone())
        
        # xsec correctiobns DOWN
        #histo_xsecdown = rt.TH2D("wHisto_xsecdown_%s" %box,"wHisto_xsecdown_%s" %box, nbinx, binedgex, nbiny, binedgey)
        #for yieldByProcessThisBox in yieldByProcess:
        #    if yieldByProcessThisBox[0] != box: continue
        #    for iProc in range(0,len(yieldByProcessThisBox[1])):
        #        thishisto = plotByProcessXSECDOWN[iProc]
        #        thisYield = yieldByProcessThisBox[1][iProc]
        #        if thishisto.Integral() != 0. : histo_xsecdown.Add(thishisto,thisYield/thishisto.Integral())   
        #wHisto_xsecdown.append(histo_xsecdown.Clone())
        
        # PDF central (new nominal) and error (for systematics)
        histo_pdfCEN = rt.TH2D("wHisto_pdfCEN_%s" %box,"wHisto_pdfCEN_%s" %box, nbinx, binedgex, nbiny, binedgey)
        for yieldByProcessThisBox in yieldByProcess:
            if yieldByProcessThisBox[0] != box: continue
            for iProc in range(0,len(yieldByProcessThisBox[1])):
                thishisto = plotByProcessPDFCEN[iProc]
                thisYield = yieldByProcessThisBox[1][iProc]
                if thishisto.Integral() != 0. : histo_pdfCEN.Add(thishisto,thisYield/thishisto.Integral()) 
        # to compute the error, we need some gymnastic
        histo_pdfSYS = rt.TH2D("wHisto_pdfSYS_%s" %box,"wHisto_pdfSYS_%s" %box, nbinx, binedgex, nbiny, binedgey)
        for ix in range(1,nbinx+1):
            for iy in range(1,nbiny+1):
                error = 0.
                for iProc in range(0,len(yieldByProcessThisBox[1])):
                    thisYield = yieldByProcessThisBox[1][iProc]
                    thishisto = plotByProcessPDFCEN[iProc]                                    
                    if thishisto.Integral() != 0. : error = error + pow(plotByProcessPDFSYS[iProc].GetBinContent(ix,iy)*thisYield/thishisto.Integral(),2.)
                histo_pdfSYS.SetBinContent(ix,iy,math.sqrt(error))
        wHisto_pdfCEN.append(histo_pdfCEN.Clone())
        wHisto_pdfSYS.append(histo_pdfSYS.Clone())

    del workspace
    del variables
    del args
    del mRmin
    del mRmax

    # random number generator
    pid = os.getpid()
    now = rt.TDatime()
    today = now.GetDate()
    clock = now.GetTime()
    seed = today+clock+pid+137
    gRnd = rt.TRandom3(seed)

    for i in xrange(nToys):
        # correlated systematics: LUMI 2.2% MULTIPLICATIVE sumInQuadrature  sumInQuadrature RvsMR trigger 2% = 4.9%
        lumiFactor = math.pow((1.022), gRnd.Gaus(0., 1.))
        # triggerLepton 3% per trigger set
        muTriggerFactor =  math.pow(1.03,gRnd.Gaus(0.,1.))
        eleTriggerFactor =  math.pow(1.03,gRnd.Gaus(0.,1.))
        #take 5% for the trigger being poorly emulated in MC
        dataMCTriggerFactor =  math.pow(1.05,gRnd.Gaus(0.,1.))
        btagFactor = math.pow(1.06,gRnd.Gaus(0.,1.))
        #up to 30% of events have a gen level tau, where the 
        tauFactor = math.pow(1.02,gRnd.Gaus(0.,1.))
        # correlated systematics: Btag scale factor corrections ADDITIVE (scaled bin by bin)
        btagSFFactor  = gRnd.Gaus(0., 1.) 
        # correlated systematics: xsection ADDITIVE (scaled bin by bin)
        #xsecFactor = gRnd.Gaus(0., 1.)
        for ibox in range(0,len(boxes)):
            box = boxes[ibox]
            workspace = rt.RooWorkspace(box)
            variables = config.getVariablesRange(box,"variables",workspace)
            args = workspace.allVars()

            #we cut away events outside our MR window
            mRmin = args['MR'].getMin()
            mRmax = args['MR'].getMax()

            #we cut away events outside our Rsq window
            rsqMin = args['Rsq'].getMin()
            rsqMax = args['Rsq'].getMax()
            rMin = rt.TMath.Sqrt(rsqMin)
            rMax = rt.TMath.Sqrt(rsqMax)

            # set the same binning for the RooRealVars
            MR =  workspace.var("MR")
            Rsq =  workspace.var("Rsq")

            box = boxes[ibox]

            #write the nominal only once
            if i == 0: 
                histo = cutFitRegion(wHisto[ibox],box,config,doMultijet)
                data = [histo.Clone()]
                rooDataHist = rt.RooDataHist("RMRHistTree_%s" %box,"RMRHistTree_%s" %box,rt.RooArgList(rt.RooArgSet(MR,Rsq)),histo)
                data.append(rooDataHist)
                data.append(wHisto_JESup[ibox])
                if write: writeTree2DataSet(data, outputFile, "%s.root" %box, rMin, mRmin)
                
            # create a copy of the histogram
            wHisto_i = rt.TH2D("wHisto_%s_%i" %(box, i),"wHisto_%s_%i" %(box, i), nbinx, binedgex, nbiny, binedgey)
            for ix in range(1,nbinx+1):
                for iy in range(1,nbiny+1):
                    # uncorrelated systematics: lepton efficiency data/MC 1%
                    lepFactor = math.pow(1.01,gRnd.Gaus(0., 1.))
                    # uncorrelated systematics: JES corrections ADDITIVE (scaled bin by bin)
                    jesFactor  = gRnd.Gaus(0., 1.)
                    # compute the total
                    # starting value
                    nominal = wHisto[ibox].GetBinContent(ix,iy)
                    if nominal != 0:
                        # add lumi systematics
                        newvalue = nominal*lumiFactor
                        if doMultijet: newvalue = newvalue*dataMCTriggerFactor
                        if not doMultijet and box == "BJet" or btag> 0.0: newvalue = newvalue*btagFactor
                        # add the lep trigger eff
                        if box == "MuMu" or box == "MuEle" or box == "Mu": newvalue = newvalue*muTriggerFactor
                        if box == "EleEle" or box == "Ele": newvalue = newvalue*eleTriggerFactor
                        #ignore leptons for the hadronic boxes
                        if box not in ['Had','BJet']: newvalue = newvalue*lepFactor
                        if doMultijet:
                            newvalue = newvalue*tauFactor
                        if doMultijet:
                            # add btag systematics
                            mBtag, sBtag = getMeanSigma(nominal, wHisto_btagup[ibox].GetBinContent(ix,iy), wHisto_btagdown[ibox].GetBinContent(ix,iy))
                            if mBtag > 0: newvalue = newvalue*mBtag/nominal*math.pow(1.+sBtag/mBtag, btagSFFactor)
                        
                        # add xsec systematics
                        #if xsecFactor > 0: newvalue = newvalue* + xsecFactor*(wHisto_xsecup[ibox].GetBinContent(ix,iy)-nominal)
                        #else: newvalue = newvalue*math.pow( + xsecFactor*(wHisto_xsecdown[ibox].GetBinContent(ix,iy)-nominal))
                        #mXsec, sXsec = getMeanSigma(nominal, wHisto_xsecup[ibox].GetBinContent(ix,iy), wHisto_xsecdown[ibox].GetBinContent(ix,iy))
                        #if mXsec > 0: newvalue = newvalue*mXsec/nominal*math.pow(1.+sXsec/mXsec, xsecFactor)
                        # add jes systematics
                        #if jesFactor > 0: newvalue = newvalue + jesFactor*(wHisto_JESup[ibox].GetBinContent(ix,iy)-nominal)
                        #else: newvalue = newvalue + jesFactor*(wHisto_JESdown[ibox].GetBinContent(ix,iy)-nominal)               
                        mJES, sJES = getMeanSigma(nominal, wHisto_JESup[ibox].GetBinContent(ix,iy), wHisto_JESdown[ibox].GetBinContent(ix,iy))
                        if mJES > 0: newvalue = newvalue*mJES/nominal*math.pow(1.+sJES/mJES, jesFactor)
                        # apply the systematic correction to the pdf value
                        # the pdf code return the efficiency in each bin, with an error
                        # that includes the systematic effect. We use this to get a
                        # new value for the content of the bin
                        mPDF = wHisto_pdfCEN[ibox].GetBinContent(ix,iy)
                        sPDF = wHisto_pdfSYS[ibox].GetBinContent(ix,iy)
                        if 1+mPDF > 0.: newvalue = newvalue*(1+mPDF)*math.pow(1+sPDF/(1+mPDF),gRnd.Gaus(0.,1.))
                        #newvalue = newvalue *(1+ gRnd.Gaus(wHisto_pdfCEN[ibox].GetBinContent(ix,iy), wHisto_pdfSYS[ibox].GetBinContent(ix,iy)))
                        # add a 20% systematics due to filling procedure
                        #byProcFactor = math.pow(1.50,gRnd.Gaus(0., 1.))
                        byProcFactor = 1
                        # fill histogram
                        wHisto_i.SetBinContent(ix,iy,max(0.,newvalue*byProcFactor))
                    else:
                        wHisto_i.SetBinContent(ix,iy,max(0.,nominal))
            # check the fit-region edge
            wHisto_i = cutFitRegion(wHisto_i,box,config,doMultijet)

            data = [wHisto_i,rt.RooDataHist("RMRHistTree_%s_%i" %(box,i),"RMRHistTree_%s_%i" %(box,i),rt.RooArgList(rt.RooArgSet(MR,Rsq)),wHisto_i)]            
            if write: writeTree2DataSet(data, outputFile, "%s.root" %box, rMin, mRmin)
            del wHisto_i
            del data
            del workspace
            del variables
            del args
            del mRmin
            del mRmax
            del MR
            del Rsq

    return []


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default=None,
                  help="Name of the config file to use")
    parser.add_option('--max',dest="max",type="int",default=-1,
                  help="The last event to take from the input Dataset")
    parser.add_option('--min',dest="min",type="int",default=0,
                  help="The first event to take from the input Dataset")  
    parser.add_option('-e','--eff',dest="eff",default=False,action='store_true',
                  help="Calculate the MC efficiencies")
    parser.add_option('-f','--flavour',dest="flavour",default='TTj',
                  help="The flavour of MC used as input")
    parser.add_option('-r','--run',dest="run",default=-1,type=float,
                  help="The minimum run number")
    parser.add_option('-d','--dir',dest="outdir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-t','--toys',dest="toys",type="int",default=1000,
                  help="Number of toys")
    parser.add_option('-v','--VariableBinning', dest="varbin", type="int",default=1,
                                        help="Use Variable Binning for mR")
    parser.add_option('--sms',dest="doSMS",action="store_true", default=False,
                      help="Run PDF filling for SMS models")
    parser.add_option('--multijet',dest="doMultijet",action="store_true", default=False,
                      help="We are in the multijet analysis")    
    parser.add_option('-b','--btag',dest="btag",type="int",default=-1,
                  help="The minimum number of Btags to allow")     
    parser.add_option('--xsecUp',dest="xsecUp",action="store_true", default=False,
                      help="Run with theory cross section up by one sigma")
    parser.add_option('--xsecDown',dest="xsecDown",action="store_true", default=False,
                      help="Run with theory cross section down by one sigma")
    parser.add_option('--tree_name',dest="tree_name",type="string",default='EVENTS',
                  help="The name of the TTree to look at")
    
    (options,args) = parser.parse_args()

    doXsec = 0
    if options.xsecUp:   doXsec = 1
    if options.xsecDown: doXsec = -1
    
    if options.config is None:
        import inspect, os
        topDir = os.path.abspath(os.path.dirname(inspect.getsourcefile(convertTree2Dataset)))
        options.config = os.path.join(topDir,'boxConfig.cfg')    
    cfg = Config.Config(options.config)
    
    print 'Input files are %s' % ', '.join(args)
    for f in args:
        if f.lower().endswith('.root'):
            input = rt.TFile.Open(f)

            decorator = options.outdir+"/"+os.path.basename(f)[:-5]
            btag = -10000000000000.
            if options.btag > 0:
                # we do only one btag
                btag = 3.3
            convertTree2Dataset(input.Get(options.tree_name), decorator, cfg,options.min,options.max,btag,options.toys,options.varbin,options.doSMS,doXsec,options.doMultijet)

        else:
            "File '%s' of unknown type. Looking for .root files only" % f
