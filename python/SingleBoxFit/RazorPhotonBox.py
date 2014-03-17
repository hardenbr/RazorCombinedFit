from RazorCombinedFit.Framework import Box
import RootTools
import RazorBox
import ROOT as rt
from makeBluePlot import *
from array import *
#this is global, to be reused in the plot making
#def getBinning(boxName, varName):
#    if varName == "MR" : return [110, 150.,...]
#    if varName == "Rsq" : return [0.03,0.23,0.28,0.32,0.37,0.41,0.46,0.51,0.57,0.63,0.70,0.80,1.5]

class RazorPhotonBox(RazorBox.RazorBox):
    
    def __init__(self, name, variables):
        super(RazorPhotonBox,self).__init__(name, variables)
        
        # now we switch off the redundant Znn component in the Had box
        self.zeros = {'QCD':['Mu','Ele','MuMu','EleEle','MuEle']}
        self.cut = 'MR >= 0.0'
        self.fitMode = "2D"
        self.btag = "NoBtag"

    def switchOff(self, species) :
        print "NAME OF SPECIIES" + str(species)
        self.workspace.var("Ntot_"+species).setVal(0.)
        self.workspace.var("Ntot_"+species).setConstant(rt.kTRUE)
        self.workspace.var("f2_"+species).setConstant(rt.kTRUE)

    #add penalty terms and float
    def float1stComponent(self,flavour):
#        self.fixParsExact("MR0_%s" % flavour, False)
#        self.fixParsExact("R0_%s" % flavour, False)
        self.fixParsExact("b_%s" % flavour, False)
        self.fixParsExact("n_%s" % flavour, False)
#        self.fixParsExact("RI_%s" % flavour, False)        

    def float2ndComponent(self,flavour):
        self.fixParsExact("MR02nd_%s" % flavour, False)
        self.fixParsExact("R02nd_%s" % flavour, False)
        self.fixParsExact("b2nd_%s" % flavour, False)
        self.fixParsExact("n2nd_%s" % flavour, False)        

    def float1stComponentWithPenalty(self,flavour, alsoB):
        self.fixParsPenalty("MR01st_%s" % flavour)
        self.fixParsPenalty("R01st_%s" % flavour)
        self.fixPars("MR01st_%s_s" % flavour)
        self.fixPars("R01st_%s_s" % flavour)

        if alsoB == True:
            self.fixParsPenalty("b1st_%s" % flavour)
            self.fixPars("b1st_%s_s" % flavour)
#            self.fixParsPenalty("n1st_%s" % flavour)
#            self.fixPars("n1st_%s_s" % flavour)
#            self.fixParsExact("n_%s" % flavour, False)
#            self.fixParsPenalty("n1st_%s" % flavour)
#            self.fixPars("n1st_%s_s" % flavour)
        else:
            self.fixParsExact("b_%s" % flavour, False)
            self.fixParsExact("n_%s" % flavour, False)

    def float2ndComponentWithPenalty(self,flavour, alsoB):
        self.fixParsPenalty("MR02nd_%s" % flavour)
        self.fixParsPenalty("R02nd_%s" % flavour)
        self.fixPars("MR02nd_%s_s" % flavour)
        self.fixPars("R02nd_%s_s" % flavour)
#        self.fixParsExact("n2nd_%s" % flavour, False)        

        if alsoB == True:
            self.fixParsPenalty("b2nd_%s" % flavour)
            self.fixPars("b2nd_%s_s" % flavour)
#            self.fixParsExact("n_%s" % flavour, False)
            self.fixParsPenalty("n2nd_%s" % flavour)
            self.fixPars("n2nd_%s_s" % flavour)
        else:
            self.fixParsExact("b_%s" % flavour, False)
            self.fixParsExact("n_%s" % flavour, False)
            

    def floatYield(self,flavour):
        self.fixParsExact("Ntot_%s" % flavour, False)
#        self.fixParsExact("n_%s" % flavour, False)
        
    def fixComponent(self,flavour):
        self.fixParsExact("MR0_%s" % flavour, True)
        self.fixParsExact("R0_%s" % flavour, True)
        self.fixParsExact("b_%s" % flavour, True)        
        self.fixParsExact("n_%s" % flavour, True)        
    
    def define(self, inputFile):
        
        #define the ranges
        mR  = self.workspace.var("MR")
        Rsq = self.workspace.var("Rsq")
        
        # add the different components:
        # - QCD
        self.addTailPdf("QCD",True)

        # build the total PDF
        #myPDFlist = rt.RooArgList(self.workspace.pdf("ePDF1st_QCD"), self.workspace.pdf("ePDF2nd_QCD"))
        myPDFlist = rt.RooArgList(self.workspace.pdf("ePDF_QCD"))
        model = rt.RooAddPdf(self.fitmodel, self.fitmodel, myPDFlist)        
        
        # import the model in the workspace.
        self.importToWS(model)
        #print the workspace
        self.workspace.Print()

        # fix all pdf parameters to the initial value
        self.fixPars("QCD")

        def floatSomething(z):
            """Switch on or off whatever you want here"""
            # to float the parameters of the 1st component
            self.float1stComponent(z)
#            # to float the parameters of the 2nd component
#            self.float2ndComponent(z) 
            # OR 
            # to float the parameters of the 1st component with penalty terms            
#            self.float1stComponentWithPenalty(z, True)
#            self.float2ndComponentWithPenalty(z, True)
#            self.floatComponent(z)
            # float the total yield 
            self.floatYield(z)
            # float the 2nd-component fraction
#            self.floatFraction(z)
#            self.floatFractionWithPenalty(z)

        # switch off not-needed components (box by box)
        fixed = []
        for z in self.zeros:
            if self.name in self.zeros[z]:
                #floatSomething(z)
                self.fixPars(z)
                self.switchOff(z)
            else:
                if not z in fixed:
                    floatSomething(z)
                    fixed.append(z)

    def plot(self, inputFile,store, box, data=None,fitmodel=None,frName=None):
        [store.store(s, dir=box) for s in self.plot1DHistoAllComponents(inputFile, "MR", 80, ranges=['FULL'],data=data)]
        [store.store(s, dir=box) for s in self.plot1DHistoAllComponents(inputFile, "Rsq", 30, ranges=['FULL'],data=data)]

    def plot1DHistoAllComponents(self, inputFile, xvarname, nbins = 25, ranges=None, data = None):
        rangeNone = False
        if ranges is None:
            rangeNone = True
            ranges = ['']
            
        #before I find a better way
        rangeCut = self.getVarRangeCutNamed(ranges=ranges)
        if data is None:
            data = RootTools.getDataSet(inputFile,'RMRTree', self.cut)
            data = data.reduce(rangeCut)

        #GENERATE 50 times the statistics in the form of toys
        toyData = self.workspace.pdf(self.fitmodel).generate(self.workspace.set('variables'), 50*data.numEntries())

        print toyData
        print ranges
        
        toyData = toyData.reduce(self.getVarRangeCutNamed(ranges=ranges))

        xmin = min([self.workspace.var(xvarname).getMin(r) for r in ranges])
        xmax = max([self.workspace.var(xvarname).getMax(r) for r in ranges])
        
        # define 1D histograms
        bintuple = Binning("Had", True)
        
        rsq_bins = bintuple[1]
        print rsq_bins
        rsq_array = array("d",rsq_bins)

        mr_bins = bintuple[0]
        print mr_bins
        mr_array = array("d",mr_bins)

        

        if xmin < .4: #THIS IS AN RSQ
            histoData = rt.TH1D("histoData", "histoData",len(rsq_bins)-1,rsq_array)
            histoToy = rt.TH1D("histoToy", "histoToy",len(rsq_bins)-1, rsq_array)
        else: #THIS IS AN MR HIST
            histoData = rt.TH1D("histoData", "histoData",len(mr_bins)-1,mr_array)
            histoToy = rt.TH1D("histoToy", "histoToy",len(mr_bins)-1, mr_array)

        
        def setName(h, name):
            h.SetName('%s_%s_%s_ALLCOMPONENTS' % (h.GetName(),name,'_'.join(ranges)) )
            h.GetXaxis().SetTitle(name)
            # x axis
            if name == "MR": h.GetXaxis().SetTitle("M_{R} [GeV]")
            elif name == "Rsq": h.GetXaxis().SetTitle("R^{2}")
            # y axis
            if name == "MR": h.GetYaxis().SetTitle("N Events")
            elif name == "Rsq": h.GetYaxis().SetTitle("N Events")
            # axis labels
            h.GetXaxis().SetTitleSize(0.05)
            h.GetYaxis().SetTitleSize(0.05)
            h.GetXaxis().SetLabelSize(0.05)
            h.GetYaxis().SetLabelSize(0.05)
            h.GetXaxis().SetTitleOffset(0.90)
            h.GetYaxis().SetTitleOffset(0.93)
                                                                                                    
            
        def SetErrors(histo, nbins):
            for i in range(1, nbins+1):
                histo.SetBinError(i,rt.TMath.Sqrt(histo.GetBinContent(i)))

        # project the data on the histograms
        #data.tree().Project("histoData",xvarname)
        data.fillHistogram(histoData,rt.RooArgList(self.workspace.var(xvarname)))
        toyData.fillHistogram(histoToy,rt.RooArgList(self.workspace.var(xvarname)))
        print "histData " + str(histoData.Integral())
        print "histToy " + str(histoToy.Integral())
        toyData.Print()
        if histoToy.Integral() > 0:
            scaleFactor = histoData.Integral()/histoToy.Integral()
        else:
            scaleFactor = 1
            print "no toys! Integral = 0"        

        histoToy.Scale(scaleFactor)
        SetErrors(histoToy, nbins)
        setName(histoData,xvarname)
        setName(histoToy,xvarname)
        histoData.SetMarkerStyle(20)
        histoToy.SetLineColor(rt.kBlue)
        histoToy.SetLineWidth(2)

        c = rt.TCanvas()
        c.SetName('DataMC_%s_%s_ALLCOMPONENTS' % (xvarname,'_'.join(ranges)) )

        histoData.Draw("pe")
        histoToy.DrawCopy('histsame')
        histoToy.SetFillColor(rt.kBlue)
        histoToy.SetFillStyle(3018)
        histoToy.Draw('e2same')

        #histoData.Draw("pesame")

        #leg = rt.TLegend(0.6,0.6,0.9,0.9)
        #leg.SetFillColor(0)
        #leg.AddEntry(histoToyWln.GetName(),"W+jets","l")
        #leg.AddEntry(histoToyTTj.GetName(),"t#bar{t}","l")
        #leg.AddEntry(histoToy.GetName(),"Total","l")
        #leg.Draw()

        histToReturn = [histoToy, histoData, c]

        return histToReturn
