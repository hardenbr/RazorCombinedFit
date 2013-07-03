from RazorCombinedFit.Framework import Box
import RootTools
import ROOT as rt

class RazorBoxSignal(Box.Box):
    
    def __init__(self, name, variables):
        super(RazorBoxSignal,self).__init__(name, variables)
        
        self.zeros = {'Signal':[],'TTj':[],'Wln':['MuMu','EleEle','MuEle'],'Zll':['MuEle','Mu','Ele','Had'],'Znn':['Mu','Ele','MuMu','EleEle','MuEle','Had'],'QCD':['MuEle','MuMu','EleEle']}
        #self.cut = 'MR <= 750 && Rsq <= 0.2'

        self.cut = 'MR >= 0.0'

    def addTailPdf(self, flavour):
        
        label = '_%s' % flavour

        #define a flavour specific yield
        self.yieldToCrossSection(flavour)
        # define the two components
        self.workspace.factory("RooRazor2DTail::PDF1st"+label+"(MR,Rsq,MR01st"+label+",R01st"+label+",b1st"+label+")")
        self.workspace.factory("RooRazor2DTail::PDF2nd"+label+"(MR,Rsq,MR02nd"+label+",R02nd"+label+",b2nd"+label+")")
        #define the two yields
        self.workspace.factory("expr::N_1st"+label+"('@0*(1-@1)',Ntot"+label+",f2"+label+")")
        self.workspace.factory("expr::N_2nd"+label+"('@0*@1',Ntot"+label+",f2"+label+")")
        #associate the yields to the pdfs through extended PDFs
        self.workspace.factory("RooExtendPdf::ePDF1st"+label+"(PDF1st"+label+", N_1st"+label+")")
        self.workspace.factory("RooExtendPdf::ePDF2nd"+label+"(PDF2nd"+label+", N_2nd"+label+")")
        #float the efficiency with a penalty term if a sigma is provided


    def switchOff(self, species) :
        self.workspace.var("Epsilon_"+species).setVal(0.)
        self.workspace.var("Epsilon_"+species).setConstant(rt.kTRUE)
        self.workspace.var("f2_"+species).setConstant(rt.kTRUE)

    def define(self, inputFile):
        
        #create the dataset
        data = RootTools.getDataSet(inputFile,'RMRTree', self.cut)


        #create the signal PDF using RooHistPDF
        signal = RootTools.getDataSet('LM6_MR300.0_R0.3_'+self.name+'.root','RMRTree',self.cut)
        signalRDH = rt.RooDataHist("signal", "signal",rt.RooArgSet(self.workspace.var("MR"), self.workspace.var("Rsq")),signal)
        self.workspace.factory("expr::Nsig('%i*@0*@1/@2', Lumi, SigmaModel, Ngen)" % signal.numEntries())
        self.importToWS(signalRDH)
        signalModel = rt.RooHistPdf("signalModel", "signalModel", rt.RooArgSet(self.workspace.var("MR"),self.workspace.var("Rsq")), signalRDH)
        self.importToWS(signalModel)
        self.workspace.factory("RooExtendPdf::eBinPDF_Signal(signalModel, Nsig)")


        #import the dataset to the workspace
        self.importToWS(data)
        #self.workspace.Print('V')

        # add the different components:
        # - W+jets
        # - Zll+jets
        # - Znn+jets
        # - ttbar+jets
        self.addTailPdf("Wln")    
        self.addTailPdf("Zll")
        self.addTailPdf("Znn")
        self.addTailPdf("TTj")
        self.addTailPdf("QCD")

        # build the total PDF
        myPDFlist = rt.RooArgList(self.workspace.pdf("ePDF1st_Wln"),self.workspace.pdf("ePDF2nd_Wln"),
                                  self.workspace.pdf("ePDF1st_Zll"),self.workspace.pdf("ePDF2nd_Zll"),
                                                                    self.workspace.pdf("ePDF1st_Znn"),self.workspace.pdf("ePDF2nd_Znn"),
                                                                    self.workspace.pdf("ePDF1st_TTj"),self.workspace.pdf("ePDF2nd_TTj"))
        myPDFlist.add(self.workspace.pdf("ePDF1st_QCD"))
        myPDFlist.add(self.workspace.pdf("ePDF2nd_QCD"))    
        myPDFlist.add(self.workspace.pdf("eBinPDF_Signal"))    
        model = rt.RooAddPdf(self.fitmodel, self.fitmodel, myPDFlist)        
        
        # import the model in the workspace.
        self.importToWS(model)
        #print the workspace
        self.workspace.Print()

        ##### THIS IS A SIMPLIFIED FIT
        # fix all pdf parameters to the initial value
        self.fixPars("Zll")
        self.fixPars("Znn")
        self.fixPars("Wln")
        self.fixPars("TTj")
        self.fixPars("QCD")
        self.workspace.var("SigmaModel").setConstant(rt.kFALSE)
        self.workspace.var("Epsilon_TTj").setConstant(rt.kFALSE)
        self.workspace.var("Epsilon_Wln").setConstant(rt.kFALSE)
        self.workspace.var("Epsilon_QCD").setConstant(rt.kFALSE)
        
        
        
        
        #add penalty terms and float
        def float1stComponentWithPenalty(flavour):
            self.fixParsPenalty("MR01st_%s" % flavour)
            self.fixParsPenalty("R01st_%s" % flavour)
            self.fixParsPenalty("b1st_%s" % flavour)
        def float2ndComponentWithPenalty(flavour):
            self.fixParsPenalty("MR02nd_%s" % flavour)
            self.fixParsPenalty("R02nd_%s" % flavour)
            self.fixParsPenalty("b2nd_%s" % flavour)
        def floatFractionWithPenalty(flavour):
            self.fixParsPenalty("Epsilon_%s" % flavour, floatIfNoPenalty = True)
            self.fixPars("Epsilon_%s_s" % flavour)
            self.fixParsPenalty("f2_%s" % flavour)
            self.fixPars("f2_%s_s" % flavour)

        # switch off not-needed components (box by box)
        fixed = []
        for z in self.zeros:
            if self.name in self.zeros[z]:
                self.switchOff(z)
            #else:
                #if not z in fixed:
                    #float1stComponentWithPenalty(z)
                    #float2ndComponentWithPenalty(z)
                    #floatFractionWithPenalty(z)
                    #fixed.append(z)
        
    def plot(self, inputFile, store, box):
        super(RazorBoxSignal,self).plot(inputFile, store, box)
        store.store(self.plot1D(inputFile, "MR", 50, 300., 1500.), dir=box)
        store.store(self.plot1D(inputFile, "Rsq",50, 0.09, .5), dir=box)
        store.store(self.plot2D(inputFile, "MR", "Rsq"), dir=box)
            
    def plot1D(self, inputFile, varname, nbin=200, xmin=-99, xmax=-99):
        # set the integral precision
        rt.RooAbsReal.defaultIntegratorConfig().setEpsAbs(1e-10) ;
        rt.RooAbsReal.defaultIntegratorConfig().setEpsRel(1e-10) ;
        # get the max and min (if different than default)
        if xmax==xmin:
            xmin = self.workspace.var(varname).getMin()
            xmax = self.workspace.var(varname).getMax()
        # project the data on the variable
        frameMR = self.workspace.var(varname).frame(xmin, xmax, nbin)
        frameMR.SetName(varname+"plot")
        frameMR.SetTitle(varname+"plot")
        data = RootTools.getDataSet(inputFile,'RMRTree')
        data_cut = data.reduce(self.cut)
        data.plotOn(frameMR, rt.RooFit.LineColor(rt.kRed),rt.RooFit.MarkerColor(rt.kRed))
        data_cut.plotOn(frameMR)
        # project the full PDF on the data
        self.workspace.pdf(self.fitmodel).plotOn(frameMR, rt.RooFit.LineColor(rt.kBlue))

        # plot each individual component: Wln
        N1_Wln = self.workspace.function("Ntot_Wln").getVal()*(1-self.workspace.var("f2_Wln").getVal())
        N2_Wln = self.workspace.function("Ntot_Wln").getVal()*self.workspace.var("f2_Wln").getVal()
        # plot each individual component: Zll
        N1_Zll = self.workspace.function("Ntot_Zll").getVal()*(1-self.workspace.var("f2_Zll").getVal())
        N2_Zll = self.workspace.function("Ntot_Zll").getVal()*self.workspace.var("f2_Zll").getVal()
        # plot each individual component: Znn
        N1_Znn = self.workspace.function("Ntot_Znn").getVal()*(1-self.workspace.var("f2_Znn").getVal())
        N2_Znn = self.workspace.function("Ntot_Znn").getVal()*self.workspace.var("f2_Znn").getVal()
        # plot each individual component: TTj
        N1_TTj = self.workspace.function("Ntot_TTj").getVal()*(1-self.workspace.var("f2_TTj").getVal())
        N2_TTj = self.workspace.function("Ntot_TTj").getVal()*self.workspace.var("f2_TTj").getVal()
        # plot each individual component: QCD
        N1_QCD = self.workspace.function("Ntot_QCD").getVal()*(1-self.workspace.var("f2_QCD").getVal())
        N2_QCD = self.workspace.function("Ntot_QCD").getVal()*self.workspace.var("f2_QCD").getVal()        
        # plot each individual component: Signal
        Nsig = self.workspace.function("Nsig").getVal()

        Ntot = N1_Wln+N2_Wln+N1_Zll+N2_Zll+N1_Znn+N2_Znn+N1_TTj+N2_TTj+N1_QCD+N2_QCD+Nsig

        if N1_Wln+N2_Wln >0:
            # project the first component: Wln
            self.workspace.pdf("PDF1st_Wln").plotOn(frameMR, rt.RooFit.LineColor(rt.kRed), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1_Wln/Ntot))
            # project the second component: Wln
            self.workspace.pdf("PDF2nd_Wln").plotOn(frameMR, rt.RooFit.LineColor(rt.kRed), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2_Wln/Ntot))
        if N1_Zll+N2_Zll >0:
            # project the first component: Zll
            self.workspace.pdf("PDF1st_Zll").plotOn(frameMR, rt.RooFit.LineColor(rt.kMagenta), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1_Zll/Ntot))
            # project the second component: Zll
            self.workspace.pdf("PDF2nd_Zll").plotOn(frameMR, rt.RooFit.LineColor(rt.kMagenta), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2_Zll/Ntot))
        if N1_Znn+N2_Znn >0:
            # project the first component: Znn
            self.workspace.pdf("PDF1st_Znn").plotOn(frameMR, rt.RooFit.LineColor(rt.kGreen), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1_Znn/Ntot))
            # project the second component: Znn
            self.workspace.pdf("PDF2nd_Znn").plotOn(frameMR, rt.RooFit.LineColor(rt.kGreen), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2_Znn/Ntot))
        if N1_TTj+N2_TTj >0:
            # project the first component: TTj
            self.workspace.pdf("PDF1st_TTj").plotOn(frameMR, rt.RooFit.LineColor(rt.kOrange), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1_TTj/Ntot))
            # project the second component: TTj
            self.workspace.pdf("PDF2nd_TTj").plotOn(frameMR, rt.RooFit.LineColor(rt.kOrange), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2_TTj/Ntot))
        if N1_QCD+N2_QCD >0:
            # project the first component: QCD
            self.workspace.pdf("PDF1st_QCD").plotOn(frameMR, rt.RooFit.LineColor(rt.kViolet), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1_QCD/Ntot))
            # project the second component: QCD
            self.workspace.pdf("PDF2nd_QCD").plotOn(frameMR, rt.RooFit.LineColor(rt.kViolet), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2_QCD/Ntot))            
        if Nsig>0:
            # project the first component: Signal
            self.workspace.pdf("signalModel").plotOn(frameMR, rt.RooFit.LineColor(rt.kBlack), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(Nsig/Ntot)) 
            
        #leg = rt.TLegend("leg", "leg", 0.6, 0.6, 0.9, 0.9)
        #leg.AddEntry("PDF1st_Wln", "W+jets 1st")
        #leg.AddEntry("PDF2nd_Wln", "W+jets 2nd")
        #leg.AddEntry("PDF1st_Zll", "Z(ll)+jets 1st")
        #leg.AddEntry("PDF2nd_Zll", "Z(ll)+jets 2nd")
        #leg.AddEntry("PDF1st_Znn", "Z(#nu#nu)+jets 1st")
        #leg.AddEntry("PDF2nd_Znn", "Z(#nu#nu)+jets 2nd")
        #leg.AddEntry("PDF1st_TTj", "t#bar{t}+jets 1st")
        #leg.AddEntry("PDF2nd_TTj", "t#bar{t}+jets 2nd")
        
        return frameMR

    def plot2D(self, inputFile, xvarname, yvarname):
        #before I find a better way
        data = RootTools.getDataSet(inputFile,'RMRTree')
        toyData = self.workspace.pdf("fitmodel").generate(rt.RooArgSet(self.workspace.argSet(xvarname+","+yvarname)), 20*data.numEntries())
        toyData = toyData.reduce(self.getVarRangeCut())

        # define 2D histograms
        histoData = rt.TH2D("histoData", "histoData",
                            100, self.workspace.var(xvarname).getMin(), self.workspace.var(xvarname).getMax(), 
                            100, self.workspace.var(yvarname).getMin(), self.workspace.var(yvarname).getMax())
        histoToy = rt.TH2D("histoToy", "histoToy",
                            100, self.workspace.var(xvarname).getMin(), self.workspace.var(xvarname).getMax(), 
                            100, self.workspace.var(yvarname).getMin(), self.workspace.var(yvarname).getMax())
        # project the data on the histograms
        data.tree().Project("histoData",yvarname+":"+xvarname)
        toyData.tree().Project("histoToy",yvarname+":"+xvarname)
        histoToy.Scale(histoData.Integral()/histoToy.Integral())
        histoData.Add(histoToy, -1)
        return histoData
