from RazorCombinedFit.Framework import Box
import RootTools
import ROOT as rt

class DalglishBox(Box.Box):
    
    def __init__(self, name, variables):
        super(DalglishBox,self).__init__(name, variables)

    def define(self, inputFile):
        
        #create the dataset
        data = RootTools.getDataSet(inputFile,'RMRTree')
        #import the dataset to the workspace
        self.importToWS(data)
        print 'Reduced dataset'
        #data.Print("V")

        # define the two components
        self.workspace.factory("RooDalglish::PDF1st(MR,Rsq,Xhat1,Yhat1,b1,MR01,R01,C1,S1,Sx1,Sy1,Xoff1,Yoff1)")
        self.workspace.factory("RooDalglish::PDF2nd(MR,Rsq,Xhat2,Yhat2,b2,MR02,R02,C2,S2,Sx2,Sy2,Xoff2,Yoff2)")
        #define the two yields
        self.workspace.factory("expr::N_1st('@0*(1-@1)',Ntot,f2)")
        self.workspace.factory("expr::N_2nd('@0*@1',Ntot,f2)")
        #associate the yields to the pdfs through extended PDFs
        self.workspace.factory("RooExtendPdf::ePDF1st(PDF1st, N_1st)")
        self.workspace.factory("RooExtendPdf::ePDF2nd(PDF2nd, N_2nd)")
        # build the total PDF
        model = rt.RooAddPdf("fitmodel", "fitmodel", rt.RooArgList(self.workspace.pdf("ePDF1st"),self.workspace.pdf("ePDF2nd")))        
        # import the model in the workspace.
        self.importToWS(model)
        #print the workspace
        self.workspace.Print()
        
    def plot(self, inputFile, store, box):
        super(DalglishBox,self).plot(inputFile, store, box)
        store.store(self.plotMR(inputFile), dir=box)
        store.store(self.plotRsq(inputFile), dir=box)
        store.store(self.plotRsqMR(inputFile), dir=box)
            
    def plotMR(self, inputFile):
        # project the data on R
        frameMR = self.workspace.var("MR").frame(self.workspace.var("MR").getMin(), 3000., 200)
        frameMR.SetName("MRplot")
        frameMR.SetTitle("MRplot")
        #        data = rt.RooDataSet(self.workspace.genobj("RMRTree"))
        #before I find a better way
        data = RootTools.getDataSet(inputFile,'RMRTree')
        data.plotOn(frameMR)
        # project the full PDF on the data
        self.workspace.pdf("fitmodel").plotOn(frameMR, rt.RooFit.LineColor(rt.kBlue))

        N1 = self.workspace.var("Ntot").getVal()*(1-self.workspace.var("f2").getVal())
        N2 = self.workspace.var("Ntot").getVal()*self.workspace.var("f2").getVal()

        # project the first component
        self.workspace.pdf("PDF1st").plotOn(frameMR, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1/(N1+N2)))
        # project the second component
        self.workspace.pdf("PDF2nd").plotOn(frameMR, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2/(N1+N2)))
        return frameMR

    def plotRsq(self, inputFile):
        # project the data on Rsq
        frameRsq = self.workspace.var("Rsq").frame(self.workspace.var("Rsq").getMin(), 1.5, 200)
        frameRsq.SetName("Rsqplot")
        frameRsq.SetTitle("Rsqplot")
        #before I find a better way
        data = RootTools.getDataSet(inputFile,'RMRTree')
        data.plotOn(frameRsq)

        N1 = self.workspace.var("Ntot").getVal()*(1-self.workspace.var("f2").getVal())
        N2 = self.workspace.var("Ntot").getVal()*self.workspace.var("f2").getVal()
        
        # project the full PDF
        self.workspace.pdf("fitmodel").plotOn(frameRsq, rt.RooFit.LineColor(rt.kBlue)) 
        # project the first component
        self.workspace.pdf("PDF1st").plotOn(frameRsq, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(8), rt.RooFit.Normalization(N1/(N1+N2)))
        # project the second component
        self.workspace.pdf("PDF2nd").plotOn(frameRsq, rt.RooFit.LineColor(rt.kBlue), rt.RooFit.LineStyle(9), rt.RooFit.Normalization(N2/(N1+N2)))

        return frameRsq

    def plotRsqMR(self, inputFile):
        #before I find a better way
        data = RootTools.getDataSet(inputFile,'RMRTree')
        toyData = self.workspace.pdf("fitmodel").generate(rt.RooArgSet(self.workspace.argSet("MR,Rsq")), 10*data.numEntries())

        # define 2D histograms
        histoData = rt.TH2D("histoData", "histoData",
                            100, self.workspace.var("MR").getMin(), 3500.,
                            100, self.workspace.var("Rsq").getMin(), 1.)
        histoToy = rt.TH2D("histoToy", "histoToy",
                            100, self.workspace.var("MR").getMin(), 3500.,
                            100, self.workspace.var("Rsq").getMin(), 1.)
        # project the data on the histograms
        data.tree().Project("histoData","Rsq:MR")
        toyData.tree().Project("histoToy","Rsq:MR")
        histoToy.Scale(histoData.Integral()/histoToy.Integral())
        histoData.Add(histoToy, -1)
        return histoData
