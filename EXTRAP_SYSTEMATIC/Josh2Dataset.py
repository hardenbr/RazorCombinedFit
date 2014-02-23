from optparse import OptionParser

import ROOT as rt
import RootTools
from RazorCombinedFit.Framework import Config
import math

def convertTree2Dataset(tree, outputFile, outputBox, config, box, min, max, scale_factor):
    """This defines the format of the RooDataSet"""
    
    workspace = rt.RooWorkspace(box)
    variables = config.getVariables(box,"variables")
    for v in variables:
        workspace.factory(v)

    args = workspace.allVars()
    data = rt.RooDataSet('RMRTree','Selected R and MR',args)
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin() * scale_factor
    print mRmin
    mRmax = args['MR'].getMax() * scale_factor
    #we cut away events outside our R window
    rMin = math.sqrt(args['Rsq'].getMin())
    rMax = math.sqrt(args['Rsq'].getMax())
    print "R MIN: " + str(rMin)

    #rcuts = config.getRCuts(box)
    #rcuts.sort()
    #rMin = rcuts[0]

    #iterate over selected entries in the input tree    
    if options.scale_factor == 1:
        tree.Draw('>>elist','PFMR  >= %f && PFMR <= %f && PFR >= %f && PFR <= %f' % (mRmin,mRmax,rMin,rMax),'entrylist')    
    else:
        tree.Draw('>>elist','PhotonPFCiC.sieie[1] > .0001 && PhotonPFCiC.sieie[0] > .0001 && PFMR  >= %f && PFMR <= %f && PFR >= %f && PFR <= %f &&iSamp==%i' % (mRmin,mRmax,rMin,rMax,options.samp),'entrylist')

    elist = rt.gDirectory.Get('elist')
    
    entry = -1;
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)

        #set the RooArgSet and save
        a = rt.RooArgSet(args)
        a.setRealValue('MR',tree.PFMR/scale_factor)
        a.setRealValue('R', (rMin + rMax) / 2)
        a.setRealValue('nBtag', 0.)
        a.setRealValue('Rsq',((rMin + rMax) / 2)*((rMin + rMax) / 2))
        data.add(a)
    numEntries = data.numEntries()
    if min < 0: min = 0
    if max < 0: max = numEntries
    
    rdata = data.reduce(rt.RooFit.EventRange(min,max))
    output = None
    
    if options.outfile != None:
        output = rt.TFile(options.outfile,"RECREATE")
        print "Writing", options.outfile
    else:
        output = rt.TFile.Open(outputFile+"_MR"+str(mRmin)+"_Rsq"+str(rMin*rMin)+"_"+str(rMax*rMax)+"_"+"Samp_"+str(options.samp)+"_"+outputBox,'RECREATE')
        print 'Writing',outputFile+"_MR"+str(mRmin)+"_Rsq"+str(rMin*rMin)+"_"+str(rMax*rMax)+"_"+outputBox

    rdata.Write()
    output.Close()
    
    return rdata.numEntries()
    
if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default=None,
                  help="Name of the config file to use")
    parser.add_option('--max',dest="max",type="int",default=-1,
                  help="The last event to take from the input Dataset")
    parser.add_option('--min',dest="min",type="int",default=0,
                  help="The first event to take from the input Dataset")
    parser.add_option('--scale',dest="scale_factor",type="int",default=1000,
                  help="Scale M_R by this amount i.e. GeV-->TeV use 1000")
    parser.add_option('--sample',dest="samp",type="int",default=1,
                      help="samp=1(fake) samp=0(tight)")
    parser.add_option('-o','--outfile',dest="outfile",type="string",default=None,
                      help="Name the roodataset file will take")
    (options,args) = parser.parse_args()


    if options.config is None:
        import inspect, os
        topDir = os.path.abspath(os.path.dirname(inspect.getsourcefile(convertTree2Dataset)))
        options.config = os.path.join(topDir,'boxConfig.cfg')    
    cfg = Config.Config(options.config)
    
    print 'Input files are %s' % ', '.join(args)
    for f in args:
        if f.lower().endswith('.root'):
            input = rt.TFile(f)
            decorator = f[:-5]
            
            #dump the trees for the different datasets
            convertTree2Dataset(input.Get('HggOutput'), decorator, 'Had.root', cfg,'Had',options.min,options.max,options.scale_factor)
            #convertTree2Dataset(input.Get('outTreeGG'), decorator, 'GaGa.root', cfg,'GaGa',options.min,options.max)
            #convertTree2Dataset(input.Get('outTreeGE'), decorator, 'Ele.root', cfg,'Ele',options.min,options.max)
            #convertTree2Dataset(input.Get('outTreeGM'), decorator, 'Mu.root', cfg,'Mu',options.min,options.max)
            #convertTree2Dataset(input.Get('outTreeGMM'), decorator, 'MuMu.root', cfg,'MuMu',options.min,options.max)
            #convertTree2Dataset(input.Get('outTreeMuEle'), decorator, 'MuEle.root', cfg,'MuEle',options.min,options.max)
            #convertTree2Dataset(input.Get('outTreeGEE'), decorator, 'EleEle.root', cfg,'EleEle',options.min,options.max)
            
        else:
            "File '%s' of unknown type. Looking for .root files only" % f
