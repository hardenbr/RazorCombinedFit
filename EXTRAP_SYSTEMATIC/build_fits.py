import ROOT as rt
import numpy, itertools
from  optparse  import OptionParser

rt.gROOT.SetBatch(True)

mr_max = 7

def make_cfg(filename, mr_min, rsq_min, rsq_max, n_tot):
    outfile = open(filename,"w")
    outfile.write("[Had]\n")
    variables_tuple = (mr_min, mr_min, mr_max, rsq_min, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables = ['MR[%f, %f, %f]', 'Rsq[%f, %f, %f], 'nBtag[%i,%i,%i] \n" % variables_tuple)
    variables_range_tuple = (mr_min, mr_max, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables_range = ['MR_FULL[%f, %f]','Rsq_FULL[%f, %f]','nBtag_FULL[%i,%i,%i]']\n" % variables_range_tuple)
    outfile.write("pdf_QCD = ['MR0_QCD[-.4,-1,0]','R0_QCD[1, 0, 1.5]','R1_QCD[.0,0,1]','b_QCD[7,1,100]','RI_QCD[.000,0,.003]','n_QCD[3., 1.2, 200]']\n")
    outfile.write("others_QCD = ['Lumi[4980]','Ntot_QCD[%i,%i,%i]']\n" % (n_tot/2., 0, n_tot*1.2))
    outfile.close()
    
cmd_build_dataset = "python scripts/Josh2Dataset.py -c gg_razor_fit_config/SingleBoxFit_Diphoton /afs/cern.ch/work/h/hardenbr/2013/RAZOR_DIPHOTON/DATA/2012ABCD" 

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                                    help="HggOutput.root file name to analyze",
                                    action="store",type="string")

parser.add_option("-c", "--cfg", dest="cfg_loc",
                                    help="location to score the configs",
                                    action="store",type="string")
(options, args) = parser.parse_args()

file = rt.TFile(options.filename,"READ")
tree = file.Get("HggOutput")

#set the cut values to probe
mr_cuts = [100 * x for x in xrange(2,8)]
low_rsq_cut = [.005 * x for x in xrange(3)]
low_rsq_width = [.005 * x for x in xrange(3)]

mr_ar= numpy.array(mr_cuts)
low_rsq_ar =  numpy.array(low_rsq_cut)
width_rsq_ar = numpy.array(low_rsq_width)

#take the direct product of the possible cuts
grid = list(itertools.product(mr_ar, low_rsq_ar, width_rsq_ar))
config_names  = []

#build a config for each of the points
for point in grid:
    #parse the grid point
    (mr_min, rsq_min, rsq_width) = point
    rsq_max = rsq_min + rsq_width
    n_events = tree.GetEntries("PFMR > %f && PFR^2 > %f && PFR^2 < %f && iSamp == 1" % (mr_min*1000, rsq_min*1000, rsq_max*1000))
    #build the config name
    path = options.cfg_loc
    name = path + "/" + "config_mr%.1f_rsq%.3f_rsqw%.3f.cfg" % point
    config_names.append(name)
    
    #build the config
    make_cfg(name, mr_min, rsq_min, rsq_max, n_events)


#use each config to convert the dataset to a roodataset


#use each roodataset to complete the fit
