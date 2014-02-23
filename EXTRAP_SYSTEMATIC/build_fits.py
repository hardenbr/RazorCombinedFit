import os
import ROOT as rt
import numpy, itertools
from  optparse  import OptionParser

rt.gROOT.SetBatch(True)

mr_max = 7

def make_cfg(filename, mr_min, rsq_min, rsq_max, n_tot):
    outfile = open(filename,"w")
    outfile.write("[Had]\n")
    variables_tuple = (mr_min, mr_min, mr_max, rsq_min, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables = ['MR[%f, %f, %f]', 'Rsq[%f, %f, %f]', 'nBtag[%i,%i,%i]'] \n" % variables_tuple)
    variables_range_tuple = (mr_min, mr_max, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables_range = ['MR_FULL[%f, %f]','Rsq_FULL[%f, %f]','nBtag_FULL[%i,%i,%i]']\n" % variables_range_tuple)
    outfile.write("pdf_QCD = ['MR0_QCD[-.4,-1,0]','R0_QCD[1, 0, 1.5]','R1_QCD[.0,0,1]','b_QCD[7,1,100]','RI_QCD[.000,0,.003]','n_QCD[3., 1.2, 70]']\n")
    outfile.write("others_QCD = ['Lumi[4980]','Ntot_QCD[%f,%f,%f]']\n" % (n_tot/2., 0, n_tot*1.2))
    outfile.close()

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                                    help="HggOutput.root file name to analyze",
                                    action="store",type="string")

parser.add_option("-c", "--cfg", dest="cfg_loc",
                                    help="location to score the configs",
                                    action="store",type="string")

parser.add_option("-s", "--sig", dest="sig",
                                    help="signal file to check for contamination",
                                    action="store",type="string")


(options, args) = parser.parse_args()

file = rt.TFile(options.filename,"READ")
tree = file.Get("HggOutput")

sig_file = rt.TFile(options.sig,"READ")
sig_tree = sig_file.Get("HggOutput")

#set the cut values to probe
mr_cuts = [.2 * x for x in xrange(1,8)]
low_rsq_cut = [.005 * x for x in xrange(1,10)]
low_rsq_width = [.005 * x for x in xrange(1,4)]

print "MR SCAN", mr_cuts
print "LOW RSQ SCAN", low_rsq_cut
print "LOW BAND WIDTH SCAN", low_rsq_width

mr_ar= numpy.array(mr_cuts)
low_rsq_ar =  numpy.array(low_rsq_cut)
width_rsq_ar = numpy.array(low_rsq_width)

#take the direct product of the possible cuts
grid = list(itertools.product(mr_ar, low_rsq_ar, width_rsq_ar))
config_names  = []

grid = grid

#build a config for each of the points
for point in grid:
    print "Scanning point....", grid.index(point) + 1 , " of ", len(grid)

    #parse the grid point
    (mr_min, rsq_min, rsq_width) = point
    rsq_max = rsq_min + rsq_width

    n_events_low = tree.GetEntries("PFMR > %f && PFR^2 > %f && PFR^2 < %f && iSamp == 1" % (mr_min * 1000., rsq_min, rsq_max))
    n_events_high = tree.GetEntries("PFMR > %f && PFR^2 > %f  && iSamp == 1" % (mr_min * 1000., rsq_max))
    n_events_low_sig = sig_tree.GetEntries("PFMR > %f && PFR^2 > %f && PFR^2 < %f && iSamp == 0" % (mr_min * 1000., rsq_min, rsq_max))
    perc_contam = float(n_events_low_sig) / 10000.

    print "Number of events low...", n_events_low
    print "Number of events high...", n_events_high
    print "Number of signal contaminated...", n_events_low_sig
    print "Percent of signal contaminated...", perc_contam
    
    #build the config name
    path = options.cfg_loc
    name = path + "/" + "config_mr%.1f_rsq%.4f_rsqw%.4f.cfg" % point
    print name
    
    #build the config
    if n_events_high > 500 and n_events_low < 10000 and n_events_low > 1000 and perc_contam < .05:
        make_cfg(name, mr_min, rsq_min, rsq_max, n_events_low)
        config_names.append(name)

#get a list of all the fit files in the dataset directory
#path_elem = ((options.filename).split("/"))[:-1]
#path = ""
#for elem in path_elem: path+= elem + "/"
#files = os.listdir(path)
#converted_files = []
#get a list of all of the converted files
#for file in files:
#    if "MR" in file and "Had" in file:
#        print "Appending....", path+"/"+file
#        converted_files.append(path+"/"+file)


#use each config to convert the dataset to a roodataset
for name in config_names:
    outdata = name[:-4] + "_rooDat_Had.root"    
    cmd = "python Josh2Dataset.py -c %s -o %s  %s" % (name, outdata, options.filename)
    print "Converting dataset....", config_names.index(name) + 1, " of ", len(config_names)
    print cmd
    os.system(cmd)


#use each roodataset to complete the fit
for name in config_names:
    config = name
    outdata = name[:-4] + "_rooDat_Had.root"        
    cmd = "python ../scripts/runAnalysis.py --fitmode 2D -a SingleBoxFit --gamma -c %s %s" % (config, outdata)    
    fit_output  = name[:-4] + "_Fit.root"        

    print cmd

    os.system(cmd)
    os.system("mv razor_output.root %s" % fit_output)
