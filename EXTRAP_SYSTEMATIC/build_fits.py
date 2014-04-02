import os
import ROOT as rt
import numpy, itertools
from  optparse  import OptionParser

rt.gROOT.SetBatch(True)

mr_max = 7
#from "working" fit razor_output_samp1 and 0.root
#(m0_0, m0_min, m0_max) = (-.4,-1,0)
#(n_0, n_min, n_max) = (3, 1.2, 100)
#(b_0, b_min, b_max) = (7, 1, 100)

#stable samp0
#(m0_0, m0_min, m0_max) = (.1,-.3,.3)
#(n_0, n_min, n_max) = (15, 1.3, 1000)
#(b_0, b_min, b_max) = (7, 1, 100)

#stable samp1
(m0_0, m0_min, m0_max) = (0,-1,1)
(n_0, n_min, n_max) = (6, 4, 100)
(b_0, b_min, b_max) = (8, 1, 10)

#test
#(m0_0, m0_min, m0_max) = (.1,-.3,.3)
#(n_0, n_min, n_max) = (15, 1.3, 1000)
#(b_0, b_min, b_max) = (7, 1, 100)

def make_cfg(filename, mr_min, rsq_min, rsq_max, n_tot):
    outfile = open(filename,"w")
    outfile.write("[Had]\n")
    variables_tuple = (mr_min, mr_min, mr_max, rsq_min, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables = ['MR[%f, %f, %f]', 'Rsq[%f, %f, %f]', 'nBtag[%i,%i,%i]'] \n" % variables_tuple)
    variables_range_tuple = (mr_min, mr_max, rsq_min, rsq_max, 0,0,4)
    outfile.write("variables_range = ['MR_FULL[%f, %f]','Rsq_FULL[%f, %f]','nBtag_FULL[%i,%i,%i]']\n" % variables_range_tuple)
    outfile.write("pdf_QCD = ['MR0_QCD[%f,%f,%f]','R0_QCD[1, 0, 1.5]','R1_QCD[.0,0,1]','b_QCD[%f,%f,%f]','RI_QCD[.000,0,.003]','n_QCD[%f, %f, %f]']\n" % (m0_0,m0_min,m0_max, b_0, b_min, b_max, n_0, n_min, n_max))
    outfile.write("others_QCD = ['Lumi[4980]','Ntot_QCD[%f,%f,%f]']\n" % (n_tot/2., 0, n_tot*2))
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

parser.add_option("-d", "--roodata", dest="roodata",
                  help="folder containing roodatasets. If this is specified the datasets will not be build, only configs.",
                  action="store",default=None, type="string")

parser.add_option("--iSamp", dest="iSamp",
                  help="sample to use",
                  action="store",default=1, type="int")



(options, args) = parser.parse_args()

file = rt.TFile(options.filename,"READ")
tree = file.Get("HggOutput")

sig_file = rt.TFile(options.sig,"READ")
sig_tree = sig_file.Get("HggOutput")

#set the cut values to probe
mr_cuts = [.1 * x for x in xrange(6,7)]
low_rsq_cut = [.005 * x for x in xrange(2,3)]
low_rsq_width = [.005 * x for x in xrange(2,3)]

print "MR SCAN", mr_cuts
print "LOW RSQ SCAN", low_rsq_cut
print "LOW BAND WIDTH SCAN", low_rsq_width

mr_ar= numpy.array(mr_cuts)
low_rsq_ar =  numpy.array(low_rsq_cut)
width_rsq_ar = numpy.array(low_rsq_width)

#take the direct product of the possible cuts
grid = list(itertools.product(mr_ar, low_rsq_ar, width_rsq_ar))
config_names  = []

grid = grid[:86]

#build a config for each of the points
results = []
for point in grid:
    print "Scanning point....", grid.index(point) + 1 , " of ", len(grid)

    #parse the grid point
    (mr_min, rsq_min, rsq_width) = point
    rsq_max = rsq_min + rsq_width

    if options.roodata != None:
        #get a list of all the fit files in the dataset directory
        path = options.roodata
        files = os.listdir(path)
        #get a list of all of the config files
        for file in files:
            if "mr%.1f_rsq%.4f_rsqw%.4f_rooDat_Had.root" % (mr_min, rsq_min, rsq_width) in file:
                path = options.cfg_loc
                name = path + "/" + "config_mr%.1f_rsq%.4f_rsqw%.4f.cfg" % point
                make_cfg(name, mr_min, rsq_min, rsq_max, 10000)                
                config_names.append(name)
        continue

    n_events_low = tree.GetEntries("PhotonPFCiC.sieie[1] > .0001 && PhotonPFCiC.sieie[0] > .0001 && PFMR > %f && PFR^2 > %f && PFR^2 < %f && iSamp == %i" % (mr_min * 1000., rsq_min, rsq_max, options.iSamp))
    n_events_high = tree.GetEntries("PhotonPFCiC.sieie[1] > .0001 && PhotonPFCiC.sieie[0] > .0001 && PFMR > %f && PFR^2 > %f  && iSamp == %i" % (mr_min * 1000., rsq_max, options.iSamp))
    n_events_low_sig = sig_tree.GetEntries("PFMR > %f && PFR^2 > %f && PFR^2 < %f && iSamp == 0" % (mr_min * 1000., rsq_min, rsq_max))
    n_events_high_sig = sig_tree.GetEntries("PFMR > %f && PFR^2 > %f  && iSamp == 0" % (mr_min * 1000., rsq_max))
    n_events_sig_total = sig_tree.GetEntries("iSamp==0 && PFMR>0")

    perc_contam = float(n_events_low_sig) / 10000.
    perc_efficient = float(n_events_high_sig) / float(n_events_sig_total)
    
    print "Number of events low...", n_events_low
    print "Number of events high...", n_events_high
    print "\n"
    print "Number of signal contaminated...", n_events_low_sig
    print "Percent of signal contaminated...", perc_contam
    print "Percent of signal in high rsq bin...", perc_efficient
    #build the config name
    path = options.cfg_loc
    name = path + "/" + "config_mr%.1f_rsq%.4f_rsqw%.4f.cfg" % point
    print name
    
    #build the config
    if n_events_high > 500 and n_events_low < 10000 and n_events_low > 1000 and perc_contam < .05 and perc_efficient > .8:
        print "\t\t@@@POINT PASSES GENERAL REQUIREMENTS@@@"
        make_cfg(name, mr_min, rsq_min, rsq_max, n_events_low)
        config_names.append(name)
        result = (mr_min, rsq_min, rsq_width, n_events_high, n_events_low, n_events_sig_total, n_events_low_sig, n_events_high_sig, perc_contam, perc_efficient)
        results.append(result)


#use each config to convert the dataset to a roodataset
for name in config_names:
    if options.roodata != None: break #if roodata has already be generated dont convert the datasets    
    outdata = name[:-4] + "_rooDat_Had.root"    
    cmd = "python Josh2Dataset.py -c %s -o %s --samp %i %s" % (name, outdata, options.iSamp, options.filename)
    print "Converting dataset....", config_names.index(name) + 1, " of ", len(config_names)
    print cmd
    os.system(cmd)

fit_results = []
warnings = []
#use each roodataset to complete the fit
for config in config_names:
    outdata = config[:-4] + "_rooDat_Had.root"        
    cmd = "python ../scripts/runAnalysis.py --fitmode 2D -a SingleBoxFit --gamma -c %s %s" % (config, outdata)    
    fit_output  = config[:-4] + "_Fit.root"        

    print cmd
    os.system(cmd)
    
    fit_res = rt.TFile("razor_output.root").Get("Had/independentFR")
    fr_list = fit_res.floatParsFinal()

    b_val = fr_list[2].getVal()
    n_val = fr_list[3].getVal()
    norm_val = fr_list[1].getVal()
    m0_val = fr_list[0].getVal()

    b_error = fr_list[2].getError()
    n_error = fr_list[3].getError()
    norm_error = fr_list[1].getError()
    m0_error = fr_list[0].getError()

    bdiff = min(b_val - b_min, b_max - b_val ) / b_val
    ndiff = min(n_val - n_min, n_max - n_val ) / n_val
    m0diff = 0
    
    if m0_val != 0:
        m0diff = min(m0_val - m0_min, m0_max - m0_val) / m0_val
    else:
        m0diff = min(m0_val - m0_min, m0_max - m0_val) 

    if abs(bdiff) < .05:
        warnings.append("\t\t\t WARNING B PARAM NEAR LIMIT "  + config)
    if abs(ndiff) < .05:
        warnings.append( "\t\t\t WARNING N PARAM NEAR LIMIT " + config)
    if abs(m0diff) < .05:
        pass #warnings.append( "\t\t\t WARNING m0 PRAM  NEAR LIMIT " + config)

    fit_result = (norm_val, m0_val, m0_error, n_val, n_error,  b_val, b_error,  config)

    fit_results.append(fit_result)

    os.system("mv razor_output.root %s" % fit_output)

#print the regional information
print  "mr_min\trsq_min\trsq_w\tn_high\tn_low\tns_tot\tns_low\tns_hi\tps_cont\tps_eff "
for result in results:
    for ii in result:
        if ii > 10:
            print "%.1f \t" % ii,
        else:
            print "%.4f \t" % ii,
    print "\n"

print "\n"
#print the parameters of the fit
print "norm\t\tm0\tm0_e\tn\tn_e\tb\tb_e"
for param in fit_results:
    print param[-1]
    for val in param[:-1]: print "%2.3f\t" % val,
    print "\n"

#print the warnings
for warning in warnings: print warning
    
