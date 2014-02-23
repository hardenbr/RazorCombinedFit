from  optparse  import OptionParser
import ROOT as rt
import RootTools
import array, os
import random
import rootlogon
rt.gROOT.SetBatch(True);
DO_DEBUG = False

rootlogon.style()

def smallest_68(hist):
    total_events = hist.Integral()
    n_bins = hist.GetNbinsX()
    
    sum = 0.0
    for bin in range(n_bins):
        count = hist.GetBinContent(bin+1)
        sum += float(count)
        
        if (sum / total_events) > .68:
            return hist.GetBinCenter(bin+1)

    return 0

#exponentially increasing bin sizes
def makebins(start_,end_,inc_,inc_inc_):
    bin = start_
    inc = inc_
    list = []
    while bin < end_:
        list.append(bin)
        bin+=inc
        inc*=(1+inc_inc_)
    return list

def Gamma(a,x):
    val = rt.TMath.Gamma(a)*rt.Math.inc_gamma_c(a,x)
    if DO_DEBUG: print "gamma", val
    return val

def Gfun(m,r,b,n):
    val =  n/pow(b*n,n)*(Gamma(n,b*n*pow(m*r,1/n)))
    if DO_DEBUG: print "gfun", val
    return val

def integral(r, m1, m2, b, M0, n):
    return (Gfun(m1-M0,r,b,n) - Gfun(m2-M0,r,b,n))

def integral_fit(m1, m2, b, M0 ,n):
    return (Gfun(m1-M0,1,b,n) - Gfun(m2-M0,1,b,n))

def meanstdv(x):
    from math import sqrt
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
    std = sqrt(std / float(n-1))
    return mean, std

def get_sigma(low, high):
    sigmas_low = []
    sigmas_high = []

    for ii in range(len(low)):
        (mean_low,sigma_low) = meanstdv(low[ii])
        (mean_high,sigma_high) = meanstdv(high[ii])

        sigmas_low.append(sigma_low)
        sigmas_high.append(sigma_high)
        
    return (sigmas_low, sigmas_high)

def get_xsec(xsec_file, msq, mgl):
    file = open(xsec_file,"r")
    lines = file.readlines()
    lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

    #scan the lines
    for line in lines:
        if msq in line and mgl in line:
            #check the ordering is correct
            split = line.split()
            idx1 =  split.index(msq)
            idx2 =  split.index(mgl)
            if idx1 < idx2:
                #parse the values
                lo_xsec = float(split[6])
                lo_xsec_plus = float(split[8])
                lo_xsec_minus = float(split[10])

                return (lo_xsec,lo_xsec_plus,lo_xsec_minus)
            
    print "\n\nERROR:CROSS SECTION WAS NOT FOUND\n\n"
    return (-1,-1,-1)

#generates a randomized signal sample based on the cross section file, signal
#strength, and signal file
def build_mixed_model(signal_file, xsec_file, mu, mr_min, rsq_min, outfile_name):

    mr_min = mr_min*1000
    
    signal_name = signal_file.split("/")[-1]
    msq = signal_name.split("_")[1]
    mgl = signal_name.split("_")[2]    

    (xsec,xsec_plus,xsec_minus) = get_xsec(xsec_file, msq, mgl)

    signal = rt.TFile(signal_file)
    signal_tree = signal.Get("HggOutput")
    
    eff = signal_tree.GetEntries("PFMR > %f & PFR^2 > %f" % (mr_min, rsq_min)) / 10000.

    #calculate the number of signal events expected
    n_events = int(eff * (xsec * 1000) * 20 * mu)
    n_events_plus = int(eff * ((xsec+xsec_plus) * 1000) * 20 * mu)    
    n_events_minus = int(eff * ((xsec-xsec_minus) * 1000) * 20 * mu)

    #make the iteration list for signal
    signal_tree.Draw(">>iterlist","","entrylist")
    itlist = rt.gDirectory.Get("iterlist")

    signal_mr_rsq = []
    #make a list of all the values for that file
    for event in range(signal_tree.GetEntries()):
        if event % 5000 == 0: print "scanning signal point....", event
        entry = itlist.Next()        
        signal_tree.GetEntry(entry)

        if signal_tree.PFMR > 0:
            signal_mr_rsq.append((signal_tree.PFMR,signal_tree.PFR))

    #randomize the entries
    random.shuffle(signal_mr_rsq)

    #output variables
    rt.gROOT.ProcessLine("struct myStruct{\
    Float_t PFMR;\
    Float_t PFR;\
    };")

    outfile = rt.TFile(outfile_name,"RECREATE")
    s = rt.myStruct()

    #generate the output tree
    outTree = rt.TTree("HggOutput","HggOutput")
    outTree.Branch("PFMR",rt.AddressOf(s,"PFMR"),"PFMR/F")
    outTree.Branch("PFR",rt.AddressOf(s,"PFR"),"PFR/F")

    print "Filling result signal..."

    #fill the tree
    for ii in signal_mr_rsq[:n_events]:
        print ii[0],ii[1]*ii[1]
        s.PFMR = float(ii[0])/ 1000.
        s.PFR = float(ii[1])
        outTree.Fill()

    outTree.Write()
    outfile.Close()
    
    return (eff,n_events, n_events_plus, n_events_minus)


def get_bin_errors(fr, norm):

    high_list = []
    low_list = []
    hists_low = []
    hists_high = []

    for ii in range(len(bins)-1):
        high_list.append([])
        low_list.append([])
        name_high = "high_bin%i" % ii
        name_low = "low_bin%i" % ii

        if ii < 4:
            hist_high = rt.TH1D(name_high,name_high, 100000,0,100000)
            hist_low = rt.TH1D(name_low,name_low, 10000,0,100000)
        else:
            hist_high = rt.TH1D(name_high,name_high, 50,0,50)
            hist_low = rt.TH1D(name_low,name_low, 1000000,0,100000)
            
        hists_low.append(hist_low)
        hists_high.append(hist_high)
        
    for ii in range(options.n_toys):

        if ii % 1000 == 0: print "Running Toy....%i" % ii

        #sample from the covariance matrix
        list = fr.randomizePars() 

        #extract parameters from the fit
        M0_toy = list[0].getVal()
        Ntot_toy = list[1].getVal()
        b_toy = list[2].getVal()
        n_toy = list[3].getVal()

        #cant be evaluated
#        if n_toy > 80: continue
        
        #derive the high rsq region parameters
        b_high_toy = b_toy #* pow(r1*r1,1/n_toy)
        n_high_toy = n_toy
        M0_high_toy = M0_toy

        if DO_DEBUG: print "b", b_high_toy, "n",n_high_toy

        #find the fit normalization
        if DO_DEBUG: print "\ndoing low integral.."
        int_low = integral_fit(mr_min_low, mr_cutoff, b_toy, M0_toy, n_toy)
        if DO_DEBUG: print "\ndoing high integral.."
        int_high = integral_fit(mr_min, mr_cutoff, b_high_toy, M0_high_toy, n_high_toy)

        if int_low ==0 or int_high == 0:continue


        #do the integral for each bin        
        for jj in range(len(bins) - 1):

            m1 = bins[jj]
            m2 = bins[jj+1]

            #do the integral for the bin
            high = integral_fit(m1, m2, b_high_toy, M0_high_toy, n_high_toy) * norm / int_high 
            low = integral_fit(m1, m2, b_toy, M0_toy, n_toy) * Ntot_toy / int_low
            
            pois_high = trandom.PoissonD(high)
            pois_low = trandom.PoissonD(low)

            hists_low[jj].Fill(pois_low)
            hists_high[jj].Fill(pois_high)

            low_list[jj].append(pois_low)            
            high_list[jj].append(pois_high)

    return (low_list,high_list,hists_low,hists_high)

def assign_values_from_fit(hist_low_rsq, hist_high_rsq, toy_hists_low, toy_hists_high ):
    for ii in range(len(bins) - 1):
        m1 = bins[ii]
        m2 = bins[ii+1]
        
        #normalization integrals
        if DO_DEBUG: print "doing low integral.."    
        int_low = integral_fit(mr_min, mr_cutoff, b_low, M0, n)
        if DO_DEBUG: print "doing high integral.."
        int_high = integral_fit(mr_min, mr_cutoff, b_high, M0_high, n_high) 
        
        #fit integrals per bin
    
        low = integral_fit(m1, m2, b_low, M0, n)         
        high = integral_fit(m1, m2, b_high, M0_high, n_high) 
        
        #scaling
        low *= 1. / int_low * (n_low_rsq)
        high *= 1. / int_high * (n_high_rsq)

        #        print "bin: (%f,%f) low: %f high: %f" % (m1, m2, low, high)

        #fill the histograms
        hist_low_rsq.SetBinContent(ii+1,low)    
        hist_high_rsq.SetBinContent(ii+1,high)

        fit_low_hist = toy_hists_low[ii]
        fit_high_hist = toy_hists_high[ii]
        
        low_func = fit_low_hist.GetFunction("gaus")
        high_func = fit_high_hist.GetFunction("gaus")
        
        low_error = high_error = 1
        
        if fit_low_hist.Integral() > 0: low_error = low_func.GetParameter(2)
        if fit_high_hist.Integral() > 0: high_error = high_func.GetParameter(2)
        
        if low_func.GetParameter(1) < 0: low_error = smallest_68(fit_low_hist)
        if high_func.GetParameter(1) < 0: high_error = smallest_68(fit_high_hist)
        
        #        print "\t\t low_error: %f high_error: %f" % (low_error,high_error)
        
        hist_low_rsq.SetBinError(ii+1, low_error)
        hist_high_rsq.SetBinError(ii+1, high_error)

    return (hist_low_rsq, hist_high_rsq)

def build_tgraph_errors(mu, bin_result):

    mid_points = []
    delta_points = []
    delta_p_points = []
    delta_m_points = []
    zeros = []

    for result in bin_result:
        print result
        bin_l = result[0]
        bin_r = result[1]
        delta = result[2]
        delta_p = result[3]
        delta_m = result[4]
        
        mid  = (bin_r + bin_l) / 2.0

        #append all the values
        mid_points.append(mid)

        delta_points.append(delta)
        delta_p_points.append(abs(delta_p - delta))
        delta_m_points.append(abs(delta - delta_m))
        zeros.append(0)
        
    mid_points_ar = array.array("d",mid_points)
    delta_points_ar = array.array("d",delta_points)
    delta_p_points_ar = array.array("d",delta_p_points)
    delta_m_points_ar = array.array("d",delta_m_points)
    zeros_ar = array.array("d",zeros)

    #build tehg raph
    graph = rt.TGraphAsymmErrors(len(mid_points), mid_points_ar, delta_points_ar, zeros_ar, zeros_ar, delta_m_points_ar, delta_p_points_ar)
    graph.SetTitle("#mu = %.2f" % mu)

    return graph

def average_bins_for_given_mu_results(common_mu_result):

    windows = []
    bin_averaged = []
    #build a blank array for each window
    for ii in bin_windows: windows.append([])

    #go through each result and split them by window
    for ii in common_mu_result:
        window = (ii[0],ii[1])
        index = bin_windows.index(window)
        windows[index].append(ii)

    for ww in windows:
        
        avg_bl = 0
        avg_br = 0
        avg_delta = 0
        avg_delta_p = 0
        avg_delta_m = 0
        
        for toy in ww:
            avg_bl = toy[0]
            avg_br = toy[1]
            avg_delta += toy[2]
            avg_delta_p += toy[3]
            avg_delta_m += toy[4]

        avg_delta = avg_delta / len(ww)
        avg_delta_p = avg_delta_p / len(ww)
        avg_delta_m = avg_delta_m / len(ww)

        bin_averaged.append((avg_bl, avg_br, avg_delta, avg_delta_p, avg_delta_m))

    return bin_averaged

trandom = rt.TRandom()

bins = makebins(.6, 5., .1, .3)

print "USING BINS:", bins

bin_windows = []

#build a list of windows
for bb in range(len(bins)-1): bin_windows.append((bins[bb],bins[bb+1]))

mu_scan = [0] + [.5] * 100 + [1]*100 + [2]*100

#mu_scan = [0] + [.5] + [1] + [2]

parser = OptionParser()

parser.add_option("-f", "--file", dest="filename",
                  help="fit.root file name to analyze FILE",
                  default="razor_output.root",
                  action="store",type="string")

parser.add_option("-d","--data",dest="datafile",
                  help="HggOutput.root data file to use",
                  action="store",type="string")

parser.add_option("-t", "--toys", dest="n_toys",
                  help="Number of Toys to Throw",
                  default="100",
                  action="store",type="int")

parser.add_option("-m", "--mix", dest="mix_file",
                 help="path to signal file to perform mix with",
                 default="no_file",
                 action="store",type="string")

parser.add_option("--outmix", dest="out_mix_file",
                 help="path to where mixed signal file will be placed",
                 default="no_file",
                 action="store",type="string")

parser.add_option("-x", "--xsec", dest="xsec_file",
                 help="file containing signal cross sections",
                 default="no_xsec",
                 action="store",type="string")

parser.add_option("-s", "--samp", dest="samp",
                 help="Sample Selection",
                 default=1,
                 action="store",type="int")


(options, args) = parser.parse_args()

parser.print_help()

mr_min_low = .6
mr_min = .6
mr_cutoff = 5
rsq_cutoff = 5

r0_cut = .01
r1_cut = .02

r0 = r0_cut
r1 = r1_cut

#BUILD THE MIX MODELS (AND CORR. FILE)
mu_signal_pairs = []
file_num = 0
for mu in mu_scan:
    print "Scanning mu %f, file number: %i" % (mu, file_num)

    #generate one file per signal scanned containing the histogram
    name = options.out_mix_file.rstrip(".root")+"mu_%f_%i.root" % (mu,file_num)
    (eff,ns_events, ns_events_plus, ns_events_minus) = build_mixed_model(options.mix_file, options.xsec_file, mu, mr_min, r0_cut, name)

    mu_info = (mu, ns_events, ns_events_plus, ns_events_minus)    
    mu_signal_pairs.append((mu_info,name))

    #increment file number
    file_num +=1
    
#NOW OPEN THE FILES WE WILL USE TO THROW TOYS
print "DATAFILE NAME", options.datafile
data_file = rt.TFile(options.datafile,"READ")

#OPEN THE FIT FILE 
fit_file = rt.TFile(options.filename)
fr = fit_file.Get("Had/independentFR")
fr.Print()
fr_list = fr.floatParsFinal()

#PARSE THE PARAMETERS FROM THE FIT
b_low = fr_list[2].getVal()
n = fr_list[3].getVal()
M0 = fr_list[0].getVal()

b_high = b_low
n_high = n
M0_high = M0

#THE TREE CONTAINING THE DATA
tree = data_file.Get("HggOutput")

#MAKE AND FILL THE LOW AND HIGH RSQ HISTS
#outputfile
out_put_file = rt.TFile("weight_hist.root","RECREATE")
#signal
hist_signal = rt.TH1F("hist_signal","Signal",len(bins)-1,array.array("d",bins))
#estimation
hist_low_rsq = rt.TH1F("hist_low","Low Rsq",len(bins)-1,array.array("d",bins))
hist_high_rsq = rt.TH1F("hist_high","High Rsq",len(bins)-1,array.array("d",bins))
#data
hist_low_rsq_data = rt.TH1F("hist_low_data","Low Rsq Data",len(bins)-1,array.array("d",bins))
hist_high_rsq_data = rt.TH1F("hist_high_data","High Rsq Data",len(bins)-1,array.array("d",bins))

#define the cuts to extract the samples
low_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && PFR^2 < %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r0_cut, r1_cut,options.samp)
high_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r1_cut,options.samp)

low_cut_signal = "(PFMR)>%f && PFR^2 > %f && PFR^2 < %f " % (mr_min,r0_cut, r1_cut)
high_cut_signal = "(PFMR)>%f && PFR^2 > %f" % (mr_min,r1_cut)

#fill data
tree.Draw("PFMR/1000>>hist_low_data",low_cut_data)
tree.Draw("PFMR/1000>>hist_high_data",high_cut_data)

#add and subtract something from the signal histogram
tree.Draw("PFMR/1000>>hist_signal",high_cut_data)
tree.Draw("PFMR/1000>>+hist_signal","-1*("+high_cut_data+")")


n_low_rsq = tree.GetEntries(low_cut_data)
n_high_rsq = tree.GetEntries(high_cut_data)

#final output object
mu_output_tuples = []

#loop over the values of mu
for pair in mu_signal_pairs:    

    print "\n\n\n PERFORMING MU SCAN # %i / %i" % (mu_signal_pairs.index(pair), len(mu_signal_pairs))

    #extract the signal mu info
    (mu, ns_events, ns_events_plus, ns_events_minus) = pair[0]
    name_signal  = pair[1]
    
    signal_file = rt.TFile(name_signal,"READ")
    tree_signal = signal_file.Get("HggOutput")

    #calculate normalizations    
    n_high_rsq_data = hist_high_rsq_data.Integral() + tree_signal.GetEntries(high_cut_signal)    
    
    #get the toys for each bin
    (low_list, high_list, hists_low, hists_high) = get_bin_errors(fr, n_high_rsq_data)
    
    #fit the toy histograms    
    for ii in hists_high: ii.Fit("gaus","q")
    for ii in hists_low: ii.Fit("gaus","q")
    
    #use the fit + toys to determine the errors and values in the two regions
    (hist_low_rsq, hist_high_rsq) = assign_values_from_fit(hist_low_rsq, hist_high_rsq, hists_low, hists_high)    
    
    #loop over the bins in the high rsq_hist
    nbins = hist_high_rsq.GetNbinsX()

    #array of bin results
    bin_results = []
    
    for bin in range(nbins):
        #expected number of events
        exp_events = hist_high_rsq.GetBinContent(bin+1)
        exp_error = hist_high_rsq.GetBinError(bin+1)

        #observed in the data (no signal injection)
        obs_events = hist_high_rsq_data.GetBinContent(bin+1)

        #scale factor to reweight the baseline cross section
        scale_factor_up = scale_factor_down = 0
        if ns_events > 0:
            scale_factor_up = float(ns_events_plus) / ns_events
            scale_factor_down = float(ns_events_minus) / ns_events

        #get the amount of signal in the given bin"
        bin_cut = "(PFMR > %f && PFMR < %f)&&" % (bins[bin], bins[bin+1])
        n_signal_bin = tree_signal.GetEntries(bin_cut+high_cut_signal)

        #observation scaled by the proper error on the cross section
        obs_events_0 = obs_events + n_signal_bin
        obs_events_p = obs_events + n_signal_bin * scale_factor_up
        obs_events_m = obs_events + n_signal_bin * scale_factor_down

        #the various deviations from the expectation up and down
        z = exp_events - obs_events_0
        z_p = exp_events - obs_events_p
        z_m = exp_events - obs_events_m
        
        delta = delta_p = delta_m = -1

        if exp_error > 0:
            delta = z / exp_error
            delta_p = z_p / exp_error
            delta_m = z_m / exp_error
        else:
            "\n\n\nERROR ON BIN IS ZERO\n\n\n" 

        #calculate delta r
        print "BIN: \t\t %f" % bin
        print "EXP EVENTS:\t %f" % exp_events        
        print "EXP ERROR:\t %f" % exp_error
        print "OBS EVENTS:\t %f" % obs_events_0
        print "DELTA: \t %f" % delta
        print "\n\n\n"
        
        #subtract off the signal when you're done with this mu
        tree_signal.Draw("PFMR>>+hist_signal", "-1*("+high_cut_signal+")")
            
        #bin range with middle value and error on delta
        output_tuple = (bins[bin], bins[bin+1], delta, delta_p, delta_m)
        
        bin_results.append(output_tuple)
        
    mu_output_tuples.append((mu, bin_results))

final_output = rt.TFile("tgraph_scan.root", "RECREATE")

graphs = []
tmultigraph = rt.TMultiGraph("total_graph","total")

#perfrom an average of repeated mu values
zipped_mu = zip(*mu_output_tuples)
#parse the mus and bin_results
mus = zipped_mu[0]
bin_res = zipped_mu[1]

print "BIN RES AFTER ZIP:", bin_res
print "BIN REST ELEMENT 0:", bin_res[0]
print "MUS AFTER ZIP:", mus
#checked mu values
mu_checked = []

averaged_result = None
average_result_list = []
for ii in range(len(mus)):
    mu = mus[ii]
    averaged_result = None


    #if there is a repeated
    if mus.count(mu) > 1 and mu not in mu_checked:
        print "Found repeated mu value....", mu, " with ", mus.count(mu), " occurances"

        #add it to the checked list
        mu_checked.append(mu)

        common_mu_results = []
        #check every element of the list and keep the bin results
        for jj in mus:
            if mu == jj:
                for indv_bin in bin_res[mus.index(jj)]:
                    common_mu_results.append(indv_bin)

        
        averaged_result = average_bins_for_given_mu_results(common_mu_results)
                
    #no repeated values    
    elif mu not in mu_checked:
        averaged_result = bin_res[ii]    
    else:
        continue
    
    average_result_list.append((mu,averaged_result, mus.count(mu)))

for result in average_result_list:
    (mu,averaged_result, counts) = result
    
    canvas = rt.TCanvas("averaged_canvas_%.2f_npoints_%i" % (mu, counts))
    graph = build_tgraph_errors(mu, averaged_result)

    xaxis = graph.GetXaxis()
    yaxis = graph.GetYaxis()

    xaxis.SetTitle("M_{R} Bin Center [TeV]")
    yaxis.SetTitle("#Delta = (N_{exp} - N_{obs}) / #sigma_{exp}")
    graph.SetMarkerStyle(21)
    graph.SetMarkerColor(rt.kBlack)
    graph.Draw("ALP")
    
    canvas.Write()
    graph.Write()

    graphs.append(graph)    

canvas_tot = rt.TCanvas("canvas_total")
#draw the first

colors = [rt.kBlack, rt.kRed+1, rt.kBlue, rt.kGreen+3, rt.kMagenta+1, rt.kCyan, rt.kOrange]
#draw the rest
for graph in graphs:
    color = colors[graphs.index(graph)]
    graph.SetLineColor(color)
    graph.SetLineWidth(2)
    graph.SetFillColor(0)
    graph.SetMarkerColor(color)
    tmultigraph.Add(graph)

#parse the title and xsec again
signal_name = options.mix_file.split("/")[-1]
msq = signal_name.split("_")[1]
mgl = signal_name.split("_")[2]   
(xsec,xsec_plus,xsec_minus) = get_xsec(options.xsec_file, msq, mgl)
title = "Signal Injection m_{sq} = %i GeV m_{gl} = %i GeV m_{#tilde{#chi}} = 375 GeV #sigma_{LO} = %2.2f fb" % (int(msq), int(mgl), float(xsec)*1000)

tmultigraph.SetTitle(title)
tmultigraph.Draw("ALP")

leg = canvas_tot.BuildLegend()
leg.SetFillColor(0)
canvas_tot.SetGridy(1)
xaxis = tmultigraph.GetXaxis()
yaxis = tmultigraph.GetYaxis()
xaxis.SetTitle("M_{R} Bin Center [TeV]")
yaxis.SetTitle("#Delta = (N_{exp} - N_{obs}) / #sigma_{exp}")


tmultigraph.Write()
canvas_tot.Write()
final_output.Close()


