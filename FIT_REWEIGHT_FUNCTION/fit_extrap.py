from  optparse  import OptionParser
import ROOT as rt
import RootTools
import array, os
import random

#speed up the drawing by supressing x11
rt.gROOT.SetBatch(True);

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
def build_mixed_model(signal_file, xsec_file, mu, mr_min, rsq_min, outfile):

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

    outfile = rt.TFile(outfile,"RECREATE")
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

trandom = rt.TRandom()


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


parser.add_option("-p", "--print", dest="print_pdf",
                 help="print pdfs",
                 default=False,
                 action="store_true")

parser.add_option("-n", "--nodraw", dest="draw",
                 help="dont draw the result",
                 default=False,
                 action="store_true")

parser.add_option("-e", "--debug", dest="debug",
                 help="debug",
                 default=False,
                 action="store_true")

parser.add_option("-r", "--rioff", dest="ri",
                 help="offset for Rsq",
                 default=0.0,
                 action="store",type="float")

parser.add_option("-s", "--samp", dest="samp",
                 help="Sample Selection",
                 default=1,
                 action="store",type="int")

parser.add_option("--mrmin", dest="mrmin",
                 help="minimum mr for the low rsq fit",
                 default=.6,
                 action="store",type="float")

parser.add_option("--rsq1", dest="rsq1",
                 help="minimum rsq for the low rsq fit",
                 default=.01,
                 action="store",type="float")

parser.add_option("--rsq2", dest="rsq2",
                 help="max rsq for the low rsq fit",
                 default=.02,
                 action="store",type="float")

parser.add_option("--mu", dest="mu",
                 help="if mix is being performed include value of mu being scanned",
                 default=1.0,
                 action="store",type="float")

(options, args) = parser.parse_args()

parser.print_help()

bins = makebins(options.mrmin, 5., .1, .3)

#FIRST SET THE BOUNDS IN MR AND RSQ
mr_min_low = options.mrmin
mr_min = bins[0]
mr_cutoff = 5
rsq_cutoff = 5

r0_cut = options.rsq1
r1_cut = options.rsq2 

r0 = r0_cut - options.ri
r1 = r1_cut - options.ri



#IF WE ARE GENREATING A MIX OPEN THAT FILE FIRST AND MAKE IT
do_mix =  options.mix_file != "no_file"
tree_signal = None
if do_mix:
    (eff,ns_events, ns_events_plus, ns_events_minus) = build_mixed_model(options.mix_file, options.xsec_file, options.mu, mr_min, r0_cut, options.out_mix_file)
    mix_file_output = rt.TFile(options.out_mix_file,"READ")
    tree_signal = mix_file_output.Get("HggOutput")

#NOW OPEN THE FILES WE WILL USE TO THROW TOYS

print "DATAFILE NAME", options.datafile
data_file = rt.TFile(options.datafile,"READ")

#OPEN THE FIT FILE OUTPUT LAST
fit_file = rt.TFile(options.filename)

fr = fit_file.Get("Had/independentFR")
fr.Print()
fr_list = fr.floatParsFinal()
    
b_low = fr_list[2].getVal()#9.63 #* 1 /( pow(r0, 1/n))# + pow(r1, 1/n))
n = fr_list[3].getVal()#58.291
M0 = fr_list[0].getVal()#-.3477


tree = data_file.Get("HggOutput")
print "TREE CHECK: ", tree.GetEntries()


#hist_high_pred = rt.TH1F("hist_high_pred","High Rsq Pred",len(bins)-1,array.array("d",bins))
hist_signal = rt.TH1F("hist_signal","Signal",len(bins)-1,array.array("d",bins))
hist_low_rsq_data = rt.TH1F("hist_low_data","Low Rsq Data",len(bins)-1,array.array("d",bins))
hist_high_rsq_data = rt.TH1F("hist_high_data","High Rsq Data",len(bins)-1,array.array("d",bins))
hist_high_rsq_signal = rt.TH1F("hist_high_signal","High Rsq Signal",len(bins)-1,array.array("d",bins))

low_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && PFR^2 < %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r0_cut, r1_cut,options.samp)
high_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r1_cut,options.samp)

low_cut_signal = "(PFMR)>%f && PFR^2 > %f && PFR^2 < %f " % (mr_min,r0_cut, r1_cut)
high_cut_signal = "(PFMR)>%f && PFR^2 > %f" % (mr_min,r1_cut)

print low_cut_data
print high_cut_data

n_low_rsq = tree.GetEntries(low_cut_data)
n_high_rsq = tree.GetEntries(high_cut_data)

if do_mix:
    n_low_rsq += tree_signal.GetEntries(low_cut_signal)
    n_high_rsq += tree_signal.GetEntries(high_cut_signal)

print "nlow: ", n_low_rsq
print "nhigh: ", n_high_rsq

b_high = b_low #*pow(r1,1/n)
n_high = n
M0_high = M0


def Gamma(a,x):
    val = rt.TMath.Gamma(a)*rt.Math.inc_gamma_c(a,x)
    if options.debug: print "gamma", val
    return val

def Gfun(m,r,b,n):
    val =  n/pow(b*n,n)*(Gamma(n,b*n*pow(m*r,1/n)))
    if options.debug: print "gfun", val
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
        b_high_toy = b_toy #* pow(r1,1/n_toy)
        n_high_toy = n_toy
        M0_high_toy = M0_toy

        if options.debug: print "b", b_high_toy, "n",n_high_toy

        #find the fit normalization
        if options.debug: print "\ndoing low integral.."
        int_low = integral_fit(mr_min_low, mr_cutoff, b_toy, M0_toy, n_toy)
        if options.debug: print "\ndoing high integral.."
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

tree.Draw("PFMR/1000>>hist_low_data",low_cut_data)
tree.Draw("PFMR/1000>>hist_high_data",high_cut_data)

if do_mix:
    tree_signal.Draw("PFMR>>hist_signal",high_cut_signal)
    #add the signal to the high rsq events
    tree_signal.Draw("PFMR>>+hist_high_data",high_cut_signal)



file = rt.TFile("weight_hist.root","RECREATE")
    
hist_low_rsq = rt.TH1F("hist_low","Low Rsq",len(bins)-1,array.array("d",bins))
hist_high_rsq = rt.TH1F("hist_high","High Rsq",len(bins)-1,array.array("d",bins))
hist_ratio = rt.TH1F("hist_ratio","Ratio",len(bins)-1,array.array("d",bins))
hist_difference = rt.TH1F("hist_difference","Difference %",len(bins)-1,array.array("d",bins))

n_high_rsq_data = hist_high_rsq_data.Integral()

print "N_HIGH_RSQ_DATA", n_high_rsq_data

#get the errors
(low_list, high_list, hists_low, hists_high) = get_bin_errors(fr, n_high_rsq_data)
(low_errors, high_errors) = get_sigma(low_list,high_list)

#fit the histograms
for ii in hists_low:    
    print "\n",ii
    ii.Fit("gaus")    
    print ii.GetFunction("gaus")
for ii in hists_high:
    ii.Fit("gaus")

#integrate each bin and assign errors
for ii in range(len(bins) - 1):
    m1 = bins[ii]
    m2 = bins[ii+1]

    #normalization integrals
    if options.debug: print "doing low integral.."    
    int_low = integral_fit(mr_min, mr_cutoff, b_low, M0, n)
    if options.debug: print "doing high integral.."
    int_high = integral_fit(mr_min, mr_cutoff, b_high, M0_high, n_high) 
    
    #fit integrals per bin

    low = integral_fit(m1, m2, b_low, M0, n)         
    high = integral_fit(m1, m2, b_high, M0_high, n_high) 

    #scaling
    low *= 1. / int_low * (n_low_rsq)
    high *= 1. / int_high * (n_high_rsq)

    print "bin: (%f,%f) low: %f high: %f" % (m1,m2,low,high)

    if low != 0:
        ratio = (high/low)
    else:
        ratio = -999
    #fill the histograms
    hist_low_rsq.SetBinContent(ii+1,low)    
    hist_high_rsq.SetBinContent(ii+1,high)
    hist_ratio.SetBinContent(ii+1,ratio)

    fit_low_hist = hists_low[ii]
    fit_high_hist = hists_high[ii]
    
    low_func = fit_low_hist.GetFunction("gaus")
    high_func = fit_high_hist.GetFunction("gaus")

    low_error = high_error = 1
    
    if fit_low_hist.Integral() > 0: low_error = low_func.GetParameter(2)
    if fit_high_hist.Integral() > 0: high_error = high_func.GetParameter(2)

    if low_func.GetParameter(1) < 0: low_error = smallest_68(fit_low_hist)
    if high_func.GetParameter(1) < 0: high_error = smallest_68(fit_high_hist)

    print "\t\t low_error: %f high_error: %f" % (low_error,high_error)

    hist_low_rsq.SetBinError(ii+1, low_error)
    hist_high_rsq.SetBinError(ii+1, high_error)


                     
#clone the low rsq
hist_ratio = hist_low_rsq.Clone("hist_ratio")
hist_low_pred = hist_low_rsq.Clone("hist_low_pred")
hist_high_pred = hist_high_rsq.Clone("hist_high_pred")
hist_difference = hist_high_rsq_data.Clone("hist_difference")

#there are no errors on the number of observed data
for ii in range(len(bins)-1):
    hist_difference.SetBinError(ii+1,0)

#Set the errors of the ratios
#hist_ratio.Sumw2()
#hist_high_pred.Sumw2()

#hist_high_rsq.Sumw2()
#hist_low_rsq.Sumw2()
hist_difference.Sumw2()

#calculate the prediction
hist_ratio.Divide(hist_high_rsq)
#hist_high_pred.Divide(hist_ratio)
hist_difference.Add(hist_high_pred,-1.)
        

hist_ratio.SetFillColor(rt.kOrange+1)
hist_ratio.SetLineColor(rt.kOrange+1)
hist_ratio.SetFillStyle(3001)
hist_ratio.GetXaxis().SetTitle("M_{R} [TeV]")
hist_ratio.GetYaxis().SetTitle("Low R^{2} / High R^{2}")

hist_ratio.GetXaxis().SetTitleSize(.07)
hist_ratio.GetYaxis().SetTitleSize(.07)


hist_high_pred.GetXaxis().SetTitle("M_{R} [TeV]")
hist_high_pred.GetYaxis().SetTitle("N Events")

hist_high_pred.GetXaxis().SetTitleSize(.07)
hist_high_pred.GetYaxis().SetTitleSize(.07)

hist_high_rsq_data.SetMarkerStyle(14)
hist_high_rsq_data.GetYaxis().SetTitle("N Events")

hist_difference.Write()
hist_low_pred.Write()
hist_low_rsq_data.Write()
hist_high_pred.Write()
hist_high_rsq_data.Write()

hist_low_rsq.Write()
hist_high_rsq.Write()
hist_ratio.Write()
hist_signal.Write()

#loop over all the low bins and setup the histograms
for ii in range(len(hists_low)):
    max_val = max(low_list[ii])
    min_val = min(low_list[ii])
    
    hists_low[ii].GetXaxis().SetTitle("N Events %2.1f < M_{R} (TeV) < %2.1f" % (bins[ii],bins[ii+1]))
    hists_low[ii].GetYaxis().SetTitle("N Toys")
    hists_low[ii].GetXaxis().SetTitleSize(.07)
    hists_low[ii].GetYaxis().SetTitleSize(.07)
    hists_low[ii].Write()

#loop over high bins and setup histograms
for ii in range(len(hists_high)):
    max_val = max(high_list[ii])
    min_val = min(high_list[ii])

    hists_high[ii].GetXaxis().SetTitle("N Events %2.1f < M_{R} (TeV) < %2.1f" % (bins[ii],bins[ii+1]))

    hists_high[ii].GetYaxis().SetTitle("N Toys")
    
    hists_high[ii].GetXaxis().SetTitleSize(.07)
    hists_high[ii].GetYaxis().SetTitleSize(.07)
    
    hists_high[ii].Write()


for ii in range(len(hists_high)):

    gaus = hists_high[ii].GetFunction("gaus")

    gaus_mean = 0
    gaus_sigma = 1

    if hists_high[ii].Integral() > 0 :
        gaus_mean = gaus.GetParameter(1)
        gaus_sigma = gaus.GetParameter(2)

    error = gaus_sigma

    true_val = hist_high_rsq_data.GetBinContent(ii+1)

    rt.gStyle.SetLabelSize(0.05,"xyz");
  
    rt.gStyle.SetPadTopMargin(0.1);
    rt.gStyle.SetPadRightMargin(0.10);
    rt.gStyle.SetPadBottomMargin(0.2);
    rt.gStyle.SetPadLeftMargin(0.15);

    canvas = rt.TCanvas("canvas_bin%i"%ii,"Bin %i" %ii, 1024,768)

    max_val = max(high_list[ii])
    min_val = min(high_list[ii])
    

    p1 = gaus_mean + error
    p2 = gaus_mean + 2 * error
    p3 = gaus_mean + 3 * error
    
    m1 = gaus_mean - error
    m2 = gaus_mean - 2 * error
    m3 = gaus_mean - 3 * error


    absmin = min([min_val, m3])
    absmax = max([max_val, max_val, p3])

    max_hist = hists_high[ii].GetMaximum()

    graph1s = None
    graph2s = None
    
    left_win1s = max(0,m1)
    right_win1s = p1

    left_win2s = max(0,m2)
    right_win2s = p2

    left_win3s = max(0,m3)
    right_win3s = p3
        
    array1s_x = array.array("d",[left_win1s,left_win1s,right_win1s,right_win1s])
    array2s_x = array.array("d",[left_win2s,left_win2s,right_win2s,right_win2s])
    array3s_x = array.array("d",[left_win3s,left_win3s,right_win3s,right_win3s])

    array1s_y = array.array("d",[2*max_hist/10.,3.*max_hist/10.,3.*max_hist/10.,2*max_hist/10.])
    array2s_y = array.array("d",[max_hist/10.,2.*max_hist/10.,2.*max_hist/10.,max_hist/10.])
    array3s_y = array.array("d",[0,max_hist/10.,max_hist/10.,0])



    fill_1s = rt.TGraph(4,array1s_x,array1s_y)
    fill_2s = rt.TGraph(4,array2s_x,array2s_y)
    fill_3s = rt.TGraph(4,array3s_x,array3s_y)

    fill_1s.SetFillColor(rt.kRed)
    fill_2s.SetFillColor(rt.kBlue)    
    fill_3s.SetFillColor(rt.kGreen)
    
    fill_1s.SetFillStyle(3004)
    fill_2s.SetFillStyle(3005)
    fill_3s.SetFillStyle(3004)    

    print absmin,absmax
    
    hists_high[ii].GetXaxis().SetRangeUser(absmin, absmax)
    #    hists_high[ii].SetMarkerStyle(20)
    #    hists_high[ii].SetFillStyle(3004)


    hists_high[ii].Draw()
    #pois = rt.TF1("poiss2","%f * TMath::Poisson(x,%f)" % (norm,true_val),absmin,absmax);
    #pois.SetLineColor(rt.kBlack)
    #pois.SetLineStyle(9)
    
    #    pois.Draw("same")

    line = rt.TLine(true_val, 0, true_val, max_hist)
    line.SetLineWidth(3)
    line.Draw("same")

    fill_1s.Draw("same")
    fill_2s.Draw("same")
    fill_3s.Draw("same")

    #top of the gaussian
    highest_point = 0
    if gaus_mean != 0: highest_point = gaus.Eval(gaus_mean)
    
    line_mean = rt.TLine(gaus_mean, 0, gaus_mean, highest_point)
    line_mean.SetLineColor(rt.kRed)
    line_mean.SetLineWidth(3)
    line_mean.Draw("same")

    #    hists_high[ii].SetLineStyle(2)
    hists_high[ii].SetFillColor(rt.kBlue)
    hists_high[ii].SetFillStyle(3005)
    hists_high[ii].SetLineWidth(1)
    hists_high[ii].SetLineColor(rt.kBlue)

    legend = rt.TLegend(.415,.394,.86,.738)
#    legend.AddEntry(pois,"Poisson Distribution on Data","l")
    legend.AddEntry(hists_high[ii],"Toy Distribution","f")
    legend.AddEntry(gaus,"Gaussian Fit","l")
    legend.AddEntry(line, "N Observed","l")
#    legend.AddEntry(fill_1s,"1 #sigma confidence","f")
#    legend.AddEntry(fill_2s,"2 #sigma confidence","f")
#    legend.AddEntry(fill_3s,"99.7% confidence","f")

    legend.SetFillColor(0)
    legend.SetLineColor(0)

    if not options.print_pdf: legend.Draw("same")

    delta_text= rt.TPaveText(.125,1,.9,.9,"NDC");

    delta = -9
    if gaus_mean > 0 :
        delta = abs(gaus_mean - true_val) / gaus_sigma
    else:
        delta = float(true_val) / float(smallest_68(hists_low[ii]))

    delta_text.SetTextFont(42);
    delta_text.SetTextSize(0.04);
    delta_text.SetFillColor(0); 
    delta_text.SetTextAlign(12);
    if gaus_mean > 0:
        delta_text.AddText("#Delta = (N data - #mu_{gaus}) / #sigma_{gaus} = (%2.0f - %2.0f) / %2.2f = %2.2f" % (max(gaus_mean,0),true_val, gaus_sigma,delta));
    else:
        delta_text.AddText("#Delta = N_{obs} /  (.68 range) = %2.2f" % delta)

    delta_text.Draw("same")
    
    canvas.Write()

    if options.print_pdf:
        canvas.SaveAs("canvas_%i.pdf" % ii)

for ii in range(len(hists_low)):
    max_val = max(low_list[ii])
    min_val = min(low_list[ii])

    hists_low[ii].GetXaxis().SetTitle("N Events %2.1f < M_{R} (TeV) < %2.1f" % (bins[ii],bins[ii+1]))

    hists_low[ii].GetYaxis().SetTitle("N Toys")
    
    hists_low[ii].GetXaxis().SetTitleSize(.07)
    hists_low[ii].GetYaxis().SetTitleSize(.07)
    
    hists_low[ii].Write()


for ii in range(len(hists_low)):

    gaus = hists_low[ii].GetFunction("gaus")

    gaus_mean = 0
    gaus_sigma = 1

    if hists_low[ii].Integral() > 0 :
        gaus_mean = gaus.GetParameter(1)
        gaus_sigma = gaus.GetParameter(2)

    error = gaus_sigma

    true_val = hist_low_rsq_data.GetBinContent(ii+1)

    rt.gStyle.SetLabelSize(0.05,"xyz");
  
    rt.gStyle.SetPadTopMargin(0.1);
    rt.gStyle.SetPadRightMargin(0.10);
    rt.gStyle.SetPadBottomMargin(0.2);
    rt.gStyle.SetPadLeftMargin(0.15);

    canvas = rt.TCanvas("low_canvas_bin%i"%ii,"Bin %i" %ii, 1024,768)

    max_val = max(low_list[ii])
    min_val = min(low_list[ii])
    

    p1 = gaus_mean + error
    p2 = gaus_mean + 2 * error
    p3 = gaus_mean + 3 * error
    
    m1 = gaus_mean - error
    m2 = gaus_mean - 2 * error
    m3 = gaus_mean - 3 * error

    absmin = min([min_val, m3])
    absmax = max([max_val, max_val, p3])

    max_hist = hists_low[ii].GetMaximum()

    graph1s = None
    graph2s = None
    
    left_win1s = max(0,m1)
    right_win1s = p1

    left_win2s = max(0,m2)
    right_win2s = p2

    left_win3s = max(0,m3)
    right_win3s = p3
        
    array1s_x = array.array("d",[left_win1s,left_win1s,right_win1s,right_win1s])
    array2s_x = array.array("d",[left_win2s,left_win2s,right_win2s,right_win2s])
    array3s_x = array.array("d",[left_win3s,left_win3s,right_win3s,right_win3s])

    array1s_y = array.array("d",[2*max_hist/10.,3.*max_hist/10.,3.*max_hist/10.,2*max_hist/10.])
    array2s_y = array.array("d",[max_hist/10.,2.*max_hist/10.,2.*max_hist/10.,max_hist/10.])
    array3s_y = array.array("d",[0,max_hist/10.,max_hist/10.,0])



    fill_1s = rt.TGraph(4,array1s_x,array1s_y)
    fill_2s = rt.TGraph(4,array2s_x,array2s_y)
    fill_3s = rt.TGraph(4,array3s_x,array3s_y)

    fill_1s.SetFillColor(rt.kRed)
    fill_2s.SetFillColor(rt.kBlue)    
    fill_3s.SetFillColor(rt.kGreen)
    
    fill_1s.SetFillStyle(3004)
    fill_2s.SetFillStyle(3005)
    fill_3s.SetFillStyle(3004)    

    print absmin,absmax
    
    hists_low[ii].GetXaxis().SetRangeUser(absmin, absmax)
    #    hists_low[ii].SetMarkerStyle(20)
    #    hists_low[ii].SetFillStyle(3004)


    hists_low[ii].Draw()

    line = rt.TLine(true_val, 0, true_val, max_hist)
    line.SetLineWidth(3)
    line.Draw("same")

    fill_1s.Draw("same")
    fill_2s.Draw("same")
    fill_3s.Draw("same")

    #top of the gaussian
    lowest_point = 0
    if gaus_mean != 0: lowest_point = gaus.Eval(gaus_mean)
    
    line_mean = rt.TLine(gaus_mean, 0, gaus_mean, lowest_point)
    line_mean.SetLineColor(rt.kRed)
    line_mean.SetLineWidth(3)
    line_mean.Draw("same")

    #    hists_low[ii].SetLineStyle(2)
    hists_low[ii].SetFillColor(rt.kBlue)
    hists_low[ii].SetFillStyle(3005)
    hists_low[ii].SetLineWidth(1)
    hists_low[ii].SetLineColor(rt.kBlue)

    legend = rt.TLegend(.415,.394,.86,.738)
#    legend.AddEntry(pois,"Poisson Distribution on Data","l")
    legend.AddEntry(hists_low[ii],"Toy Distribution","f")
    legend.AddEntry(gaus,"Gaussian Fit","l")
    legend.AddEntry(line, "N Observed","l")
#    legend.AddEntry(fill_1s,"1 #sigma confidence","f")
#    legend.AddEntry(fill_2s,"2 #sigma confidence","f")
#    legend.AddEntry(fill_3s,"99.7% confidence","f")

    legend.SetFillColor(0)
    legend.SetLineColor(0)

    if not options.print_pdf: legend.Draw("same")

    delta_text= rt.TPaveText(.125,1,.9,.9,"NDC");

    print gaus_mean , true_val, gaus_sigma , "delta", delta

    delta = -9
    if gaus_mean > 0 :
        delta = abs(gaus_mean - true_val) / gaus_sigma
    else:
        delta = float(true_val) / float(smallest_68(hists_low[ii]))

    delta_text.SetTextFont(42);
    delta_text.SetTextSize(0.04);
    delta_text.SetFillColor(0); 
    delta_text.SetTextAlign(12);

    if gaus_mean > 0:
        delta_text.AddText("#Delta = (N data - #mu_{gaus}) / #sigma_{gaus} = (%2.0f - %2.0f) / %2.2f = %2.2f" % (max(gaus_mean,0),true_val, gaus_sigma,delta));
    else:
        delta_text.AddText("#Delta = N_{obs} /  (.68 containment window) = %2.2f" % delta)

    delta_text.Draw("same")
    
    canvas.Write()

    if options.print_pdf:
        canvas.SaveAs("canvas_%i.pdf" % ii)


file.Close()

if do_mix:
    os.system("cp weight_hist.root weight_hist_mix.root")
if not options.draw:
    os.system("root -l -c ~/rootlogon.C ratio_extrap.C")
