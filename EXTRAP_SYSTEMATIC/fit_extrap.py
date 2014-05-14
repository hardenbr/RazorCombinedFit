from optparse  import OptionParser
import ROOT as rt
import RootTools
import array, os
import random, math


class error_summary:
    def __init__(self, list_of_toy_dists):
        self.left_vals = []
        self.right_vals = []
        self.nobs = []
        self.nexp = []
        self.delta = []

        for ii in list_of_toy_dists:
            self.add_toy_dist(ii)

    def add_toy_dist(self, dist):
        self.left_vals.append(dist.get_l68())
        self.right_vals.append(dist.get_r68())
        self.nobs.append(dist.get_nobs())
        self.nexp.append(dist.get_nexp())
        self.delta.append(dist.get_delta())

    def get_l68(self): return self.left_vals
    def get_r68(self): return self.right_vals
    def get_nobs(self): return self.nobs
    def get_nexp(self): return self.nexp
    def get_delta(self): return self.delta

    def clear_vals(self):
        self.left_vals = []
        self.right_vals = []
        self.nobs = []
        self.nexp = []
        self.delta = []

class fit_toy_summary:

    def __init__(self):
        self.mr0 = []
        self.n = []
        self.b = []
        self.ntot = []

    #add a single toy
    def add_toy_dist(self, mr0, n, b, ntot):
        self.mr0.append(mr0)
        self.n.append(n)
        self.b.append(b)
        self.ntot.append(ntot)

    #retrieve a list of fit values
    def get_mr0(self): return self.mr0
    def get_n(self): return self.n
    def get_b(self): return self.b
    def get_ntot(self): return self.ntot

    def get_ntoys(self): return len(self.mr0)

#container class for a single distribution
class toy_distribution:
    def __init__(self, l, r, ne, no, d):
        (self.l, self.r, self.ne, self.no, self.d)  = (l,r,ne,no,d)

    #retreive values of the container
    def get_l68(self):return self.l
    def get_r68(self):return self.r
    def get_nexp(self):return self.ne
    def get_nobs(self):return self.no
    def get_delta(self):return self.d

#build graph and exp_up and exp_down
def build68graph(error_summary, bins, name):

    nobs = error_summary.get_nobs()
    x = error_summary.get_nexp()
    r_x = error_summary.get_r68()
    l_x = error_summary.get_l68()

    #build the sigma up and sigma down histograms
    hist_exp_up = rt.TH1F("hist_%s_up" % name,"68% pseudo-experiment Window",len(bins)-1,array.array("d",bins))
    hist_exp_down = rt.TH1F("hist_%s_down" % name,"expectation sigma down",len(bins)-1,array.array("d",bins))

    hist_err_up = rt.TH1F("hist_%s_err_up" % name, "", len(bins)-1, array.array("d",bins))
    hist_err_down = rt.TH1F("hist_%s_err_down" % name, "", len(bins)-1, array.array("d",bins))

    for bin in range(len(r_x)): hist_exp_up.SetBinContent(bin+1, r_x[bin])
    for bin in range(len(l_x)): hist_exp_down.SetBinContent(bin+1, l_x[bin])

    for bin in range(len(x)): hist_err_up.SetBinContent(bin+1, (r_x[bin] - x[bin]) / x[bin])
    for bin in range(len(x)): hist_err_down.SetBinContent(bin+1, (x[bin] - l_x[bin]) / x[bin])

    print "nobs",nobs,"\n"
    print "x", x,"\n"
    print "r_x",r_x,"\n"
    print "l_x",l_x,"\n"
    
    width = []
    y = [0] * len(x)
    
    # get the width sizes
    for ii in range(len(bins)-1):
        w = bins[ii+1] - bins[ii]
        if w < 0: print "ERROR: BINS MUST BE INCREASING WIDTH"
        width.append(w / 2.0)

    ex_up = []
    ex_down = []
    obs_val = []
    delta = []

    for ii in range(len(x)):
        win_size_r = (r_x[ii] - x[ii])
        win_size_l = (x[ii] - l_x[ii])
        win_size = win_size_r + win_size_l + 1

#        print r_x[ii], l_x[ii]
#        print "win_size", win_size

#        if win_size_r == 0 or win_size_l == 0: win_size = 1

        if x[ii] > 1:

            exp_up_val = 2*(r_x[ii] - x[ii]) / win_size
            exp_down_val = 2*(x[ii] - l_x[ii]) / win_size
            obs_value = 2*(nobs[ii] - x[ii]) / win_size
            if nobs[ii] > x[ii]:
                obs_value = 2*(nobs[ii] - x[ii]) / win_size
                
            ex_up.append(exp_up_val)
            ex_down.append(exp_down_val)
            obs_val.append(obs_value)
            delta.append(obs_val[ii])

            print "----WINDOW %s BIN %i----" % (name,ii)
            print "win_size_l", win_size_l
            print "win_size_r", win_size_r
            print "ex_up", ex_up[ii]
            print "ex_down", ex_down[ii]
            print "obs_val", obs_val[ii]
            print "delta", delta[ii]
            print "\n"

        else:
            ex_down.append(0)
            ex_up.append(1)
            obs_value = (nobs[ii] - x[ii]) / (win_size - 1)
            if nobs[ii] == 0:
                obs_value = 0
            
            obs_val.append(obs_value)
            
            delta.append((nobs[ii] - x[ii]) / (win_size - 1))


    adjusted_bins = []

    for ii in bins: adjusted_bins.append(ii)
    
    for ii in range(len(width)):
        adjusted_bins[ii] += width[ii]

    bins_ar = array.array("f",adjusted_bins)
    x_ar = array.array("f", x)
    y_ar = array.array("f", y)
    ey_u_ar = array.array("f", ex_up)
    ey_d_ar = array.array("f", ex_down)
    width_ar = array.array("f", width)
    obs_val_ar = array.array("f", obs_val) 


    mg = rt.TMultiGraph(name,name)
    
    gr = rt.TGraphAsymmErrors(len(x),  bins_ar, y_ar, width_ar, width_ar, ey_d_ar, ey_u_ar)
    gr.SetLineColor(rt.kBlack)
    gr.SetLineWidth(2)
    gr.SetFillColor(rt.kBlue-9)
#    gr.SetFillStyle(3000)
    
    gr2 = rt.TGraphAsymmErrors(len(x), bins_ar, obs_val_ar, y_ar, y_ar, y_ar,y_ar)
    gr2.SetMarkerStyle(20)
    gr2.SetMarkerSize(1.4)
    

    
    return (hist_exp_up, hist_exp_down, gr, gr2, hist_err_up, hist_err_down)

#exponentially increasing bin sizes
def makebins(start_,end_,inc_,inc_inc_):
    bin = start_
    inc = inc_
    list = []
    while True:
        list.append(bin)
        bin+=inc
        inc*= (1+inc_inc_)

        if bin > end_:
            list.append(bin)
            return list


#speed up the drawing by supressing x11
rt.gROOT.SetBatch(True);

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

parser.add_option("--t5gg", dest="t5gg",
                 help="flag to run t5gg limit",
                 action="store_true",default=False)

(options, args) = parser.parse_args()

parser.print_help()

#FIRST SET THE BOUNDS IN MR AND RSQ
mr_min_low = options.mrmin
mr_min = mr_min_low


r0_cut = options.rsq1
r1_cut = options.rsq2 

r0 = r0_cut
r1 = r1_cut

#bins = [0.6, 0.7, 0.82, 0.964, 1.1368,  2.2499084799999998, 2.679890176, 3.1958682111999996, 3.8150418534399995, 4.558050224127999] #combination bin
#bins = [0.6, 0.7, 0.82, 0.964, 1.1368, 1.34416, 1.592992, 1.8915904, 2.2499084799999998,  4.558050224127999]
bins = makebins(mr_min, 4.5, .1, .2)
print bins
fit_toy_summary = fit_toy_summary()
mr_cutoff = bins[-1]
rsq_cutoff = 1

def lowest_filled_bin(hist):
    bin = 1
    while True:
        if hist.GetBinContent(bin) != 0:
            return bin
        bin += 1

def highest_filled_bin(hist):
    bin = hist.GetNbinsX()    
    while True:
        if hist.GetBinContent(bin) != 0: return bin        
        bin -= 1

def get_smallest68_hist(exp, hist):
    width_left = 0
    width_right = 0

    total_events = hist.Integral()

    low_bin = lowest_filled_bin(hist)
    high_bin = highest_filled_bin(hist)
    exp_bin = hist.GetXaxis().FindBin(exp)#hist.GetMaximumBin()###
    zero_bin = hist.GetXaxis().FindBin(0)

    counter = 0 
    while True:
        if counter > 300: return (100,0,100)
        counter+=1
        left_edge = exp_bin - width_left
        right_edge = exp_bin + width_right

        if left_edge < 1 : left_edge = 1  #cant go below zero
        if right_edge > high_bin: right_edge = high_bin #cant go outside bounds

        integral = float(hist.Integral(left_edge, right_edge))
        
        if left_edge == right_edge:
            integral = hist.GetBinContent(left_edge)

        containment = integral / total_events

        if containment > .68:

            left_val = hist.GetBinLowEdge(left_edge) 
            right_val = hist.GetBinLowEdge(right_edge) + 1

            #interpolate the exact place in the bin
            if hist.GetBinLowEdge(right_edge) < 10:
            
                integral_minus_bin = float(hist.Integral(left_edge, right_edge - 1))
                containment_minus_this_bin =  integral_minus_bin / total_events 
                right_val_minus_1 = hist.GetBinLowEdge(right_edge-1) + 1                     
                #correct the right value
                print "contain minus", containment_minus_this_bin
                print "right val", right_val
                print "right val minus", right_val_minus_1
                
                right_val = (.68 - containment_minus_this_bin) * (right_val - right_val_minus_1) / (containment - containment_minus_this_bin) + right_val_minus_1

                if right_edge == 1:
                    right_val = (.68 / containment) 

            #once we have containment return a tuple ( sigma_eff, left value, right_value)            
            if left_val == right_val:
                if exp > .1:
                    return (1, hist.GetBinLowEdge(left_edge),hist.GetBinLowEdge(left_edge+1)) #special case of 68 existing in one bin
                else:
                    return (1, hist.GetBinLowEdge(left_edge),hist.GetBinLowEdge(left_edge)) #special case of 68 existing in one bin
            else:
                return (float(right_val - left_val + 1) / 2., left_val, right_val)

        move_left = hist.GetBinContent(left_edge-1)
        if left_edge == 1:
            move_left = 0
            
        move_right = hist.GetBinContent(right_edge+1)
        
        if left_edge == 1:
            width_right += 1
        elif move_left > move_right:
            width_left += 1
        else:
            width_right += 1

def get_xsec_t5gg(xsec_file,mgl):
    mgl = str(mgl)

    file = open(xsec_file,"r")
    lines = file.readlines()
    lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

    #scan the lines
    for line in lines:
        if mgl in line:
            #check the ordering is correct
            split = line.split()

            xsec = float(split[1])
            xsec_plus = (float(split[2])) * xsec
            xsec_minus = (float(split[2])) * xsec

            return (xsec, xsec_plus, xsec_minus)

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
                lo_xsec = float(split[12])
                lo_xsec_plus = float(split[14])
                lo_xsec_minus = float(split[16])

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


    (xsec,xsec_plus,xsec_minus) = (0,0,0)

    if options.t5gg:
        (xsec,xsec_plus,xsec_minus) = get_xsec_t5gg(xsec_file, msq)
    else:
        (xsec,xsec_plus,xsec_minus) = get_xsec(xsec_file, msq, mgl)

    signal = rt.TFile(signal_file)
    signal_tree = signal.Get("HggOutput")
    
    eff = signal_tree.GetEntries("PFMR > %f & PFR^2 > %f" % (mr_min, rsq_min)) / 10000.

    #calculate the number of signal events expected
    n_events = int(eff * (xsec * 1000) * 19.7 * mu)
    n_events_plus = int(eff * ((xsec+xsec_plus) * 1000) * 19.7 * mu)    
    n_events_minus = int(eff * ((xsec-xsec_minus) * 1000) * 19.7 * mu)

    print "EFFICIENCY FOR SIGNAL INJECTION", eff
    print "XSEC FOR SIGNAL INJECTION", xsec

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



#IF WE ARE GENREATING A MIX OPEN THAT FILE FIRST AND MAKE IT
do_mix =  options.mix_file != "no_file"
tree_signal = None
if do_mix:
    (eff,ns_events, ns_events_plus, ns_events_minus) = build_mixed_model(options.mix_file, options.xsec_file, options.mu, mr_min, r1_cut, options.out_mix_file)
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

name_sample_string = "#gamma#gamma"
if options.samp == 1:
    name_sample_string = "ff"
if options.mix_file != "no_file":
    name_sample_string = "ff + signal"
    
#study_sample_string = "Razor #gamma#gamma + #geq 1 jet: "
#if options.samp == 1:
#    study_sample_string+="Low R^{2} Fit"
#if options.mix_file != "no_file":
#    study_sample_string+=""

#HIGH AND LOW RSQ DATA + SIGNAL HISTOGRAMS
hist_signal = rt.TH1F("hist_signal","Signal",len(bins)-1,array.array("d",bins))
hist_low_rsq_data = rt.TH1F("hist_low_data","Low R^{2} Data (%s)" % name_sample_string, len(bins)-1,array.array("d",bins))
hist_high_rsq_data = rt.TH1F("hist_high_data","High R^{2} Data (%s)" % name_sample_string, len(bins)-1,array.array("d",bins))
hist_high_rsq_signal = rt.TH1F("hist_high_signal","High R^{2} Signal", len(bins)-1,array.array("d",bins))

#JES HISTOGRAMS
hist_high_rsq_jes_up = rt.TH1F("hist_high_data_up","High R^{2} Data (%s) JES UP" % name_sample_string, len(bins)-1,array.array("d",bins))
hist_high_rsq_jes_down = rt.TH1F("hist_high_data_down","High R^{2} Data (%s) JES DOWN" % name_sample_string, len(bins)-1,array.array("d",bins))
hist_low_rsq_jes_up = rt.TH1F("hist_low_data_up","High R^{2} Data (%s) JES UP" % name_sample_string, len(bins)-1,array.array("d",bins))
hist_low_rsq_jes_down = rt.TH1F("hist_low_data_down","High R^{2} Data (%s) JES DOWN" % name_sample_string, len(bins)-1,array.array("d",bins))


low_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && PFR^2 < %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r0_cut, r1_cut,options.samp)
high_cut_data = "(PFMR/1000.)>%f && PFR^2 > %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r1_cut,options.samp)

#JES CUTS
low_cut_data_up = "(PFMR_UP/1000.)>%f && PFR_UP^2 > %f && PFR_UP^2 < %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min,r0_cut, r1_cut,options.samp)
high_cut_data_up = "(PFMR_UP/1000.)>%f && PFR_UP^2 > %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min, r1_cut, options.samp)

low_cut_data_down = "(PFMR_DOWN/1000.)>%f && PFR_DOWN^2 > %f && PFR_DOWN^2 < %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min, r0_cut, r1_cut,options.samp)
high_cut_data_down = "(PFMR_DOWN/1000.)>%f && PFR_DOWN^2 > %f && iSamp==%i && PhotonPFCiC.sieie[0] > .0001 && PhotonPFCiC.sieie[1] > .0001" % (mr_min, r1_cut, options.samp)

low_cut_signal = "(PFMR)>%f && PFR^2 > %f && PFR^2 < %f " % (mr_min, r0_cut, r1_cut)
high_cut_signal = "(PFMR)>%f && PFR^2 > %f" % (mr_min, r1_cut)

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
            hist_high = rt.TH1D(name_high,name_high, 5000,0,10000)
            hist_low = rt.TH1D(name_low,name_low, 50000,0,100000)
        else:
            hist_high = rt.TH1D(name_high,name_high, 100,0,100)
            hist_low = rt.TH1D(name_low,name_low, 10000,0,10000)
            
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
#        if n_toy > 95: continue
        
        #derive the high rsq region parameters
        b_high_toy = b_toy #* pow(r1,1/n_toy)
        n_high_toy = n_toy
        M0_high_toy = M0_toy

        if options.debug: print "b", b_high_toy, "n",n_high_toy
        int_low = int_high = 0

        #find the fit normalization
        try:
            if options.debug: print "\ndoing low integral.."
            int_low = integral_fit(mr_min_low, mr_cutoff, b_toy, M0_toy, n_toy)
            if options.debug: print "\ndoing high integral.."
            int_high = integral_fit(mr_min, mr_cutoff, b_high_toy, M0_high_toy, n_high_toy)
        except OverflowError:
            print "OVERFLOW ERROR n: %f b:%f mr0: %f: " % (n_toy, b_toy, M0_toy)
            continue

        if int_low ==0 or int_high == 0:continue

        fit_toy_summary.add_toy_dist(M0_toy, n_toy, b_toy, Ntot_toy)

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

    return (low_list, high_list, hists_low, hists_high)

tree.Draw("PFMR/1000>>hist_low_data",low_cut_data)
tree.Draw("PFMR/1000>>hist_high_data",high_cut_data)

#JES HISTOGRAMS
tree.Draw("PFMR_UP/1000>>hist_high_data_up", high_cut_data_up)
tree.Draw("PFMR_DOWN/1000>>hist_high_data_down", high_cut_data_down)
tree.Draw("PFMR_UP/1000>>hist_low_data_up", low_cut_data_up)
tree.Draw("PFMR_DOWN/1000>>hist_low_data_down", low_cut_data_down)


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

#fit the histograms with gaussians
for ii in hists_low:    
    print "\n",ii
    ii.Fit("gaus","0q")    
for ii in hists_high:
    ii.Fit("gaus","q0")

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

    high_data = float(hist_high_rsq_data.GetBinContent(ii+1))
    low_data = float(hist_low_rsq_data.GetBinContent(ii+1))
    
    if high_data != 0:
        ratio = (low_data/high_data)
        print "RATIO VALUE:", ratio
    else:
        ratio = -1
        
    #fill the histograms
    hist_low_rsq.SetBinContent(ii+1,low)    
    hist_high_rsq.SetBinContent(ii+1,high)
    hist_ratio.SetBinContent(ii+1,ratio)

    fit_low_hist = hists_low[ii]
    fit_high_hist = hists_high[ii]
    
    low_func = fit_low_hist.GetFunction("gaus")
    high_func = fit_high_hist.GetFunction("gaus")

    low_error = high_error = 1
    
    if fit_low_hist.Integral() > 0:
        low_error = get_smallest68_hist(low,fit_low_hist)[0]
    if fit_high_hist.Integral() > 0:
        high_error = get_smallest68_hist(high,fit_high_hist)[0]        

    #compute the ratio error

    ratio_error = 0
    if low != 0:
        ratio_error = math.sqrt((1/low)*(1/low)*low_error*low_error + (high/low)*(high/low)*high_error*high_error)
    hist_ratio.SetBinError(ii+1, ratio_error)
        
    print "\t\t low_error: %f high_error: %f" % (low_error,high_error)

#    hist_low_rsq.SetBinError(ii+1, low_error)
    #hist_high_rsq.SetBinError(ii+1, high_error)
                     
#clone the low rsq
hist_low_pred = hist_low_rsq.Clone("hist_low_pred")
hist_high_pred = hist_high_rsq.Clone("hist_high_pred")
hist_difference = hist_high_rsq_data.Clone("hist_difference")

#there are no errors on the number of observed data
for ii in range(len(bins)-1):
    hist_difference.SetBinError(ii+1,0)
    #hist_high_rsq_data.SetBinError(ii+1,0)
    
#hist_difference.Sumw2()

#calculate the prediction
hist_difference.Add(hist_high_pred,-1.)
        

hist_ratio.SetFillColor(rt.kOrange+1)
hist_ratio.SetLineColor(rt.kOrange+1)
#hist_ratio.SetFillStyle(3001)
hist_ratio.GetXaxis().SetTitle("M_{R} [TeV]")
hist_ratio.GetYaxis().SetTitle("High R^{2} / Low R^{2}")

hist_ratio.GetXaxis().SetTitleSize(.07)
hist_ratio.GetYaxis().SetTitleSize(.07)

hist_high_pred.GetXaxis().SetTitle("M_{R} [TeV]")
hist_high_pred.GetYaxis().SetTitle("N Events")
hist_high_pred.GetXaxis().SetTitleSize(.07)
hist_high_pred.GetYaxis().SetLabelSize(.13)

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

#END BUILDING HIST DATA
###############################################################
#MAKING CANVASES
################################################################

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

def get_canvases_and_68win(hists, name):
    bin_dists = []
    for ii in range(len(hists)):
        gaus = hists[ii].GetFunction("gaus")
        
        gaus_mean = 0
        gaus_sigma = 1

        if hists[ii].Integral() > 0 :
            gaus_mean = gaus.GetParameter(1)
            gaus_sigma = gaus.GetParameter(2)

        error = gaus_sigma

        true_val = 0
        exp_val = -1

        if name == "high":
            true_val = hist_high_rsq_data.GetBinContent(ii+1)
            exp_val = hist_high_rsq.GetBinContent(ii+1)
        elif name == "low":
            true_val = hist_low_rsq_data.GetBinContent(ii+1)
            exp_val = hist_low_rsq.GetBinContent(ii+1)
        else:
            print "WARNING :: GRAPH NAME NOT FOUND"
        
        rt.gStyle.SetLabelSize(0.05,"xyz");
        rt.gStyle.SetPadTopMargin(0.1);
        rt.gStyle.SetPadRightMargin(0.10);
        rt.gStyle.SetPadBottomMargin(0.2);
        rt.gStyle.SetPadLeftMargin(0.15);
        
        canvas = rt.TCanvas("%s_canvas_bin%i" % (name,ii),"Bin %i" %ii, 1024,768)
                
        #determine the size of the window
        max_val = min_val = -1 
        if name == "high":
            max_val = max(high_list[ii])
            min_val = min(high_list[ii])    
        elif name == "low":
            max_val = max(low_list[ii])
            min_val = min(low_list[ii])    
        else:
            print "WARNING :: GRAPH NAME NOT FOUND"
            

        p3 = gaus_mean + 3 * error
        m3 = gaus_mean - 3 * error
        
        absmin = min([min_val, m3])
        absmax = max([max_val, max_val, p3])
            
        max_hist = hists[ii].GetMaximum()    
        hists[ii].GetXaxis().SetRangeUser(absmin, absmax)
        
        
        hists[ii].Draw()
        
        line = rt.TLine(true_val, 0, true_val, max_hist)
        line.SetLineWidth(3)
        line.Draw("same")
        
        exp_line = rt.TLine(exp_val, 0, exp_val, max_hist)
        exp_line.SetLineWidth(3)
        exp_line.SetLineStyle(2)
        exp_line.Draw("same")
        
        window68 = get_smallest68_hist(exp_val,hists[ii])
        window_l = window68[1]
        window_r = window68[2]
        
        wleft_line = rt.TLine(window_l, 0, window_l, max_hist)
        wleft_line.SetLineWidth(4)
        #wleft_line.SetLineStyle(2)
        wleft_line.SetLineColor(rt.kRed)
        #if window_l > 0: 
        wleft_line.Draw("same")
            
        wright_line = rt.TLine(window_r, 0, window_r, max_hist)
        wright_line.SetLineWidth(4)
        #wright_line.SetLineStyle(2)
        wright_line.SetLineColor(rt.kRed)
        if window_r > 0:
            wright_line.Draw("same")

        #top of the gaussian
        highest_point = 0
        if gaus_mean != 0: highest_point = gaus.Eval(gaus_mean)
    
        line_mean = rt.TLine(gaus_mean, 0, gaus_mean, highest_point)
        line_mean.SetLineColor(rt.kRed)
        line_mean.SetLineWidth(3)
        
        hists[ii].SetFillColor(rt.kAzure-9)
#        hists[ii].SetFillStyle(3000)
        hists[ii].SetLineWidth(1)
        hists[ii].SetLineColor(rt.kAzure-9)
        
        legend = rt.TLegend(.415,.394,.86,.738)
        legend.AddEntry(hists[ii],"Pseudo-experiments","f")
        #legend.AddEntry(gaus,"Gaussian Fit","l")
        legend.AddEntry(line, "N Observed","l")
        legend.AddEntry(exp_line, "N Expected","l")
        legend.AddEntry(wright_line, "68%  window","l")
            

        legend.SetFillColor(0)
        legend.SetLineColor(0)

        if not options.print_pdf: legend.Draw("same")
        
        delta_text= rt.TPaveText(.125,1,.9,.9,"NDC");

        delta = -9
    
        # if you expect 1 or less events do a one sided window
        # otherwise do a two sided window
        if exp_val > 1:
            delta = (true_val - exp_val) / window68[0]
        else:
            delta = (true_val - exp_val) / (window68[0])
            #delta = float(true_val) / float(smallest_68(hists[ii]))

        delta_text.SetTextFont(42);
        delta_text.SetTextSize(0.04);
        delta_text.SetFillColor(0); 
        delta_text.SetTextAlign(12);

        #make the toy distribution object
        toy_dist = toy_distribution(window_l, window_r, exp_val, true_val, delta)

        #add it to the list of istirubtions
        bin_dists.append(toy_dist)

        
        if exp_val > 1:
            delta_text.AddText("#Delta = 2(N_{data} - N_{exp}) / win_{68} = 2(%2.1f - %2.1f) / %2.1f = %2.2f" % (true_val, exp_val, 2*window68[0] ,delta));
        else:
            delta_text.AddText("#Delta = (N_{data} - N_{exp}) / win_{68} = (%2.1f - %2.1f) / %2.1f = %2.2f" % (true_val, exp_val, window68[0] ,delta));
            #delta_text.AddText("#Delta = N_{obs} /  (.68 range) = %2.2f" % delta)

        delta_text.Draw("same")

        canvas.Write()

    return bin_dists

#build a tree out of the toy parameters
# Make a tree
fit_toy_tree = rt.TTree('fit_toy_tree','tree containing the toy fit parameters')

# Create a struct
rt.gROOT.ProcessLine(\
      "struct MyStruct{\
      Float_t mr0;\
      Float_t n;\
      Float_t b;\
      Float_t ntot;\
      };")

from ROOT import MyStruct

# Create branches in the tree
struct = MyStruct()
fit_toy_tree.Branch('mr0',rt.AddressOf(struct,'mr0'),'mr0/F')
fit_toy_tree.Branch('n',rt.AddressOf(struct,'n'),'n/F')
fit_toy_tree.Branch('b',rt.AddressOf(struct,'b'),'b/F')
fit_toy_tree.Branch('ntot',rt.AddressOf(struct,'ntot'),'ntot/F')

mr0_list = fit_toy_summary.get_mr0()
n_list = fit_toy_summary.get_n()
b_list = fit_toy_summary.get_b()
ntot_list = fit_toy_summary.get_ntot()

# Fill tree
for toy in range(fit_toy_summary.get_ntoys()):
    struct.mr0 = mr0_list[toy]
    struct.n = n_list[toy]
    struct.b = b_list[toy]
    struct.ntot = ntot_list[toy]
    fit_toy_tree.Fill()
    
fit_toy_tree.Write()

#draw the toy canvses and get the distributions
bin_dists_high = get_canvases_and_68win(hists_high, "high")
bin_dists_low = get_canvases_and_68win(hists_low, "low")

#build an object containing a collpased list of the information from the bin distirbutions
error_sum_high = error_summary(bin_dists_high)
error_sum_low = error_summary(bin_dists_low)

#build the pull distribution
hist_pull_high = rt.TH1F("pull_hist_high","pull distirbution of high bins",20,-3,3)
hist_pull_low = rt.TH1F("pull_hist_low","pull distirbution of low bins",20,-3,3)

#get the delta distributions from the error summary
delta_high = error_sum_high.get_delta()
delta_low = error_sum_low.get_delta()
#fill the pull distributions
for ii in delta_high: hist_pull_high.Fill(ii)
for ii in delta_low: hist_pull_low.Fill(ii)
#save pulls
hist_pull_high.Write()
hist_pull_low.Write()

#build the graphs that will go in the canvas
(high_exp_up, high_exp_down, graph68_high_exp, graph68_high_obs, hist_high_err_up, hist_high_err_down) = build68graph(error_sum_high, bins, "high")
#error_sum_high.clear_vals()
(low_exp_up, low_exp_down, graph68_low_exp, graph68_low_obs, hist_low_err_up, hist_low_err_down) =  build68graph(error_sum_low, bins, "low")

high_exp_down.GetXaxis().SetTitle("M_{R} [TeV]")
high_exp_down.GetYaxis().SetTitle("N Events")
high_exp_down.GetYaxis().SetTitleSize(.12)
high_exp_down.GetYaxis().SetLabelSize(.09)
high_exp_down.GetYaxis().SetTitleSize(.10)
high_exp_down.GetYaxis().SetTitleOffset(.65)

low_exp_down.GetXaxis().SetTitle("M_{R} [TeV]")
low_exp_down.GetYaxis().SetTitle("N Events")
low_exp_down.GetYaxis().SetTitleSize(.10)
low_exp_down.GetYaxis().SetTitleOffset(.65)
low_exp_down.GetYaxis().SetLabelSize(.09)

high_exp_up.Write()
high_exp_down.Write()
low_exp_up.Write()
low_exp_down.Write()

hist_low_err_up.Write()
hist_low_err_down.Write()

hist_high_err_up.Write()
hist_high_err_down.Write()

print "Writing graphs"

graph68_high_exp.Write("highgraph")
graph68_low_exp.Write("lowgraph")
graph68_high_obs.Write("highgraph_obs")
graph68_low_obs.Write("lowgraph_obs")
hist_high_rsq_jes_up.Write() 
hist_high_rsq_jes_down.Write() 
hist_low_rsq_jes_up.Write() 
hist_low_rsq_jes_down.Write() 
file.Close()

if do_mix:
    os.system("cp weight_hist.root weight_hist_mix.root")
if not options.draw:
    os.system("root -l -c ~/rootlogon.C ratio_extrap.C")
