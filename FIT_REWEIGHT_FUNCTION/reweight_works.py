import ROOT as rt
import array,os
#bins = [.400, .450, .500,.550,.6,.7,.800, 1.000, 1.500, 3]
bins = [ .6,.7,.8,1, 1.5, 2, 3,4,5]

data_file = rt.TFile("/Users/hardenbr/Documents/2013/RAZOR_DIPHOTON/DATA/2012ABCD_FULL_Jan22.root","READ")

do_mix = False

if do_mix:
    data_file = rt.TFile("mix.root")

tree = data_file.Get("HggOutput")

hist_low_rsq_data = rt.TH1F("hist_low_data","Low Rsq Data",len(bins)-1,array.array("d",bins))
hist_low_rsq_data_reweighted = rt.TH1F("hist_low_data_reweight","Low Rsq Data",len(bins)-1,array.array("d",bins))
hist_high_rsq_data = rt.TH1F("hist_high_data","High Rsq Data",len(bins)-1,array.array("d",bins))
hist_high_pred = rt.TH1F("hist_high_pred","High Rsq Pred",len(bins)-1,array.array("d",bins))

mr_min = bins[0]
mr_cutoff = 10
rsq_cutoff = 10
r0 = .01
r1 = .02

low_cut = "(PFMR/1000.)>%f && PFR^2 > %f && PFR^2 < %f && iSamp==1" % (mr_min,r0, r1)
high_cut = "(PFMR/1000.)>%f && PFR^2 > %f && iSamp==1" % (mr_min,r1)

if do_mix:
    low_cut = "(PFMR)>%f && PFR^2 > %f && PFR^2 < %f " % (mr_min,r0, r1)
    high_cut = "(PFMR)>%f && PFR^2 > %f" % (mr_min,r1)

print low_cut
print high_cut

n_low_rsq = tree.GetEntries(low_cut)
n_high_rsq = tree.GetEntries(high_cut) + 10 

print "nlow: ", n_low_rsq
print "nhigh: ", n_high_rsq

b_low = 10.2 #* 1 /( pow(r0, 1/n))# + pow(r1, 1/n))
R0 = 0
n = 67.82
M0 = -.4529

if do_mix:
    b_low = 10.8
    n = 59.7
    M0 = -.91

b_high = b_low * (pow((r1/r0),1/n))
R0_high = 0
n_high = n
M0_high = M0


def Gamma(a,x):
    return rt.TMath.Gamma(a)*rt.Math.inc_gamma_c(a,x)

def Gfun(m,r,b,n):
    val =  n/pow(b*n,n)*pow(10,100)*Gamma(n,b*n*pow(m*r,1/n))
    print val
    return val

def integral(r, m1, m2, b, M0, n):
    return (Gfun(m1-M0,r-R0,b,n) - Gfun(m2-M0,r-R0,b,n))

def integral_fit(r, m1, m2, b, M0 ,n):
    return (Gfun(m1-M0,1,b,n) - Gfun(m2-M0,1,b,n))*pow(10,200)



#fill the trees



if do_mix:
    tree.Draw("PFMR>>hist_low_data","PFMR>%f&&PFR^2>%f && PFR^2 < %f" % (mr_min, r0, r1))
    tree.Draw("PFMR>>hist_high_data","PFMR>%f&&PFR^2>%f" % (mr_min, r1))
else:
    tree.Draw("PFMR/1000.>>hist_low_data","iSamp==1&&PFMR>400&&PFR^2>.01 && PFR^2 < .02")
    tree.Draw("PFMR/1000.>>hist_high_data","iSamp==1&&PFMR>400&&PFR^2>.02")

file = rt.TFile("weight_hist.root","RECREATE")
hist_low_rsq = rt.TH1F("hist_low","Low Rsq",len(bins)-1,array.array("d",bins))
hist_high_rsq = rt.TH1F("hist_high","High Rsq",len(bins)-1,array.array("d",bins))
hist_ratio = rt.TH1F("hist_ratio","Ratio",len(bins)-1,array.array("d",bins))
hist_difference = rt.TH1F("hist_difference","Difference %",len(bins)-1,array.array("d",bins))

for ii in range(len(bins) - 1):
    m1 = bins[ii]
    m2 = bins[ii+1]

    int_low = integral_fit(1, mr_min, mr_cutoff, b_low, M0, n)
    int_high = integral_fit(r1, mr_min, mr_cutoff, b_high, M0_high, n_high) 
    
    print "INTLOW:", int_low
    print "INTHIGH:",  int_high
               
    low = integral_fit(r0, m1, m2, b_low, M0, n)         
    high = integral_fit(r1, m1, m2, b_high, M0_high, n_high) 

    low *= 1. / int_low * (n_low_rsq)
    high *= 1. / int_high * (n_high_rsq)

    ratio = (high/low)

    hist_low_rsq.SetBinContent(ii+1,low)
    hist_high_rsq.SetBinContent(ii+1,high)
    hist_ratio.SetBinContent(ii+1,ratio)


#clone the low rsq
#hist_ratio = hist_low_rsq.Clone("hist_ratio")
#hist_high_pred = hist_low_rsq_data.Clone("hist_high_pred")

#Set the errors of the ratios
#hist_ratio.Sumw2()
#hist_high_pred.Sumw2()

#hist_high_rsq.Sumw2()
#hist_low_rsq.Sumw2()


#calculate the prediction
#hist_ratio.Divide(hist_high_rsq)
#hist_high_pred.Divide(hist_ratio)

for ii in range(1,len(bins)):
    factor = hist_ratio.GetBinContent(ii)    
    val = hist_low_rsq_data.GetBinContent(ii)
    pred = val*factor
    truth = hist_high_rsq_data.GetBinContent(ii)

    hist_high_pred.SetBinContent(ii,pred)

    if truth > 1:
        hist_difference.SetBinContent(ii,(pred-truth)/truth)
    else:
        continue
    

hist_ratio.SetFillColor(rt.kOrange+1)
hist_ratio.SetLineColor(rt.kOrange+1)
hist_ratio.SetFillStyle(3001)
hist_ratio.GetXaxis().SetTitle("M_{R} [TeV]")
hist_ratio.GetYaxis().SetTitle("Low R^{2} / High R^{2}")

hist_ratio.GetXaxis().SetTitleSize(.06)
hist_ratio.GetYaxis().SetTitleSize(.06)


hist_difference.Write()
hist_low_rsq_data.Write()
hist_high_pred.Write()
hist_high_rsq_data.Write()


hist_low_rsq.Write()
hist_high_rsq.Write()
hist_ratio.Write()
file.Close()
    
os.system("root -l -c draw.C")
