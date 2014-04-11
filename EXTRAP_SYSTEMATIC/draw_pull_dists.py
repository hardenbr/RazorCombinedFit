from optparse  import OptionParser
import ROOT as rt
import RootTools
import array, os
import random, math
import rootlogon

rootlogon.style()

parser = OptionParser()

parser.add_option("-f", "--file", dest="filename",
                  help="fit.root file name to analyze FILE",
                  default="razor_output.root",
                  action="store",type="string")

parser.add_option("-t", "--toys", dest="toys",
                  help="weight_hist.root file containing toy tree",
                  default="weight_hist.root",
                  action="store",type="string")

(options, args) = parser.parse_args()

parser.print_help()


file = rt.TFile(options.toys)
fitfile = rt.TFile(options.filename)

tree = file.Get("fit_toy_tree")
fr = fitfile.Get("Had/independentFR")
fr_list = fr.floatParsFinal()

b = fr_list[2].getVal()
n = fr_list[3].getVal()
M0 = fr_list[0].getVal()

b_e = fr_list[2].getError()
n_e = fr_list[3].getError()
M0_e = fr_list[0].getError()

hist_m = rt.TH1F("hm","hm",100,-5,5)
hist_b = rt.TH1F("hb","hb",100,-5,5)
hist_n = rt.TH1F("hn","hn",100,-5,5)

#do the drawing
rt.gStyle.SetOptTitle(0)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(1)

tree.Draw("(mr0 - %f)/%f>>hm" % (M0,M0_e))
hist_m.GetXaxis().SetTitle("(M_{R}^{0,toy} - M_{R,fit}^{0}) / #sigma_{fit,M_{R}^{0}}")
hist_m.Fit("gaus")

c2 = rt.TCanvas()

tree.Draw("(b - %f)/%f>>hb" % (b,b_e))
hist_b.GetXaxis().SetTitle("(b_{toy} - b_{fit}) / #sigma_{fit,b}")
hist_b.Fit("gaus")

c3 = rt.TCanvas()

tree.Draw("(n - %f)/%f>>hn" % (n,n_e))
hist_n.GetXaxis().SetTitle("(n_{toy} - n_{fit}) / #sigma_{fit,n}")
hist_n.Fit("gaus")
raw_input("RAW INPUT")
