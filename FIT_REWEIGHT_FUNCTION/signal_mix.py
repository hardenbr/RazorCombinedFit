import ROOT as rt
import random


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

def build_mixed_model(data_file, signal_file, xsec_file, mu, mr_min, rsq_min, outfile):

    mr_min = mr_min*1000
    
    signal_name = signal_file.split("/")[-1]
    msq = signal_name.split("_")[1]
    mgl = signal_name.split("_")[2]    

    (xsec,xsec_plus,xsec_minus) = get_xsec(xsec_file, msq, mgl)

    signal = rt.TFile(signal_file)
    signal_tree = signal.Get("HggOutput")
    data = rt.TFile(data_file)
    data_tree = data.Get("HggOutput")
    
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
        signal_mr_rsq.append((signal_tree.PFMR,signal_tree.PFR))

    #do the same for data
    data_cut = "iSamp==1 && PFMR > %f && PFR^2 > %f" % (mr_min,rsq_min)
    data_tree.Draw(">>iterlist2", data_cut ,"entrylist")
    itlist_data = rt.gDirectory.Get("iterlist2")

    data_mr_rsq = []

    #make a list of all the values forthe fake
    n_data = data_tree.GetEntries(data_cut) 
    print "total data events...", n_data
    for event in range(n_data):
        if event % 100 == 0: print "scanning datapoint....", event
        entry = itlist_data.Next()
        data_tree.GetEntry(entry)
        data_mr_rsq.append((data_tree.PFMR,data_tree.PFR))

    random.shuffle(signal_mr_rsq)



    #output variables
    rt.gROOT.ProcessLine("struct myStruct{\
    Float_t PFMR;\
    Float_t PFR;\
    Int_t isFake;\
    };")


    outfile = rt.TFile(outfile,"RECREATE")
    s = rt.myStruct()

    #generate the output tree
    outTree = rt.TTree("HggOutput","HggOutput")
    outTree.Branch("PFMR",rt.AddressOf(s,"PFMR"),"PFMR/F")
    outTree.Branch("PFR",rt.AddressOf(s,"PFR"),"PFR/F")
    outTree.Branch("isFake",rt.AddressOf(s,"isFake"),"isFake/I")

    print "Filling result signal..."

    for ii in signal_mr_rsq[:n_events]:
        print ii
        s.PFMR = float(ii[0])/ 1000.
        s.PFR = float(ii[1])
        s.isFake = 0
        outTree.Fill()

    print "Filling result data..."
    for ii in data_mr_rsq:
        s.PFMR = ii[0]
        s.PFR = ii[1]
        s.isFake = 1

        if s.PFMR > mr_min and s.PFR*s.PFR > rsq_min:
            s.PFMR = s.PFMR / 1000.
            outTree.Fill()            

    outTree.Write()
    outfile.Close()
    
    return (n_events, n_events_plus, n_events_minus)
    
        
xsec = "Spectra_gsq_B_8TeV.xsec"
signal = "/afs/cern.ch/work/h/hardenbr/2014/RAZOR_DIPHOTON/MC/SIGNAL/selected_1400_1820_375.root"
dataset = "/afs/cern.ch/work/h/hardenbr/2014/RAZOR_DIPHOTON/DATA/2012ABCD_selected_FULL.root"
mu = 1
(n_events, n_events_plus, n_events_minus) = build_mixed_model(dataset, signal, xsec, mu, .6, .01,"test.root")

print "nevents: %f nevents_plus: %f nevents_minus %f" % (n_events, n_events_plus, n_events_minus)

