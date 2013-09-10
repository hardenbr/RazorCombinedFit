from  optparse  import OptionParser
import ROOT as rt
import RootTools


parser = OptionParser(usage="usage python JoshToy.py -t n_toys -f file_name")

parser.add_option("-f", "--file", dest="filename",
                  help="fit.root file name to analyze FILE",
                  default="razor_output.root",
                  action="store",type="string")

parser.add_option("-t", "--toys", dest="n_toys",
                  help="Number of Toys to Throw",
                  default="100",
                  action="store",type="int")


(options, args) = parser.parse_args()

parser.print_help()


file = rt.TFile(options.filename)

fr = file.Get("Had/independentFR")

toys = []

for ii in range(options.n_toys):
    list = fr.randomizePars()

    list.Print()

    MR = [list[0].getVal(),list[0].getError()]
    Ntot = [list[1].getVal(),list[1].getError()]
    b = [list[2].getVal(),list[2].getError()]
    n = [list[3].getVal(),list[3].getError()]

    toys.append([MR, Ntot, b, n])
    
    



