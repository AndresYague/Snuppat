import sys
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 {} <file1> [file2 ...]".format(sys.argv[0]))
        return 1
    
    # Get carbon and oxygen values
    c12_file = []
    for archiv in sys.argv[1:]:
        c12 = []
        with open(archiv, "r") as fread:
            # Read header
            header = fread.readline().split(",")
            #iC12 = header.index(" C12")
            iC12 = 7
            
            # Now get values
            for line in fread:
                lnlst = line.split()
                c12.append(float(lnlst[iC12]))
        
        c12_file.append(c12)

    labels = ["3", "4", "5"]
    styles = ["bo-", "g^-", "ys-"]
    for ii in range(len(sys.argv[1:])):
        c12 = c12_file[ii]
        plt.plot(range(1, len(c12) + 1), c12, styles[ii], lw = 2,
                markersize = 10, label = labels[ii] + " M$_\odot$")

    plt.legend(loc = 0, numpoints = 1, prop = {'size': 12})
    plt.ylabel("$^{12}$C mass fraction", size = 14)
    plt.xlabel("TDU number", size = 14)
    plt.tick_params(right = True, top = True)
    
    # Minor ticks
    plt.minorticks_on()
    plt.tick_params(which = "minor", right = True)
    plt.tick_params(which = "minor", bottom = False)
    
    plt.subplots_adjust(hspace = 0)
    
    # Set integer tick labels
    maxim_tick = max(plt.xticks()[0])
    plt.xticks(np.arange(0, maxim_tick, step = 1))

    plt.savefig('c12MassFractions.pdf')
    plt.show()

if __name__ == "__main__":
    main()
