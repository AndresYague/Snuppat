import sys
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 {} <file1> [file2 ...]".format(sys.argv[0]))
        return 1
    
    # Get carbon and oxygen values
    c12_file = []; o16_file = []
    for archiv in sys.argv[1:]:
        c12 = []; o16 = []
        with open(archiv, "r") as fread:
            # Read header
            header = fread.readline().split(",")
            #iC12 = header.index(" C12")
            #iO16 = header.index("O16")
            iC12 = 7
            iO16 = 8
            
            # Now get values
            for line in fread:
                lnlst = line.split()
                c12.append(float(lnlst[iC12]))
                o16.append(float(lnlst[iO16]))
        
        c12_file.append(c12)
        o16_file.append(o16)

    fig, axArr = plt.subplots(2, sharex = True)
    labels = ["3", "4", "5"]
    styles = ["bo-", "g^-", "ys-"]
    for ii in range(len(sys.argv[1:])):
        c12 = c12_file[ii]
        o16 = o16_file[ii]
        axArr[0].plot(range(1, len(c12) + 1), c12, styles[ii], lw = 2,
                markersize = 10, label = labels[ii] + " M$_\odot$")
        axArr[1].plot(range(1, len(o16) + 1), o16, styles[ii], lw = 2,
                markersize = 10)

    axArr[0].legend(loc = 0, numpoints = 1, prop = {'size': 12})
    axArr[0].set_ylabel("$^{12}$C mass fraction")
    axArr[0].minorticks_on()
    axArr[0].tick_params(right = True, top = True)
    axArr[0].tick_params(which = "minor", right = True, top = True)

    axArr[1].set_ylabel("$^{16}$O mass fraction")
    axArr[1].set_xlabel("TDU number")
    axArr[1].minorticks_on()
    axArr[1].tick_params(right = True, top = True)
    axArr[1].tick_params(which = "minor", right = True, top = True)
    
    fig.subplots_adjust(hspace = 0)
    for ax in axArr:
        ax.xaxis.label.set_fontsize(12)
        ax.yaxis.label.set_fontsize(12)

    fig.canvas.draw()


    #for ii in range(len(axArr)):
        #axLabs = axArr[ii].get_yticklabels()
        #axLabs[0] = axLabs[-1] = ""
        #if ii == 0:
            #axLabs[1] = ""
        #elif ii == 1:
            #axLabs[-2] = ""
        
        #axArr[ii].set_yticklabels(axLabs)

    plt.savefig('c12O16MassFractions.pdf')
    plt.show()

if __name__ == "__main__":
    main()
