import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main():
    if len(sys.argv) < 2:
        print "Usage: python {} profile1 <profile2 ...>".format(sys.argv[0])
        return 1
    
    subplotNum = len(sys.argv[1:])
    fig = plt.figure()
    
    colors = ["b-", "g-", "y-"]
    labels = ["$\omega = 0.10$", "$\omega = 0.12$", "$\omega = 0.14$"]
    
    # Plot all the profiles
    ii = 0; maxSaveMass = 0
    for arch in sys.argv[1:]:
        # Read event
        fread = open(arch, "r")
        
        # Add subplot
        ii += 1; isFirstPlot = True
        ax = fig.add_subplot(subplotNum, 1, ii)
        
        lastMass = 0
        maxProfileHeight = 0
        mass = []; profile = []
        maxProfile = []
        c13Count = 0; c13TotMass = 0
        while True:
            line = fread.readline()
            label = arch
            
            if "#" not in line and len(line) > 0:
                lnlst = line.split()
                if len(lnlst) == 0:
                    continue
                
                mass.append(float(lnlst[0]))
                profile.append(float(lnlst[1]))
                
            elif len(mass) > 0:
                # Count one c13 pocket more
                c13Count += 1
                
                # Add the mass to the maximum
                c13TotMass += max(mass) - min(mass)
                
                # Shift mass so profiles are consecutive
                mass0 = mass[0]
                mass = [x - mass0 + lastMass for x in mass]
                
                # Plot
                if isFirstPlot:
                    ax.plot(mass, profile, colors[ii - 1], lw = 2,
                            label = labels[ii - 1])
                    isFirstPlot = False
                    
                else:
                    ax.plot(mass, profile, colors[ii - 1], lw = 2)
                
                # Modify lastMass
                lastMass += mass[-1] - mass[0]
                
                # Restart variables
                saveMass = mass; mass = []
                maxProfileHeight = 0
                maxProfile = []; profile = []
                
                # Exit if last line
                if len(line) == 0:
                    break
                mass = []; profile = []
        
        if max(saveMass) > maxSaveMass:
            maxSaveMass = max(saveMass)
        
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1e-2))
        ax.xaxis.set_major_locator(ticker.MaxNLocator(prune = "both"))
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.legend()
        
        # Prune the yaxis
        yLim = ax.get_ylim()
        yticks = list(ax.get_yticks())
        yticks = [x for x in yticks if x >= yLim[0] and x <= yLim[1]*1.1]
        yticks.pop(0); yticks.pop()
        ax.set_yticks(yticks)
        
        fread.close()
        
        # Print the average size
        if c13Count > 0:
            avgMass = c13TotMass/c13Count
            print "Average pocket mass in {} = {}".format(arch, avgMass)
    
    # Fix axes
    ii = 0
    for axi in fig.axes:
        axi.set_xlim([0, maxSaveMass*1.1])
        
        # Set central ylabel
        if ii == 1:
            axi.set_ylabel("Mass fraction")
        
        ii += 1
    
    fig.subplots_adjust(hspace = 0)
    ax.set_xlabel("Mass (M$_\odot$)")
    
    plt.show()

if __name__ == "__main__":
    main()
