import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main():
    if len(sys.argv) < 2:
        print "Usage: python {} profile1 <profile2 ...>".format(sys.argv[0])
        return 1
    
    subplotNum = len(sys.argv[1:])
    fig = plt.figure()
    
    colors = ["b-", "y-"]
    labels = ["$\omega = 0.08$ Bounded", "$\omega = 0.12$ Boundless"]
    
    # Plot all the profiles
    ii = 0
    for arch in sys.argv[1:]:
        # Read event
        fread = open(arch, "r")
        
        # Add subplot
        ii += 1; isFirstPlot = True
        ax = fig.add_subplot(subplotNum, 1, ii)
        
        lastMass = 0
        maxProfileHeight = 0
        mass = []; profile = []
        maxMass = []; maxProfile = []
        while True:
            line = fread.readline()
            label = arch
            
            if "#" not in line and len(line) > 0:
                lnlst = line.split()
                mass.append(float(lnlst[0]))
                profile.append(float(lnlst[1]))
            
            if "# Profile" in line:
                if len(mass) > 0:
                    # Store the largest profile
                    maxHeight = max(profile)
                    if maxHeight > maxProfileHeight and profile[-1] <= 0:
                        maxProfileHeight = maxHeight
                        maxMass = mass
                        maxProfile = profile
                
                mass = []; profile = []
            
            if "# TDU" in line or len(line) == 0:
                # Recover mass and profile
                mass = maxMass
                profile = maxProfile
                
                # Shift mass so profiles are consecutive
                if len(mass) > 0:
                    mass0 = mass[0]
                    mass = [x - mass0 + lastMass for x in mass]
                    
                    # Plot
                    if isFirstPlot:
                        ax.set_ylabel("Mass fraction")
                        ax.plot(mass, profile, colors[ii - 1], lw = 2,
                                label = labels[ii - 1])
                        isFirstPlot = False
                        
                    else:
                        ax.plot(mass, profile, colors[ii - 1], lw = 2)
                    
                    # Modify lastMass
                    lastMass += mass[-1] - mass[0]
                
                # Restart variables
                saveMass = mass
                maxProfileHeight = 0
                maxMass = []; mass = []
                maxProfile = []; profile = []
                
                # Exit if last line
                if len(line) == 0:
                    break
        
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1e-2))
        ax.xaxis.set_major_locator(ticker.MaxNLocator(prune = "both"))
        ax.set_xlim([0, max(saveMass)*1.1])
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.legend(loc = 0)
        
        # Clean xaxis
        xticks = list(ax.get_xticks())
        newXticks = []
        for elem in xticks:
            if len(newXticks) == 0 or abs(elem - newXticks[-1]) >= 5.001e-5:
                newXticks.append(elem)
            
        
        ax.set_xticks(newXticks)
        
        # Prune the yaxis
        yLim = ax.get_ylim()
        yticks = list(ax.get_yticks())
        yticks = [x for x in yticks if x >= yLim[0] and x <= yLim[1]*1.1]
        yticks.pop(0); yticks.pop()
        ax.set_yticks(yticks)
        
        fread.close()
    
    ax.set_xlabel("Mass (M$_\odot$)")
    plt.show()

if __name__ == "__main__":
    main()
