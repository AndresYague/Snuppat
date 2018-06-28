'''Program to create an evolution movie for a given element'''

import sys, os, math, subprocess
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

def element_pos(elem, correct):
    # Species file
    spec = os.path.join("..", "..", "data", "species.dat")
    
    with open(spec, "r") as fread:
        ii = 1
        for line in fread:
            lnlst = line.split()
            stri = "{}{}".format(lnlst[1], lnlst[0])
            elm_mass = int(lnlst[0])
            
            # Compare
            if elem == stri:
                break
            
            ii += 1
    
    if elem != stri:
        print "{} is not in the network!".format(elem)
        return (None, None)
    
    return (ii - 1 + correct, elm_mass)

def plotElements(masses, data, listsOfData, temp, bordMass, age, firstime,
                 convRegions, modNum, elem, elms, discont, xLimit, yLimit):
    
    # Plot preprocess
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel("M/M$_\odot$")
    
    # Plot convective regions
    for region in convRegions:
        plt.fill_between(region, [yLimit[1], yLimit[1]], color = "none",
                         hatch = "/", edgecolor = "k")
    
    # Limit for the plot
    xSpan = xLimit[1] - xLimit[0]
    
    # Floating text
    plt.text(xLimit[0] + xSpan/30., yLimit[1]*0.1,
             "Time = {:.3f} kyr".format(age - firstime), size = 8)
    plt.text(xLimit[0] + xSpan/30., yLimit[1]*0.05,
             "Model number = {}".format(modNum), size = 8)
    
    # Plot chemistry
    ax.set_ylabel("Mass fraction")
    ax.set_yscale("log")
    ax.set_ylim(yLimit)
    
    capel = elem[0].upper() + elem[1:]
    lines = [ax.plot(masses, data, label = capel, lw = 2)]
    for ii in range(len(elms)):
        name = elms[ii]
        capel = name[0].upper() + name[1:]
        if name not in discont:
            lines.append(ax.plot(masses, listsOfData[ii],
                                 label = capel, lw = 2))
        else:
            lines.append(ax.plot(masses, listsOfData[ii], "o",
                                 label = capel, markersize = 5))
    
    # Plot mesh
    constLevel = [math.sqrt(yLimit[0]*yLimit[1]) for x in bordMass]
    lines.append(ax.plot(bordMass, constLevel, "r.", label = "Mesh",
                         lw = 2))
    
    # Plot temperature
    ax2 = ax.twinx()
    ax2.set_ylabel("Temperature (K)")
    ax2.set_yscale("log")
    ax2.set_ylim([1e5, 1e9])
    lines.append(ax2.plot(masses, temp, "k--", label = "Temperature",
                          lw = 2))
    
    # Legend
    lins = None
    for li in lines:
        if lins is None:
            lins = li
        else:
            lins += li
    
    labs = [l.get_label() for l in lins]
    ax.legend(lins, labs, prop = {"size": 8})
    
    # Xlimit
    plt.xlim(xLimit)
    ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5))
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset = False))
    plt.tight_layout()
    
def main():
    
    # Manage input
    if len(sys.argv) < 3:
        print "Use: python {} <simulation file> <element> ".format(sys.argv[0])
        return 1
    
    arch = sys.argv[1]
    elem = sys.argv[2].lower()
    
    # Read input file
    inpt = []; inptErr = False
    with open("input.in", "r") as fIn:
        for line in fIn:
            lnlst = line.split()
            inpt.append(lnlst[0])
    
    # Check correctness
    if len(inpt) < 8:
        inptErr = True
    else:
        outpt = inpt[0]
        strtAge = inpt[1]
        xLow = inpt[2]; xHigh = inpt[3]
        yLow = inpt[4]; yHigh = inpt[5]
        snapFreq = inpt[6]
        createVideo = inpt[7].lower()
        eraseSnapshots = inpt[8].lower()
        adjustAxes = inpt[9].lower()
        
        try:
            strtAge = float(strtAge)
            xLimit = [float(xLow), float(xHigh)]
            yLimit = [float(yLow), float(yHigh)]
            snapFreq = float(snapFreq)
        except ValueError:
            inptErr = True
        except:
            raise
        
        if createVideo not in ["y", "n"]:
            inptErr = True
        
        if eraseSnapshots not in ["y", "n"]:
            inptErr = True
        
        if adjustAxes not in ["y", "n"]:
            inptErr = True
    
    if inptErr:
        print "Wrongly formatted input file, please check Readme.txt"
        return 1
    
    createVideo = True if createVideo == "y" else False
    eraseSnapshots = True if eraseSnapshots == "y" else False
    
    # Set fontsize
    plt.rcParams.update({"font.size": 14})
    
    # Position corrector
    correct = 4
    
    # Get element position
    posit, elm_mass = element_pos(elem, correct)
    
    if posit is None:
        return 1
    
    # List of elements
    elms = ["he4", "n14", "c12", "c13", "ne22", "p1"]
    discont = []
    
    # Graphic elements
    positions = list(); mass_el = list()
    for el in elms:
        pos, mass = element_pos(el, correct)
        positions.append(pos)
        mass_el.append(mass)
    
    # Read data
    with open(arch, "r") as fread:
        nModels = -1 # Plot first model
        kk = 0; firstime = None; brokeNewModel = False
        while True:
            # Look for the next model
            while True:
                # Read next line
                if brokeNewModel:
                    line = line2
                    brokeNewModel = False
                else:
                    line = fread.readline()
                
                if len(line) == 0:
                    break
                
                # Check if new model
                if "Model" in line:
                    # Read next two lines and exit:
                    modNum = int(line.split()[-1])
                    mass = float(fread.readline().split()[-1])
                    age = 10**(float(fread.readline().split()[-1]) - 3)
                    
                    if firstime is None:
                        firstime = age
                    
                    nModels += 1
                    break
            
            if len(line) == 0:
                break
            
            # Restart loop if not in strtAge
            if (age - firstime) < strtAge:
                continue
            
            # Only plot one every snapFreq models
            if nModels%snapFreq != 0:
                continue
            
            # Initialize lists
            listsOfData = list()
            for el in elms:
                listsOfData.append(list())
            
            masses = list()
            bordMass = list()
            temp = list()
            rad = list()
            data = list()
            
            # Now let's go mass by mass
            lnlst1 = fread.readline().split()
            while True:
                line2 = fread.readline()
                lnlst2 = line2.split()
                if len(lnlst2) == 0:
                    break
                
                # Check if in new model
                if "Model" in line2:
                    brokeNewModel = True
                    break
                
                # Calculate and append values
                # Mass
                masses.append((float(lnlst1[0]) + float(lnlst2[0]))*0.5*mass)
                
                bordMass.append(float(lnlst1[0])*mass)
                
                # Temperature
                valtemp = (float(lnlst2[1]) + float(lnlst1[1]))*5e8
                temp.append(valtemp)
                
                # Nabla radiative - nabla adiabatic
                valRad = float(lnlst2[3]) + float(lnlst1[3])
                rad.append(valRad)
                
                # Elements
                for ii in range(len(elms)):
                    pos = positions[ii]
                    val = (float(lnlst2[pos]) +
                           float(lnlst1[pos]))*0.5*mass_el[ii]
                    listsOfData[ii].append(val)
                
                # Changing element
                val = (float(lnlst2[posit]) + float(lnlst1[posit]))*0.5*elm_mass
                data.append(val)
                
                lnlst1 = lnlst2
            
            # Get convective regions
            convRegions = list()
            inConvec = False; xConvLow = None; xConvHi = None
            for i in range(len(masses)):
                # Look for first and last convective
                if rad[i] > 0 and not inConvec:
                    # If more than one convective in a row
                    if rad[i + 1] > 0:
                        xConvLow = masses[i]
                        inConvec = True
                elif rad[i] <= 0 and inConvec:
                    # If more than one radiative in a row
                    if i < len(masses) - 1 and rad[i + 1] <= 0:
                        xConvHi = masses[i - 1]
                        inConvec = False
                        
                        # Add region if its wider than 1e-3 MSun
                        if (xConvHi - xConvLow) > 1e-3:
                            convRegions.append([xConvLow, xConvHi])
            
            # Add last one
            if inConvec:
                convRegions.append([xConvLow, masses[-1]])
            
            # Interactive axes
            plt.ion()
            while adjustAxes == "y":
                print "Creating plot for axes adjusting"
                print
                
                # Create plot
                plotElements(masses, data, listsOfData, temp, bordMass, age,
                             firstime, convRegions, modNum, elem, elms, discont,
                             xLimit, yLimit)
                plt.show()
                
                # Ask for adjustment
                adjustAxes = None
                while adjustAxes not in ["y", "n"]:
                    print "Do you want to adjust the axes? (y/n)"
                    adjustAxes = raw_input()
                
                # Read adjustment and apply it
                if adjustAxes == "y":
                    print "Write the axis name and the limits"
                    print "Example: x 0.53 0.61"
                    
                    axisLims = raw_input().split()
                    
                    if axisLims[0] == "x":
                        xLimit = [float(axisLims[1]), float(axisLims[2])]
                    elif axisLims[0] == "y":
                        yLimit = [float(axisLims[1]), float(axisLims[2])]
                    
                else:
                    storeLimits = None
                    while storeLimits not in ["y", "n"]:
                        print "Do you wish to save these limits? (y/n)"
                        storeLimits = raw_input()
                    
                    # Compose input file. Set interactive adjusting to n
                    if storeLimits == "y":
                        with open("input.in", "r") as fIn:
                            with open("input.tmp", "w") as fOut:
                                for line in fIn:
                                    if "Low x" in line:
                                        # Write a pretty input
                                        s = str(xLimit[0])
                                        fOut.write("{}".format(s))
                                        if len(s) < 8:
                                            for i in range(7 - len(s)):
                                                fOut.write(" ")
                                        
                                        fOut.write(" # Low x-axis bound\n")
                                    elif "High x" in line:
                                        # Write a pretty input
                                        s = str(xLimit[1])
                                        fOut.write("{}".format(s))
                                        if len(s) < 8:
                                            for i in range(7 - len(s)):
                                                fOut.write(" ")
                                        
                                        fOut.write(" # High x-axis bound\n")
                                    elif "Low y" in line:
                                        # Write a pretty input
                                        s = str(yLimit[0])
                                        fOut.write("{}".format(s))
                                        if len(s) < 8:
                                            for i in range(7 - len(s)):
                                                fOut.write(" ")
                                        
                                        fOut.write(" # Low y-ayis bound\n")
                                    elif "High y" in line:
                                        # Write a pretty input
                                        s = str(yLimit[1])
                                        fOut.write("{}".format(s))
                                        if len(s) < 8:
                                            for i in range(7 - len(s)):
                                                fOut.write(" ")
                                        
                                        fOut.write(" # High y-ayis bound\n")
                                    elif "Interactively adjust" in line:
                                        fOut.write("n       # Interactively")
                                        fOut.write(" adjust axes (y/n)\n")
                                    else:
                                        fOut.write(line)
                        
                        os.rename("input.tmp", "input.in")
                    
                    print
                    print "Creating the rest of the plots"
                
                print
                plt.close()
            
            # Stop interactive mode for plot
            plt.ioff()
            
            # Create plot
            plotElements(masses, data, listsOfData, temp, bordMass, age,
                         firstime, convRegions, modNum, elem, elms, discont,
                         xLimit, yLimit)
            
            # Save image
            plt.savefig("{}_{:05d}.png".format(elem, kk))
            plt.close()
            kk += 1
    
    # Create video
    if createVideo:
        subprocess.call(["ffmpeg", "-framerate", "10", "-i",
                    "{}_%d.png".format(elem), "{}{}.mp4".format(elem, outpt)])
    
    # Erase images
    if eraseSnapshots:
        subprocess.call("rm {}*.png".format(elem), shell = True)

if __name__ == "__main__":
    main()
