import sys, math, os
import matplotlib.pyplot as plt

def main():
    # Check that there's at least one argument
    if len(sys.argv) < 2:
        print("Usage python {}".format(sys.argv[0]), end = " ")
        print("<file1> [<file2> ...]")
        return 1
    
    # Automatically detect if decayed
    if "decayed" in sys.argv[1]:
        plotDecayed = True
    else:
        plotDecayed = False
    
    # Read input file
    fil = "finalGraph.in"
    if os.path.isfile(fil):
        with open(fil, "r") as fread:
            lstyles = fread.readline().strip().split()
            
            labs = []
            for line in fread:
                labs.append(line.strip())
    
    lowZ = 34 # Lowest z value to represent
    highZ = 56 # Highest z value to represent
    
    # Read "species.dat" and store all the values in lists
    species = "../../data/species.dat"
    atomicNum = []; atomicMass = []; namesZ = {}
    with open(species, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            # Correct special names
            if lnlst[1] == "d" or lnlst[2] == "0":
                lnlst[1] = "h"
            
            # Now relate positions with atomic numbers, atomic masses, and names
            zNum = int(lnlst[0]) - int(lnlst[2])
            
            atomicNum.append(zNum)
            atomicMass.append(int(lnlst[0]))
            namesZ[lnlst[1]] = zNum
    
    # Read all initial solar values
    solar = "../../data/solarVals.dat"
    solarValues = {}
    with open(solar, "r") as fread:
        for line in fread:
            lnlst = line.split()
            isotName = lnlst[0] + lnlst[2]
            
            # Add mass fraction value per atomic number
            key = namesZ[lnlst[0]]; val = float(lnlst[1])*float(lnlst[2])
            solarValues[key] = solarValues.get(key, 0) + val
    
    # Go file by file
    numDens = []
    for archivo in sys.argv[1:]:
        
        # Open file for reading
        dens = []
        fread = open(archivo, "r")
        
        # Each line has mass, temperature, rho, radiat
        # and elements in number fraction
        newline = None
        for line in fread:
            if "#" in line:
                continue
            
            lnlst = line.split()
            if len(lnlst) == 0:
                if plotDecayed:
                    break
                else:
                    continue
            
            if not plotDecayed:
                # Surface (newline[0] is the mass)
                prevline = newline
                newline = [float(x) for x in lnlst]
                if newline[0] > 0.85:
                    break
            
            if plotDecayed:
                dens.append(float(lnlst[1]))
        
        # Close file
        fread.close()
        
        # Calculate values of interest
        if plotDecayed:
            numDens.append(dens)
        else:
            numDens.append([(x + y)*0.5 for (x, y) in
                            zip(prevline[4:], newline[4:])])
    
    # Calculate now the agb values and print the surface mass fractions per
    # each isotope
    print("# Surface number fraction values")
    agbValues = []
    for ii in range(len(numDens)):
        dic = {}
        dens = numDens[ii]
        
        # Print the model name
        print("# {}".format(sys.argv[ii + 1]))
     
        # Add the values for each element
        for jj in range(len(atomicNum)):
            key = atomicNum[jj]
            dic[key] = dic.get(key, 0) + dens[jj]*atomicMass[jj]
            
            # Print the number fraction
            print(dens[jj])
        
        agbValues.append(dic)
        print("")
    
    # Now identify iron:
    ironNumber = namesZ["fe"]
    
    # Now divide every element by iron
    for dens in agbValues:
        ironDens = dens[ironNumber]
        for key in dens:
            dens[key] /= ironDens
    
    # Solar as well
    ironDens = solarValues[ironNumber]
    for key in solarValues:
        solarValues[key] /= ironDens
    
    # Now create the final values
    finalValues = []
    zList = [x for x in solarValues.keys()]
    zList.sort(); hFe = []
    for ii in range(len(agbValues)):
        dens = agbValues[ii]
        thisDens = []
        for key in zList:
            if key == 1:
                hFe.append(math.log10(dens[key]/solarValues[key]))
            
            if key < lowZ or key > highZ:
                continue
            
            val = math.log10(dens[key]/solarValues[key])
            
            # Modify values for [X/H]
            val -= hFe[ii]
            thisDens.append(val)
        
        finalValues.append(thisDens)
    
    # Create xaxis:
    xx = [x for x in zList if x >= lowZ and x <= highZ]
    
    # Print final values
    print("# [X/Fe] values")
    for ii in range(len(sys.argv[1:])):
        print("# {}".format(sys.argv[ii + 1]))
        print("")
        
        for jj in range(len(xx)):
            print(xx[jj], finalValues[ii][jj])
        
        print("")
    
    # From zList create contIndx. This list contains a number of
    # tuples with the first and last index of any contiguous sequence
    indx = 1; first = 0
    prevKey = None; contIndx = []
    for key in xx:
        if prevKey is None:
            prevKey = key
            continue
        
        # Check if keys are contiguous
        if key - prevKey > 1:
            contIndx.append((first, indx))
            first = indx
        
        prevKey = key
        indx += 1
    
    # Add last tuple
    contIndx.append((first, indx + 1))
    
    # Begin plot
    figure = plt.figure()
    plt.xlabel("Atomic number Z", size = 14)
    plt.ylabel("[X/H]", size = 14)
    
    # Plot values
    if labs is None:
        labs = sys.argv[1:]
    
    ii = 0
    for dens in finalValues:
        
        # Plot first range
        first, last = contIndx[0]
        if lstyles is None:
            lin, = plt.plot(xx[first:last], dens[first:last],
                    label = labs[ii], lw = 2)
        else:
            lin, = plt.plot(xx[first:last], dens[first:last], lstyles[ii],
                    label = labs[ii], lw = 2)
        
        # Get color and line style
        col, lst = lin.get_color(), lin.get_linestyle()
        colStyle = col + lst
        for elem in contIndx[1:]:
            first, last = elem
            plt.plot(xx[first:last], dens[first:last], colStyle, lw = 2)
        
        ii += 1
    
    # Set floating text
    namAtm = {"Se":34, "Kr":36, "Sr":38, "Zr":40,
            "Mo":42, "Pd":46, "Cd":48, "Sn":50, "Te":52, "Ba":56,
            "Ce":58, "Nd":60, "Sm":62, "Dy":66, "Er":68,
            "Yb":70, "Hf":72, "Os":76, "Hg":80, "Pb":82,
            "Rb":37, "Cs":55, "Xe":54, "Br": 35}
    
    rNamAtm = []
    
    for name in namAtm:
        yVal = 0
        if namAtm[name] < lowZ or namAtm[name] > highZ:
            continue
        
        for ii in range(len(xx)):
            if xx[ii] == namAtm[name]:
                yVal = finalValues[-1][ii]
                break
        
        plt.text(namAtm[name] - 0.5, yVal*1.01, name, size = 14)
        
        # Dots for the plot
        #if name in rNamAtm:
            #plt.plot(namAtm[name], yVal, "ro")
        #else:
            #plt.plot(namAtm[name], yVal, "ko")
    
    # Observations values
    # (NGC3918, NGC7027)
    # Elements: Se, Kr, Rb, Xe
    xxObs = [34, 35, 36, 37, 52, 54]
    yyErrs = [
             [0.10, "-", 0.11, 0.13, "-", 0.11]
            ,[0.17, 0.08, 0.09, 0.13, 0.12, 0.13]
             ]
    PNe = [
          [0.14, "-", 0.65, 0.25, "-", 0.38]
         ,[0.26, 0.13, 0.84, 0.70, 0.56, 0.66]
          ]
    
    gray = (0.75, 0.75, 0.75)
    mrk = ["^", "s"]
    fillstyle = ["none", "full"]
    for star_ii in range(len(PNe)):
        xxHere = []; yyHere = []; errHere = []
        for ii in range(len(xxObs)):
            if PNe[star_ii][ii] == "-":
                continue
            else:
                # Fill a region for each error barr
                plt.fill_between([xxObs[ii] - 0.5, xxObs[ii] + 0.5],
                        y1 = PNe[star_ii][ii] - yyErrs[star_ii][ii],
                        y2 = PNe[star_ii][ii] + yyErrs[star_ii][ii],
                        facecolor = gray, edgecolor = "k")
                
                # Now append things for the actual plot
                xxHere.append(xxObs[ii])
                yyHere.append(PNe[star_ii][ii])
                errHere.append(yyErrs[star_ii][ii])
        
        plt.plot(xxHere, yyHere, "k" + mrk[star_ii],
                fillstyle = fillstyle[star_ii], lw = 2, ms = 8)
    
    plt.minorticks_on()
    plt.tick_params(right = True, top = True)
    plt.tick_params(which = "minor", right = True, top = True)
    plt.legend(loc=0, ncol = 2, prop = {'size': 12})
    plt.show()

if __name__ == "__main__":
    main()
