import sys, math, os
import matplotlib.pyplot as plt

def main():
    # Check that there's at least one argument
    if len(sys.argv) < 3:
        print "Usage python {} <title>".format(sys.argv[0]),
        print "<file1> [<file2> ...]"
        return 1
    
    # Read input file
    if os.path.isfile("inputFinalGraph.in"):
        with open("inputFinalGraph.in", "r") as fread:
            lstyles = fread.readline().strip().split()
            
            labs = []
            for line in fread:
                labs.append(line.strip())
    
    title = sys.argv[1]
    lowZ = 27 # Lowest z value to represent
    
    # Read all initial solar values
    solar = "../../data/solarVals.dat"
    solarValues = {}
    with open(solar, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            # Add value per atomic number
            key = int(lnlst[2]); val = float(lnlst[1])
            solarValues[key] = val
    
    # Read "species.dat" and store all the values in lists
    species = "../../data/species.dat"
    names = []; atomicNum = []; agbSpecies = []
    with open(species, "r") as fread:
        for line in fread:
            lnlst = line.split()
            names.append(lnlst[1])
            
            # Put specific isotope in list (name + mass)
            agbSpecies.append(lnlst[1] + lnlst[0])
            
            # Now relate positions with atomic numbers
            atomicNum.append(int(lnlst[0]) - int(lnlst[2]))
    
    # Go file by file
    numDens = []
    for archivo in sys.argv[2:]:
        
        # Open file for reading
        fread = open(archivo, "r")
        
        # Each line has mass, temperature, rho, radiat
        # and elements in number fraction
        newline = None
        for line in fread:
            if "#" in line:
                continue
            
            lnlst = line.split()
            if len(lnlst) == 0:
                continue
            
            # Surface (newline[0] is the mass)
            prevline = newline
            newline = [float(x) for x in lnlst]
            if newline[0] > 0.85:
                break
        
        # Close file
        fread.close()
        
        # Calculate values of interest
        numDens.append([(x + y)*0.5 for (x, y) in
                        zip(prevline[4:], newline[4:])])
    
    # Calculate now the agb values
    agbValues = []
    for dens in numDens:
        dic = {}
        
        # First add the values we are not tracking in the AGB
        for key in solarValues:
            if key in atomicNum:
                dic[key] = solarValues[key]
        
        # Now add the rest
        for ii in range(len(atomicNum)):
            key = atomicNum[ii]
            dic[key] = dic.get(key, 0) + dens[ii]
        
        agbValues.append(dic)
    
    # Now identify iron:
    ironNumber = None
    for ii in range(len(atomicNum)):
        if names[ii].lower() == "fe":
            ironNumber = atomicNum[ii]
            break
    
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
    zList = solarValues.keys()
    zList.sort()
    for dens in agbValues:
        thisDens = []
        for key in zList:
            if key < lowZ:
                continue
            
            val = math.log10(dens[key]/solarValues[key])
            thisDens.append(val)
        
        finalValues.append(thisDens)
    
    # Create xaxis:
    xx = [x for x in zList if x >= lowZ]
    
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
    plt.title(title, size = 14)
    plt.xlabel("Atomic number Z", size = 14)
    plt.ylabel("[X/Fe]", size = 14)
    
    # Plot values
    if labs is None:
        labs = sys.argv[3:]
    
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
    namAtm = {"Co":27, "Ge":32, "Se":34, "Kr":36, "Sr":38, "Zr":40,
            "Mo":42, "Pd":46, "Cd":48, "Sn":50, "Te":52, "Ba":56,
            "Ce":58, "Nd":60, "Sm":62, "Gd":64, "Dy":66, "Er":68,
            "Yb":70, "Hf":72, "W":74, "Os":76, "Hg":80, "Pb":82}
    
    for name in namAtm:
        yVal = 0
        for ii in range(len(xx)):
            if xx[ii] == namAtm[name]:
                yVal = finalValues[-1][ii]
                break
        
        plt.text(namAtm[name] - 0.5, yVal*1.01, name, size = 14)
        plt.plot(namAtm[name], yVal, "ko")
    
    plt.legend(loc=0, ncol = 2)
    plt.show()

if __name__ == "__main__":
    main()
