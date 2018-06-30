import sys, math, os
import matplotlib.pyplot as plt

def main():
    # Check that there's at least one argument
    if len(sys.argv) < 2:
        print "Usage python {} <input>".format(sys.argv[0])
        return 1
    
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
    
    # Open file for reading
    fread = open(sys.argv[1], "r")
    
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
    numDens = [(x + y)*0.5 for (x, y) in zip(prevline[4:], newline[4:])]
    
    # Calculate now the agb values
    agbValues = {}
    
    # Add the AGB values
    for ii in range(len(atomicNum)):
        key = atomicNum[ii]
        agbValues[key] = agbValues.get(key, 0) + numDens[ii]
    
    # Now identify iron:
    ironNumber = None
    for ii in range(len(atomicNum)):
        if names[ii].lower() == "fe":
            ironNumber = atomicNum[ii]
            break
    
    # Now divide every element by iron
    ironDens = agbValues[ironNumber]
    for key in agbValues:
        agbValues[key] /= ironDens
    
    # Solar as well
    ironDens = solarValues[ironNumber]
    for key in solarValues:
        solarValues[key] /= ironDens
    
    # Now create the final values
    zList = solarValues.keys()
    zList.sort()
    finalValues = []
    for key in zList:
        if key < lowZ:
            continue
        
        val = math.log10(agbValues[key]/solarValues[key])
        finalValues.append(val)
    
    # Create xaxis:
    xx = [x for x in zList if x >= lowZ]
    
    # Set floating text
    namAtm = {"Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga": 31, "Ge":32, "As":33,
            "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40,
            "Nb":41, "Mo":42, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48,
            "In":49, "Sn":50, "Sb":51, "Te":52, "I":53, "Xe":54, "Cs":55,
            "Ba":56, "La":57, "Ce":58, "Pr":59, "Nd":60, "Sm":62, "Eu":63,
            "Gd":64, "Tb":65, "Dy":66, "Ho": 67, "Er":68, "Tm": 69, "Yb":70,
            "Lu":71, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77,
            "Pt":78, "Au":79, "Hg":80, "Tl":81, "Pb":82, "Bi":83, "Po":84}
    
    zName = dict((val, key) for key, val in namAtm.iteritems())
    
    print "# " + sys.argv[1]
    print
    for jj in range(len(xx)):
        if jj%3 <= 1:
            print "{} & {} & {:.2f} & ".format(zName[xx[jj]], xx[jj], finalValues[jj]),
        else:
            print "{} & {} & {:.2f}\\\\".format(zName[xx[jj]], xx[jj], finalValues[jj])

if __name__ == "__main__":
    main()
