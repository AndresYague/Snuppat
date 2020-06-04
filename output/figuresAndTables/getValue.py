import sys, math, os, string

def printValue(stri, dic):
    if stri == "HsLs":
        valHs = (dic["ba"] + dic["la"] + dic["ce"])/3
        valLs = (dic["sr"] + dic["y"] + dic["zr"])/3
        val = valHs - valLs
    elif stri == "PbHs":
        valHs = (dic["ba"] + dic["la"] + dic["ce"])/3
        val = dic["pb"] - valHs
    elif stri == "YCeAvg":
        val = 0.5*(dic["ce"] + dic["y"])
    elif stri == "CeY":
        val = dic["ce"] - dic["y"]
    else:
        val = dic[stri.lower()]
    
    print("{} {:.2f}".format(stri, val))

def main():
    '''Returns the mass fraction, Y and [X] values of the asked species'''
    
    # Check that there's at least one argument
    if len(sys.argv) < 3:
        print("Usage: python3 {}".format(sys.argv[0]), end = " ")
        print("<file> <species> [species2 ...]")
        return 1
    
    archivo = sys.argv[1]
    valsString = sys.argv[2:]
    valsString = [x.lower() for x in valsString]
    
    # Automatically detect if decayed
    if "decayed" in sys.argv[1]:
        plotDecayed = True
    else:
        plotDecayed = False
    
    # Read "species.dat" and store all the values in lists
    species = "../../data/species.dat"
    atomicNum = []; atomicMass = []
    specNames = []; isotNames = []
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
            specNames.append(lnlst[1])
            isotNames.append(lnlst[1] + lnlst[0])
    
    # Read all initial solar values
    solar = "../../data/solarVals.dat"
    solarValuesMass = {}; solarValuesDens = {}
    with open(solar, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            isotName = lnlst[0] + lnlst[2]
            specName = lnlst[0]
            numDens = float(lnlst[1])
            mass = float(lnlst[2])
            
            # Add mass fraction and number fraction
            solarValuesMass[specName] = solarValuesMass.get(specName, 0) +\
                                        numDens*mass
            solarValuesDens[isotName] = numDens
    
    # Store number abundances
    dens = []
    with open(archivo, "r") as fread:
        
        # Each line has mass, temperature, rho, radiat
        # and elements in number fraction, for non-decayed plot
        
        # For a decayed plot each line has simply the abundance in
        # number fraction
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
            
            if plotDecayed:
                dens.append(float(lnlst[1]))
            else:
                # Surface (newline[0] is the mass)
                prevline = newline
                newline = [float(x) for x in lnlst]
                if newline[0] > 0.85:
                    break
        
        # Calculate dens if not decayed
        if not plotDecayed:
            dens = [(x + y)*0.5 for (x, y) in zip(prevline[4:], newline[4:])]
    
    # Make dictionary with number densities
    numDensStar = dict(zip(isotNames, dens))
    
    # Get mass fractions
    massFracStar = {}
 
    # Add the values for each element
    for jj in range(len(specNames)):
        key = specNames[jj]
        key2 = isotNames[jj]
        massFracStar[key] = massFracStar.get(key, 0) + dens[jj]*atomicMass[jj]
        massFracStar[key2] = dens[jj]*atomicMass[jj]
    
    # Calculate all the [X] abundances
    zList = [x for x in solarValuesMass.keys()]
    zList.sort()
    xLog = {}
    for key in zList:
        if key in solarValuesMass:
            val = math.log10(massFracStar[key]/solarValuesMass[key])
            xLog[key] = val
    
    # Print all values
    print("# Number density:")
    for name in valsString:
        if name[-1] not in string.digits:
            continue
        
        if name not in numDensStar:
            s = "{} not in database".format(name)
        else:
            s = "Y_{} = {:.2E}".format(name, numDensStar[name])
        print(s)
    print()
    
    print("# Solar number density:")
    for name in valsString:
        if name[-1] not in string.digits:
            continue
        
        if name not in solarValuesDens:
            s = "{} not in database".format(name)
        else:
            s = "s_Y_{} = {:.2E}".format(name, solarValuesDens[name])
        print(s)
    print()
    
    print("# Mass fraction:")
    for name in valsString:
        if name not in massFracStar:
            s = "{} not in database".format(name)
        else:
            s = "X_{} = {:.2E}".format(name, massFracStar[name])
        print(s)
    print()
    
    print("# Solar mass fraction:")
    for name in valsString:
        if name not in solarValuesMass:
            s = "{} not in database".format(name)
        else:
            s = "s_X_{} = {:.2E}".format(name, solarValuesMass[name])
        print(s)
    print()
    
    print("# [X]:")
    for name in valsString:
        if name[-1] in string.digits:
            continue
        
        if name not in xLog:
            s = "{} not in database".format(name)
        else:
            s = "[{}] = {:.2E}".format(name, xLog[name])
        print(s)
    print()

if __name__ == "__main__":
    main()
