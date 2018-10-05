import sys, math, numpy
import matplotlib.pyplot as plt

# TODO add [X/Fe] possibility and option

def main():
    '''Get evolution of one element in epsilon or [X/Fe]'''
    
    # Check arguments
    if len(sys.argv) < 4:
        print "Usage python {} <(eps|xfe)> <model>".format(sys.argv[0]),
        print "<elem1> [elem2, elem3, ...]"
        return 1
    
    data = "../../data/species.dat"
    mode = sys.argv[1]
    archivo = sys.argv[2]
    elms = sys.argv[3:]
    
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
    
    # Now go model by model, calculating everything for every element
    with open(archivo, "r") as fread:
        
        # Each line has mass, temperature, rho, radiat
        # and elements in number fraction
        ages = []; evolEps = []; evolXFe = []
        newline = None; foundLine = False
        for line in fread:
            if "#" in line:
                if "Model" in line and newline is not None:
                    # Take all abundances
                    dens = [(x + y)*0.5 for (x, y) in
                            zip(prevline[4:], newline[4:])]
                    
                    epsVals = {}; xFeVals = {}
                    # Add the values for each element
                    for ii in range(len(atomicNum)):
                        key = atomicNum[ii]
                        epsVals[key] = epsVals.get(key, 0) + dens[ii]
                        xFeVals[key] = xFeVals.get(key, 0) +\
                                       dens[ii]*atomicMass[ii]
                    
                    # Now calculate values of interest
                    hydroVal = epsVals[namesZ["h"]]
                    feVal = xFeVals[namesZ["fe"]]
                    sunFeVal = solarValues[namesZ["fe"]]
                    selectedEps = []; selectedFe = []
                    for elem in elms:
                        try:
                            val = epsVals[namesZ[elem]]/hydroVal + 1e-100
                        except KeyError:
                            print "{} is not on the list".format(elem)
                        except:
                            raise
                        
                        val = math.log10(val) + 12
                        selectedEps.append(val)
                        
                        try:
                            val = xFeVals[namesZ[elem]]/feVal + 1e-100
                        except KeyError:
                            print "{} is not on the list".format(elem)
                        except:
                            raise
                        sunVal = solarValues.get(namesZ[elem], 1e-100)/sunFeVal
                        val = math.log10(val) - math.log10(sunVal)
                        selectedFe.append(val)
                    
                    evolEps.append(selectedEps)
                    evolXFe.append(selectedFe)
                    
                    # Reset foundLine
                    foundLine = False
                    
                elif "Age" in line:
                    age = 10**(float(line.split()[-1]) - 3)
                    if len(ages) == 0:
                        ages.append(age)
                    else:
                        ages.append(age - ages[0])
                
                continue
            
            if foundLine:
                continue
            
            lnlst = line.split()
            if len(lnlst) == 0:
                continue
            
            # Surface (newline[0] is the mass)
            prevline = newline
            newline = [float(x) for x in lnlst]
            if newline[0] > 0.85:
                foundLine = True
    
    # Transform age and evol values to something plottable
    ages[0] = 0
    evolEps = numpy.transpose(numpy.array(evolEps))
    evolXFe = numpy.transpose(numpy.array(evolXFe))
    
    # Now plot values
    if mode == "eps":
        for ii in range(len(elms)):
            evEps = evolEps[ii]
            minLen = min(len(ages), len(evEps))
            plt.plot(ages[0:minLen], evEps[0:minLen], label = elms[ii], lw = 2)
        
        plt.xlabel("TP-AGB time in ky")
        plt.ylabel("Log epsilon")
        plt.ylim([-2, 5])
        
    elif mode == "xfe":
        for ii in range(len(elms)):
            evXFe = evolXFe[ii]
            minLen = min(len(ages), len(evXFe))
            plt.plot(ages[0:minLen], evXFe[0:minLen], label = elms[ii], lw = 2)
        
        plt.xlabel("TP-AGB time in ky")
        plt.ylabel("[X/Fe]")
        plt.ylim([-0.2, 2])
    
    plt.legend(loc = 0)
    plt.show()
    
    return 0

if __name__ == "__main__":
    main()
