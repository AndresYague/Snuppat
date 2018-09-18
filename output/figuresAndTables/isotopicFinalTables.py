import sys, math, os
import matplotlib.pyplot as plt

def main():
    # Check that there's at least one argument
    if len(sys.argv) < 2:
        print "Usage python {}".format(sys.argv[0]),
        print "<file1> [<file2> ...]"
        return 1
    
    files = sys.argv[1:]
    
    # Read "species.dat" and store all the values in lists
    species = "../../data/species.dat"
    atomicMass = []; isotNames = []
    with open(species, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            # And names with atomic number
            if lnlst[1] == "d" or lnlst[2] == "0":
                lnlst[1] = "h"
            
            # Add isotopic name
            isotNames.append(lnlst[1] + lnlst[0])
            
            # Now relate positions with atomic mass
            atomicMass.append(int(lnlst[0]))
    
    # Go file by file
    print "Isotopical mass fraction per file"
    for archivo in files:
        
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
        
        # Calculate mass fractions for each element
        massFrac = [(x + y)*0.5*z for (x, y, z) in
                    zip(prevline[4:], newline[4:], atomicMass)]
        
        print "# {}".format(archivo)
        print ""
        
        for ii in range(len(atomicMass)):
            print isotNames[ii], massFrac[ii]

if __name__ == "__main__":
    main()
