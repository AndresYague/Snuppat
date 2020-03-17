import sys
import scipy.stats as sts
import numpy as np
import matplotlib.pyplot as plt

def plotPocketSizes(fname, plotIndep = "lambda"):
    '''Plot pocket sizes against core mass or lambda'''
    # Get omega from title
    omg = float(fname.split("Ov")[1].split(".")[0])*1e-2
    
    # Ignore first x models
    ignoreFirst = 0
    
    # Now do the rest
    pSizes = []; lambs = []; coreMass = []
    intMass = []; integC13 = []; pIntegs = []
    with open(fname, "r") as fread:
        while True:
            line = fread.readline()
            if len(line) == 0:
                break
            
            if "Lambda" in line:
                if ignoreFirst > 0:
                    ignoreFirst -= 1
                    continue
                
                lnlst = line.split(";")
                lamb = float(lnlst[0].split()[-1])
                pockSize = float(lnlst[1].split()[-1])
                cMass = float(lnlst[2].split()[-1])
                iMass = float(lnlst[3].split()[-1])
                
                if pockSize > 0:
                    lambs.append(lamb)
                    coreMass.append(cMass)
                    intMass.append(iMass)
                    pSizes.append(pockSize)
                    
                    intg = 0
                    m1 = None; c131 = None
                    while True:
                        lnlst = map(float, fread.readline().split())
                        if len(lnlst) == 0:
                            break
                        
                        m2 = lnlst[0]
                        c132 = lnlst[1]
                        
                        if m1 is None:
                            m1 = m2; c131 = c132
                            continue
                        else:
                            intg += (m2 - m1)*(c131 + c132)*0.5
                            m1 = m2; c131 = c132
                    
                    pIntegs.append(intg)
        
        # Define dependent
        dep = pSizes
        #dep = pIntegs
        
        # Define independent
        if plotIndep == "lambda":
            indep = map(lambda x: x, lambs)
        elif plotIndep == "core mass":
            indep = map(lambda x: x, coreMass)
        elif plotIndep == "intershell mass":
            indep = map(lambda x: x, intMass)
        
        mm, nn, r_val, p_val, std_err = sts.linregress(indep, dep)
        linPSiz = map(lambda x: mm*x + nn, indep)
        
        lab = "r$^2 = {:.2f}$; ".format(r_val**2)
        lab += "$\omega = {:.2f}$".format(omg)
        lin = plt.plot(indep, dep, "o")
        plt.plot(indep, linPSiz, lin[-1].get_color() + "-",
            label = lab)

def main():
    if len(sys.argv) < 2:
        print("Usage: python {} profile1 <profile2 ...>".format(sys.argv[0]))
        return 1
    
    subplotNum = len(sys.argv[1:])
    fig = plt.figure()
    
    # Choose what to plot against
    plotIndep = "lambda"
    
    # Plot all the profiles
    ii = 0; maxSaveMass = 0
    for arch in sys.argv[1:]:
        plotPocketSizes(arch, plotIndep = plotIndep)
    
    if plotIndep == "lambda":
        plt.xlabel("$\lambda$")
    elif plotIndep == "core mass":
        plt.xlabel("Core mass")
    elif plotIndep == "intershell mass":
        plt.xlabel("Intershell mass")
    
    plt.ylabel("Pocket mass")
    plt.legend(loc = 0)
    plt.show()

if __name__ == "__main__":
    main()
