'''Program to create an evolution graphic for different magnitudes'''

import sys, os, math
import scipy.stats as stats
import matplotlib.pyplot as plt

def element_position(elem, correct):
    '''Function that searchs for a given element and correction offset,
    returning a tuple of lists with the name of the element, its atomic
    mass and the position. This format is used because is possible to ask
    for all the isotopes of an element.
    
    Input can be an element name or the name and atomic mass of the
    specific isotope. Examples would be "ne22" or "fe"'''
    
    # See if there's a number in the name
    theres_num = False
    if (elem[-1].isdigit()):
        theres_num = True
    
    # Species file
    spec = os.path.join("..", "..", "data", "species.dat")
    
    # Get position
    posit = list()
    num = list()
    with open(spec, "r") as fread:
        ii = 1
        for line in fread:
            lnlst = line.split()
            
            if theres_num:
                stri = "{}{}".format(lnlst[1], lnlst[0])
            else:
                stri = "{}".format(lnlst[1])
            
            # Compare
            if theres_num:
                if elem == stri:
                    posit.append(ii + correct - 1)
                    num.append(int(lnlst[0]))
                    break
            else:
                if elem == stri:
                    posit.append(ii + correct - 1)
                    num.append(int(lnlst[0]))
            
            ii += 1
    
    # Check if element found
    if len(posit) == 0:
        raise Exception("{} not in the list".format(elem))
    
    # Capitalize name
    if theres_num:
        name = "{}".format(lnlst[1])
    else:
        name = "{}".format(elem)
    
    name = name[0].upper() + name[1:]
    
    return (name, num, posit)

def main():
    '''Create evolution graphic after selecting input'''
    
    # Manage input
    if len(sys.argv) < 3:
        print "Use: python {} <simulation file> <element> ".format(sys.argv[0])
        return 1
    
    arch = sys.argv[1]
    elem = sys.argv[2].lower()
    
    # Read input file
    inpt = []; inptErr = False
    with open("input.in") as fIn:
        for line in fIn:
            lnlst = line.split()
            inpt.append(lnlst[0])
    
    # Check correctness
    if len(inpt) < 4:
        inptErr = True
    else:
        mode = inpt[0].lower()
        compare = inpt[1].lower()
        outpt = inpt[2]
        smooth = inpt[3].lower()
        
        if mode not in ["surface", "intershell", "maximum"]:
            inptErr = True
        
        if smooth not in ["y", "n"]:
            inptErr = True
    
    if inptErr:
        print "Wrongly formatted input file, please check Readme.txt"
        return 1
    
    smooth = True if smooth == "y" else False
    
    # Define surface mass
    surf_mass = 0.95
    
    # Position correction
    correct = 4
    
    # Temperature position
    tempos = 1
    
    # Hydrogen, ne22, he4 and iron positions and mass
    name, hMass, posHid = element_position("p1", correct)
    name, neMass, posNe = element_position("ne22", correct)
    name, heMass, posHe = element_position("he4", correct)
    name, cMass, posC = element_position("c12", correct)
    name, cmpMass, posCmp = element_position(compare, correct)
    posHid = posHid[0]; posNe = posNe[0]; posHe = posHe[0]; posC = posC[0]
    hMass = hMass[0]; neMass = neMass[0]; heMass = heMass[0]; cMass = cMass[0]
    
    # Get element name and position
    name, num, posit = element_position(elem, correct)
    
    # Read data
    with open(arch, "r") as fread:
        time = list(); data = list()
        firstage = None; lnlst = None
        while True:
            # Look for the next model
            while True:
                # Read next line
                if lnlst is None or "#" not in lnlst:
                    line = fread.readline()
                    lnlst = line.split()
                
                if len(line) == 0:
                    break
                
                # Check if new model
                if len(lnlst) > 0 and lnlst[0] == "#" and lnlst[1] == "Model":
                    # Read next two lines and exit:
                    mass = float(fread.readline().split()[-1])
                    age = float(fread.readline().split()[-1])
                    
                    # Store AGB starting age
                    if firstage is None:
                        firstage = 10**(age - 6)
                    
                    break
            
            if len(line) == 0:
                break
            
            # Look for condition and store values
            lnlst = fread.readline().split()
            maxVal = None; prevHe4 = None
            while True:
                prevlst = lnlst
                lnlst = fread.readline().split()
                
                if mode == "maximum":
                    if "Model" in lnlst or len(lnlst) == 0:
                        prevlst = maxPrevlst
                        lnlst = maxlst
                        break
                    
                    valelm = 0.
                    for ii in range(len(posit)):
                        newval = float(prevlst[posit[ii]])
                        newval += float(lnlst[posit[ii]])
                        valelm += newval*0.5
                    
                    # Store the lines of interest
                    if maxVal is None or valelm > maxVal:
                        maxVal = valelm
                        maxPrevlst = prevlst
                        maxlst = lnlst
                    
                elif mode == "surface":
                    # Check if bigger than wanted mass
                    if float(lnlst[0]) > surf_mass:
                        break
                    
                elif mode == "intershell":
                    # Check where does the ne22 falls below 1e-3 (both are
                    # multiplied by 0.5 and their mass for mass fraction).
                    valc12 = (float(prevlst[posC]) + float(lnlst[posC])) \
                              *0.5*cMass
                    valhe4 = (float(prevlst[posHe]) + float(lnlst[posHe])) \
                             *0.5*heMass
                    
                    if valhe4 > valc12:
                        break
            
            # Get the value for the element:
            valelm = 0.
            for ii in range(len(posit)):
                newval = float(prevlst[posit[ii]]) + float(lnlst[posit[ii]])
                valelm += newval*0.5
            
            if valelm > 0.:
                time.append(10**(age - 6) - firstage)
                print time[-1], valelm
                
                # Choose representation
                if mode == "surface":
                    # Hydrogen value
                    valhid = (float(prevlst[posHid]) + float(lnlst[posHid]))*0.5
                    
                    # Comparison value
                    if compare != "p1":
                        valCmp = 0
                        for indCmp in posCmp:
                            newval = float(prevlst[indCmp])
                            newval += float(lnlst[indCmp])
                            valCmp += newval*0.5
                    
                    # Add the one relevant
                    if compare == "p1":
                        data.append(math.log10(valelm/valhid) + 12)
                    else:
                        if len(data) == 0:
                            firstCmpVal = valelm/valCmp
                        
                        data.append(math.log10((valelm/valCmp)/firstCmpVal))
                    
                elif mode == "intershell" or mode == "maximum":
                    data.append(valelm)
    
    # Clean input from numerical artifacts
    
    # To do that, calculate the local mean and the standard deviation
    # within 2*step points and eliminate every point above sigRej sigma.
    # Then repeat until no longer happens.
    sigRej = 2; step = len(data)/20
    while smooth:
        if step == 0:
            break
        
        smooth = False
        
        newData = []; newTime = []
        for ii in range(len(data)):
            # Define bin limits
            before = ii - step
            after = ii + step
            
            # Fringe cases
            if before < 0:
                after -= before
                before = 0
            elif after > len(data):
                before -= len(data) - after
                after = len(data)
            
            dataMean = stats.tmean(data[before:after])
            dataSigma = math.sqrt(stats.moment(data[before:after], moment = 2))
            
            if data[ii] < (dataMean + sigRej*dataSigma):
                newData.append(data[ii])
                newTime.append(time[ii])
            else:
                smooth = True
        
        data = newData
        time = newTime
    
    # Plot
    plt.rcParams.update({"font.size": 14})
    plt.plot(time, data, linewidth = 2)
    plt.tight_layout(pad=2.0, w_pad=4.0, h_pad=4.0)
    plt.xlabel("Time in AGB (Myrs)")
    
    # Get the comparison element in upper case
    if compare == "p1":
        upComp = compare
    else:
        upComp = compare[0].upper() + compare[1:]
    
    # Select scales and labels
    if mode != "surface":
        plt.yscale("log")
    
    if len(num) == 1:
        if mode == "surface":
            # Compare with hydrogen or iron
            if compare == "p1":
                plt.ylabel("Log $\epsilon(^{{{}}}{})$".format(num[0], name))
            else:
                plt.ylabel("$[^{{{}}}{}/{}]$".format(num[0], name, upComp))
            
            plt.title("{}{} over {}".format(name,num[0], upComp))
            plt.savefig("{}{}_over_{}{}.png".format(name, num[0], upComp, outpt))
            
        elif mode == "intershell" or mode == "maximum":
            plt.ylabel("$^{{{}}}${} Number fraction".format(num[0], name))
            
            plt.title("{}{}".format(name, num[0]))
            plt.savefig("{}{}{}.png".format(name, num[0], outpt))
        
        
    else:
        if mode == "surface":
            # Compare with hydrogen or iron
            if compare == "p1":
                plt.ylabel("Log $\epsilon({})$".format(name))
            else:
                plt.ylabel("[{}/{}]".format(name, upComp))
            
            plt.title("{} over {}".format(name, upComp))
            plt.savefig("{}_over_{}{}.png".format(name, upComp, outpt))
            
        elif mode == "intershell" or mode == "maximum":
            plt.ylabel("{} Number fraction".format(name))
            
            plt.title("{}".format(name))
            plt.savefig("{}{}.png".format(name, outpt))
        
    
    plt.show()

if __name__ == "__main__":
    main()
