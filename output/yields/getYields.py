import sys, os

def getMassFracs(abunds, specs):
    '''Multiply number fraction abundances by mass for mass fraction'''
    
    massFracs = {}
    for ii in range(len(specs)):
        nam, zz, mass = specs[ii]
        massFracs[zz] = massFracs.get(zz, 0) + abunds[ii]*mass
    
    return massFracs

def main():
    '''Program to extract yields from the SNUPPAT models'''
    
    # First check number of arguments
    if len(sys.argv) < 2:
        print "Usage: python {} <input file>".format(sys.argv[0])
        return 1
        
    else:
        outFil = sys.argv[1]
    
    # Species dictionary
    specFile = os.path.join("..", "..", "data", "species.dat")
    specs = {}; zzNam = {}; ii = 0
    with open(specFile) as fread:
        for line in fread:
            lnlst = line.split()
            
            zz = int(lnlst[0]) - int(lnlst[2])
            nam = lnlst[1]; mass = int(lnlst[0])
            
            if zz == 1:
                nam = "h"
            
            specs[ii] = (nam, zz, mass)
            zzNam[zz] = nam
            ii += 1
    
    firstAge = None; prevMass = None
    yields = None; prevMassFrac = None
    with open(outFil) as fread:
        for line in fread:
            lnlst = line.split()
            
            # Check first mass and age
            if "#" in line:
                if "Mass" in line:
                    totMass = float(lnlst[-1])
                    if prevMass is None:
                        prevMass = totMass
                    
                    dMass = prevMass - totMass
                    
                elif "Age" in line:
                    age = 10**(float(lnlst[-1]) - 3)
                    if firstAge is None:
                        firstAge = age
                
                prevLine = None
                foundLine = False
                continue
            
            if foundLine:
                continue
            
            # Now check each line
            if prevLine is None:
                prevLine = lnlst
                continue
            
            # Calculate mass
            mass = (float(lnlst[0]) + float(prevLine[0]))*0.5
            
            if mass > 0.85:
                # Get chemistry
                lst1 = [float(x) for x in lnlst]
                lst2 = [float(x) for x in prevLine]
                abunds = [(x + y)*0.5 for x, y in zip(lst1[4:], lst2[4:])]
                
                # Calculate mass fractions from abundances
                massFracDic = getMassFracs(abunds, specs)
                massFracs = massFracDic.values()
                
                # Calculate yields
                if yields is None:
                    # Initialize yields array and store initial
                    # mass fraction
                    yields = [0 for x in massFracs]
                    initMassFrac = massFracs
                    prevMassFrac = [0 for x in massFracs]
                else:
                    # Do the yields integral
                    massFracs = [x - y for x, y in zip(massFracs, initMassFrac)]
                    
                    # Apply trapezoidal rule
                    if prevMassFrac is not None:
                        tempFrac = [x + y for x, y in zip(massFracs, prevMassFrac)]
                        yields = [x + y*dMass*0.5 for x, y in zip(yields, tempFrac)]
                    
                    prevMassFrac = massFracs
                
                foundLine = True
                prevMass = totMass
            
            # Print results
            if foundLine:
                dt = age - firstAge
                if dt > 0:
                    print "# Mass: {} MSun. Time: {} ky".format(totMass, dt)
                    for ii in range(len(zzNam)):
                        print "{:2} {:2} {:11.4E}".format(zzNam[ii], ii, yields[ii])
                    
                    print
            
            prevLine = lnlst

if __name__ == "__main__":
    main()
