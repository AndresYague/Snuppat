import sys, os, math, numpy

class SpeciesDict(dict):
    '''Redefine __getitem__ from dict so we can get whole species. This
    dictionary always returns a list'''
    
    def __init__(self, *args, **kw):
        '''Get all the usual proprieties from a dictionary'''
        
        super(SpeciesDict, self).__init__(*args, **kw)
    
    def __setitem__(self, key, val):
        '''Extend __setitem__ creating one-element lists for isotopes and
        multi-element lists for species'''
        
        # Add exception for neutrons
        if key == "n1":
            super(SpeciesDict, self).__setitem__(key, [val])
            return
        
        # First, get species name
        # Suppose format a[bc]<numbers>
        name = ""
        for lett in key:
            if lett in "0123456789":
                break
            else:
                name += lett
        
        # Now add value to this isotope
        super(SpeciesDict, self).__setitem__(key, [val])
        
        # Now append this value to list
        if name in self:
            lst = self[name]
            lst.append(val)
            super(SpeciesDict, self).__setitem__(name, lst)
        else:
            super(SpeciesDict, self).__setitem__(name, [val])

class ReadObject(object):
    '''Object to read TDU data from simulations'''
    
    def __init__(self, inptFile, speciesDict, showSpecies):
        
        self.inptFile = open(inptFile, "r")
        self.speciesDict = speciesDict
        self.showSpecies = showSpecies
        self.nTDU = 0
        self.lastModels = None
        self.firstAge = None
    
    def getNextTDU(self):
        '''Get next TDU from the models'''
        
        # Initialize
        if self.lastModels is not None:
            TDUModels = self.lastModels[:]
            firstCore = TDUModels[0][4]
            lastCore = TDUModels[1][4]
            maxCore = TDUModels[1][4]
        else:
            TDUModels = []
            firstCore = None
            lastCore = None
            maxCore = None
        
        inTDU = False
        minLastCore = None
        maxC13DeltMass = 0; maxC13Height = 0; maxC13Area = 0
        while not inTDU:
            '''Store each model of this TDU'''
            
            nextModel = self.__getModel()
            
            # Exit if last model
            if nextModel is None:
                self.lastModels = None
                break
            
            # Get age
            if self.firstAge is None:
                self.firstAge = 10**(nextModel[2] - 6)
            
            line2 = None
            age = 10**(nextModel[2] - 6) - self.firstAge
            mass = []; temp = []; rad = []; p1 = []; he4 = []
            c13 = []; n14 = []; showSpecLst = []
            for line in nextModel[3]:
                
                # If we don't have a second line, continue
                if line2 is None:
                    line2 = line
                    continue
                
                # Store the relevant quantities
                mass.append((line[0] + line2[0])*0.5*nextModel[1])
                temp.append((line[1] + line2[1])*5e8)
                rad.append((line[3] + line2[3])*0.5)
                p1.append((line[self.speciesDict["p1"][0]] +
                    line2[self.speciesDict["p1"][0]])*0.5)
                he4.append((line[self.speciesDict["he4"][0]] +
                    line2[self.speciesDict["he4"][0]])*0.5*4)
                c13.append((line[self.speciesDict["c13"][0]] +
                    line2[self.speciesDict["c13"][0]])*0.5*13)
                n14.append((line[self.speciesDict["n14"][0]] +
                    line2[self.speciesDict["n14"][0]])*0.5*14)
                
                ii = 0
                for species in self.showSpecies:
                    val = 0
                    for indx in self.speciesDict[species]:
                        val += (line[indx] + line2[indx])*0.5
                    
                    if len(showSpecLst) < len(self.showSpecies):
                        showSpecLst.append([val])
                    else:
                        showSpecLst[ii].append(val)
                    
                    ii += 1
                
                line2 = line
            
            # Get proton envelope value
            p1Env = self.__getEnvP1(p1, rad, mass)
            
            # Get c13 pocket width and height
            c13Parameters = self.__getEffectiveC13(p1, c13, n14, mass)
            (c13Area, c13Mass, c13Profile) = c13Parameters
            if c13Mass is not None:
                c13DeltMass = c13Mass[-1] - c13Mass[0]
                c13Height = max(c13Profile)
            else:
                c13DeltMass = 0
                c13Height = 0
            
            # Get the biggest one
            if c13DeltMass > maxC13DeltMass:
                maxC13Area = c13Area
                maxC13DeltMass = c13DeltMass
                maxC13Height = c13Height
            
            # Store envelope borders
            for ii in range(len(mass)):
                if abs(p1[ii] - p1Env)/p1Env < 1e-2:
                    envBord = mass[ii]
                    break
            
            nMod = nextModel[0]
            TDUModels.append((nMod, age, mass, temp, envBord, p1, he4, c13Mass,
                c13Profile, showSpecLst))
            
            # Store firstCore and cycle
            if firstCore is None:
                firstCore = envBord
                lastCore = envBord
                maxCore = envBord
                continue
            
            # Store lastCore and maxCore
            lastCore = envBord
            if envBord > maxCore:
                maxCore = envBord
            
            # If the core has not grown enough, cycle
            if (maxCore - firstCore) < 1e-4:
                continue
            
            # Check if lastCore is shrinking
            if minLastCore is None:
                if lastCore < maxCore:
                    minLastCore = lastCore
                
            elif minLastCore >= lastCore:
                minLastCore = lastCore
            else:
                # If it stopped shrinking, exit and get lambda
                inTDU = True
                self.nTDU += 1
                self.lastModels = TDUModels[-2:]
                TDUModels.pop()
        
        # If at the end, close file
        if len(TDUModels) == 0:
            self.inptFile.close()
            return None
        
        if minLastCore is None:
            minLastCore = lastCore
        
        # Calculate lambda, dredge mass and ratio of dredged mass and
        # envelope mass
        try:
            lambd = (maxCore - minLastCore)/(maxCore - firstCore)
            dredMass = maxCore - minLastCore
            envMass = mass[-1] - minLastCore
        except ZeroDivisionError:
            lambd = 0
            dredMass = 0
            envMass = 0
        except:
            raise
        
        return (TDUModels, lambd, dredMass, minLastCore, envMass,
                maxC13DeltMass, maxC13Height)
    
    def __getModel(self):
        '''Read model'''
        
        # Read model header
        line = self.inptFile.readline()
        
        # Exit if last model
        if len(line) == 0:
            return None
        
        mass = self.inptFile.readline().split()
        age = self.inptFile.readline().split()
        
        modNum = int(line.split()[-1])
        mass = float(mass[-1])
        age = float(age[-1])
        
        # Read chemistry values
        # mass, temp, rho, radiat, abundances
        modelStore = []
        while True:
            prevPosition = self.inptFile.tell()
            line = self.inptFile.readline().split()
            
            if len(line) == 0:
                break
            elif line[0] == "#":
                self.inptFile.seek(prevPosition)
                break
            
            lnlst = map(lambda x: float(x), line)
            modelStore.append(lnlst)
        
        # Check if pulse
        return (modNum, mass, age, modelStore)
    
    def __getEnvP1(self, p1, rad, mass):
        '''With the p1 list, identify the convective envelope value'''
        
        # Store total mass
        totMass = mass[-1]
        
        # Return if inside of a convective zone and below 95% of total mass
        ii = len(p1) - 2
        while ii > 0:
            if rad[ii] > 0:
                if p1[ii] == p1[ii - 1] and p1[ii] == p1[ii + 1]:
                    if mass[ii] < 0.95*totMass:
                        return p1[ii]
            
            ii -= 1
        
        return p1[-1]
    
    def __getEffectiveC13(self, p1, c13, n14, mass):
        '''Identify the effective c13 area'''
        
        effC13Area = 0
        massRange = None
        effC13Profile = None
        
        # Look for c13 pocket
        seen = False
        for ii in range(len(mass)):
            if p1[ii] < 1e-10 and c13[ii] > 1e-4 and c13[ii]/13 > n14[ii]/14:
                if not seen:
                    low = ii - 1
                
                high = ii + 1
                seen = True
        
        # Ensure there is at least 1 point
        if seen and low == high:
            seen = False
        
        # Get Area
        if seen:
            effC13Area = quadr(mass, c13, n14, low, high)
            massRange = mass[low:high + 1]
            effC13Profile = [x - y*13./14 for (x, y) in
                    zip(c13[low:high + 1], n14[low:high + 1])]
            
            # TODO
            print max(effC13Profile)
        
        return (effC13Area, massRange, effC13Profile)

def createSpeciesRelation(species):
    '''Return a dictionary with indexes for species'''
    
    # Open and read
    correction = 4
    speciesDict = SpeciesDict()
    with open(species, "r") as fread:
        i = 0
        for line in fread:
            lnlst = line.split()
            
            # Add specie
            name = lnlst[1] + lnlst[0]
            speciesDict[name] = i + correction
            
            # Advance index
            i += 1
    
    return speciesDict

def quadr(xx, yy, yy2 = None, low = None, high = None):
    '''Quadrature method'''
    
    # Manage input
    if low is None:
        low = 0
    if high is None:
        high = len(xx) - 1
    if yy2 is None:
        yy2 = [0 for x in xx]
    
    # Now apply trapezoidal rule
    ii = low; sol = 0
    while ii < high - 1:
        dx = xx[ii + 1] - xx[ii]
        sol += (yy[ii + 1] + yy[ii] - yy2[ii + 1] - yy2[ii])*dx*0.5
        
        # Advance ii
        ii += 1
    
    return sol

def getReferenceValues(showSpecies, solarVals):
    '''Calculate reference values'''
    
    # Open solarVals and get dictionary
    di = {}
    with open(solarVals, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            name = lnlst[0]
            zz = int(lnlst[2])
            value = float(lnlst[1])
            
            di[name] = value
    
    refs = []
    for name in showSpecies:
        if name[-1] in "0123456789":
            raise Exception("{} is an isotopic value, stopping".format(name))
        elif name in di:
            val = di[name]
            refs.append([(val, val)])
        else:
            raise Exception("{} not in solar values".format(name))
    
    return refs

def extractAbund(abundances, TDUModel):
    mass = TDUModel[2]
    envBord = TDUModel[4]
    values = TDUModel[-1]
    
    # Get values within a shell of 2e-2 solar masses from the envelope border
    bef = envBord - 1e-2
    aft = envBord + 1e-2
    
    # Go element by element, this can be optimized
    for jj in range(len(values)):
        befVal = 0
        aftVal = None
        
        nVals = 0
        for ii in range(len(mass)):
            if mass[ii] < envBord and mass[ii] > bef:
                befVal += values[jj][ii]
                nVals += 1
            
            if mass[ii] > envBord and mass[ii] < aft:
                if aftVal is None:
                    aftVal = values[jj][ii]
                else:
                    aftVal = min(aftVal, values[jj][ii])
        
        # Fix befVal and aftVal
        befVal /= nVals
        if aftVal is None:
            aftVal = values[jj][ii]
        
        abundances[jj].append((befVal, aftVal))
    
    # Now calculate the last abundances knowing that we have the references in
    # the first values
    lastVals = []
    refs = [x[0][0] for x in abundances]
    for ii in range(len(refs)):
        lastVals.append((math.log10(abundances[ii][-1][0]/refs[ii]),
            math.log10(abundances[ii][-1][1]/refs[ii])))
    
    return lastVals

def getLsHs(abundances):
    '''Obtain ls and hs'''
    
    # Calculate hs/ls knowing that the 6 last values are what interest us
    cpyAbnd = abundances[-6:]
    
    lastVals = []
    refs = [x[0][0] for x in cpyAbnd]
    for ii in range(len(refs)):
        lastVals.append((math.log10(cpyAbnd[ii][-1][0]/refs[ii]),
            math.log10(cpyAbnd[ii][-1][1]/refs[ii])))
    
    ls = [reduce(lambda x, y: x + y, lst) for lst in zip(*lastVals[:3])]
    hs = [reduce(lambda x, y: x + y, lst) for lst in zip(*lastVals[3:])]
    ls = map(lambda x: x/3., ls)
    hs = map(lambda x: x/3., hs)
    
    return (ls, hs)

def main():
    '''Make the TDU tables'''
    
    # Read input
    if len(sys.argv) < 2:
        print "Usage: python {} <input file> [output file]".format(sys.argv[0])
        return 1
    
    inptFile = sys.argv[1]
    
    if len(sys.argv) > 2:
        outpFile = sys.argv[2]
        c13ProfFile = "C13Profile" + outpFile
    else:
        outpFile = None
        c13ProfFile = None
    
    # Species info
    species = os.path.join("..", "..", "data", "species.dat")
    solarVals = os.path.join("..", "..", "data", "solarVals.dat")
    speciesDict = createSpeciesRelation(species)
    
    showSpecies = ["rb", "sr", "ba", "pb"]
    
    # Add the 6 species for calculating hs/ls
    showSpecies += ["sr", "y", "zr", "ba", "la", "ce"]
    
    # Create data object
    data = ReadObject(inptFile, speciesDict, showSpecies)
    s = "TDU # | Age range | Max T9 | Lambda | Dredged mass | Core mass | "
    s += "Envelope mass | Dredged ratio | Total dredged | "
    s += "Effective c13 (width, height)"
    for sp in showSpecies[:-6]:
        s += " | {}".format(sp)
    
    # Add hs/ls
    s += " | ls | hs"
    
    # Write and print
    if outpFile is not None:
        with open(outpFile, "w") as fwrite:
            fwrite.write(s + '\n')
    
    print s
    
    # Get next TDU
    eventNum = 0; totDredge = 0
    while True:
        eventNum += 1
        
        # Get next TDU
        dat = data.getNextTDU()
        if dat is not None:
            TDUModels, lambd, dredMass, coreMass  = dat[:4]
            envMass, c13Mass, c13Height = dat[4:]
        else:
            break
        
        # Get maximum temperature
        maxTemp = 0
        for model in TDUModels:
            for ii in range(len(model[2])):
                temp = model[3][ii]
                he4Val = model[6][ii]
                
                if he4Val < 1e-4:
                    continue
                
                if temp > maxTemp:
                    maxTemp = temp
        
        # Print c13 profile
        if c13ProfFile is not None:
            flag = "w" if eventNum == 1 else "a"
            with open(c13ProfFile, flag) as fwrite:
                # Write first TDU event number
                fwrite.write("# TDU event: {}\n".format(eventNum))
                for model in TDUModels:
                    c13MassProf = model[-3]
                    c13Profile = model[-2]
                    if c13MassProf is not None:
                        fwrite.write("# Profile\n")
                        for ii in range(len(c13MassProf)):
                            fwrite.write("{} {}\n".format(c13MassProf[ii],
                                                        c13Profile[ii]))
        
        # Get reference for the abundances
        if eventNum == 1:
            abundances = getReferenceValues(showSpecies, solarVals)
        
        # Get abundances
        lastAbund = extractAbund(abundances, TDUModels[-1])
        
        # Get hs/ls
        ls, hs = getLsHs(abundances)
        
        try:
            ratioDredge = dredMass/envMass
        except ZeroDivisionError:
            ratioDredge = 0
        except:
            raise
        
        totDredge += dredMass
        
        s = "{:2d}  {:>6f} - ".format(eventNum, TDUModels[0][1])
        s += "{:<6f}  {:10.3e}  ".format(TDUModels[-1][1], maxTemp*1e-9)
        s += "{:10.3e}  {:10.3e}  ".format(lambd, dredMass)
        s += "{:10.3e}  {:10.3e}  ".format(coreMass, envMass)
        s += "{:10.3e}  {:10.3e}  ".format(ratioDredge, totDredge)
        s += "({:10.3e}, {:10.3e})".format(c13Mass, c13Height)
        for tup in lastAbund[:-6]:
            s += "  ({:10.3e}, {:10.3e})".format(*tup)
        s += "  ({:10.3e}, {:10.3e})".format(*ls)
        s += "  ({:10.3e}, {:10.3e})".format(*hs)
        
        # Write and print
        if outpFile is not None:
            with open(outpFile, "a") as fwrite:
                fwrite.write(s + '\n')
        
        print s

if __name__ == "__main__":
    main()
