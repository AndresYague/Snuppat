import struct, sys
import numpy as np

class readBinaryModels(object):
    '''Class for reading binary models'''
    
    def __init__(self, fil):
        '''Initialize'''
        
        super(readBinaryModels, self).__init__()
        self.fread = open(fil, "rb")
        self.head = None
        self.model = None
    
    def close(self):
        '''Close file'''
        
        self.fread.close()
    
    def __readHeader(self):
        '''Return header'''
        
        head = []
        byte = self.fread.read(4)
        if len(byte) == 0:
            return None
        
        head.append(*struct.unpack('i', byte))
        head.append(*struct.unpack('d', self.fread.read(8)))
        head.append(*struct.unpack('d', self.fread.read(8)))
        head.append(*struct.unpack('i', self.fread.read(4)))
        head.append(*struct.unpack('i', self.fread.read(4)))
        
        return head
    
    def nextModel(self):
        '''Calculate next model, unpacked'''
        
        # Read header
        self.head = self.__readHeader()
        if self.head is None:
            return False
        
        self.model = [self.head]
        for ii in range(self.head[3]):
            s = []
            for jj in range(self.head[4]):
                s.append(*struct.unpack('d', self.fread.read(8)))
            
            self.model.append(s)
        
        return True
    
    def readOnlyHeader(self):
        '''Look only for the header and skip the rest'''
        
        # Read header
        self.head = self.__readHeader()
        if head is None:
            return False
        
        # Skip file
        for ii in range(head[3]):
            for jj in range(head[4]):
                self.fread.read(8)
        
        return True

def getSpeciesDict(species, corrPos = 0):
    '''Return a dictionary of species names and positions'''
    
    ii = 0; speciesDict = {}
    with open(species, "r") as fread:
        for line in fread:
            weight, name = line.split()[0:2]
            name = "{}{}".format(name, weight)
            
            speciesDict[name] = (ii + corrPos, int(weight))
            ii += 1
    
    return speciesDict

def getCoreMass(model, speciesDict):
    '''Return the core mass of this model, defined
    as the point where H < 0.1*maxH'''
    
    hindx = speciesDict["p1"][0]
    mass = []; hydro = []; temp = []
    h2 = None; m2 = None; t2 = None
    for ii in range(model[0][3]):
        h1 = model[ii + 1][hindx]
        m1 = model[ii + 1][0]
        t1 = model[ii + 1][1]*1e3
        if h2 is None:
            h2 = h1; m2 = m1; t2 = t1
            continue
        
        hydro.append((h1 + h2)*0.5)
        mass.append((m1 + m2)*0.5)
        temp.append((t1 + t2)*0.5)
        h2 = h1; m2 = m1; t2 = t1
    
    maxHydr = max(hydro)
    ii = len(hydro) - 1
    while ii > -1:
        if hydro[ii] <= 1e-1*maxHydr:
            return mass[ii]*model[0][1], temp[ii]
        
        ii -= 1
    
    return None, None

def getBCEMass(model, speciesDict):
    '''Return the BCE mass of this model, defined
    as the last convective poin where H ~ maxH'''
    
    hindx = speciesDict["p1"][0]
    mass = []; hydro = []; temp = []; radiat = []
    h2 = None; m2 = None; t2 = None; rad2 = None
    for ii in range(model[0][3]):
        h1 = model[ii + 1][hindx]
        m1 = model[ii + 1][0]
        t1 = model[ii + 1][1]*1e3
        rad1 = model[ii + 1][3]
        if h2 is None:
            h2 = h1; m2 = m1; t2 = t1; rad2 = rad1
            continue
        
        hydro.append((h1 + h2)*0.5)
        mass.append((m1 + m2)*0.5)
        temp.append((t1 + t2)*0.5)
        radiat.append((rad1 + rad2)*0.5)
        h2 = h1; m2 = m1; t2 = t1; rad2 = rad1
    
    maxHydr = max(hydro)
    ii = len(hydro) - 1
    bceMass, bceTemp = None, None
    while ii > -1:
        if hydro[ii] >= (maxHydr*.9) and radiat[ii] > 0:
            bceMass = mass[ii]*model[0][1]
            bceTemp = temp[ii]
        
        ii -= 1
    
    return bceMass, bceTemp

def getMaximumTemperature(model):
    '''Return the maximum temperature'''
    
    temp = []
    t2 = None
    for ii in range(model[0][3]):
        t1 = model[ii + 1][1]*1e3
        if t2 is None:
            t2 = t1
            continue
        
        temp.append((t1 + t2)*0.5)
        t2 = t1
    
    return max(temp)

def getIntershellComposition(model, coreMass, wantedElements, speciesDict):
    '''Return the intershell wantedElements composition'''
    
    # Make sure he4 and c12 are in the list, store their index
    cpyElements = wantedElements[:]
    if "he4" not in cpyElements:
        cpyElements.append("he4")
    if "c12" not in cpyElements:
        cpyElements.append("c12")
    indxHe4 = cpyElements.index("he4")
    indxC12 = cpyElements.index("c12")
    
    # Store indices and weights
    elemIndxWeight = []
    for element in cpyElements:
        elemIndxWeight.append(speciesDict[element])
    
    # Prepare for reading
    mass = []; m2 = None
    composition = [[] for x in cpyElements]
    for ii in range(model[0][3]):
        # Store mass
        m1 = model[ii + 1][0]
        
        # Store elements
        for jj in range(len(elemIndxWeight)):
            val = model[ii + 1][elemIndxWeight[jj][0]]
            val *= elemIndxWeight[jj][1]
            composition[jj].append(val)
            
            # If not in first case, go to previous values and correct
            if ii > 0:
                composition[jj][ii - 1] = (composition[jj][ii - 1] + val)*0.5
        
        # If in first case, skip
        if ii == 0:
            m2 = m1
            continue 
       
        mass.append((m1 + m2)*0.5)
        m2 = m1
    
    # After the last step, we have to drop the last value in the composition list
    for elem in composition:
        elem.pop()
    
    # Now look for intershell value
    totMass = 0
    maxHe = max(composition[indxHe4])
    averageComposition = [0 for x in composition]
    for ii in range(len(mass)):
        he4 = composition[indxHe4][ii]
        c12 = composition[indxC12][ii]
        if abs(he4 - maxHe) < 0.4*maxHe:
            if he4 > c12 and c12 > 0.1:
                for jj in range(len(composition)):
                    averageComposition[jj] += composition[jj][ii]*mass[ii]
                totMass += mass[ii]
    
    try:
        for ii in range(len(averageComposition)):
            averageComposition[ii] /= totMass
    except ZeroDivisionError:
        pass
    except:
        raise
    
    return averageComposition

def getTDUMass(binObj, speciesDict, wantedElements, massTh, ageThreshold):
    '''Finds next core masses'''
    
    # For storage
    foundTDU = False
    inTDU = False
    totMass = None
    maxCoreMass = 0
    minCoreMass = 0
    maxTDUTemp = 0
    maxTempBCE = 0
    maxTempHyd = 0
    maxComposition = [0 for elem in wantedElements]
    prevModNum = 0
    
    # Read model
    tduAge = None
    while True:
        newModel = binObj.nextModel()
        if not newModel:
            break
        
        model = binObj.model
        
        # Age
        currAge = 10**binObj.head[2]
        
        # Get core and BCE mass
        coreMass, tempHyd = getCoreMass(model, speciesDict)
        BCEMass, tempBCE = getBCEMass(model, speciesDict)
        if coreMass is None:
            raise Exception("Core mass not found!")
        if BCEMass is None:
            raise Exception("BCE mass not found!")
        
        if coreMass > maxCoreMass and not inTDU:
            maxCoreMass = coreMass
            minCoreMass = coreMass
        elif (maxCoreMass - coreMass) > massTh:
            inTDU = True
            foundTDU = True
            tduAge = currAge
            if coreMass < minCoreMass:
                minCoreMass = coreMass
                totMass = model[0][1]
            
        elif foundTDU and (currAge - tduAge) > ageThreshold:
            inTDU = False
        
        # Get intershell composition
        intComposition = getIntershellComposition(model, coreMass,
                                                  wantedElements, speciesDict)
        if None in intComposition:
            for ii in range(len(wantedElements)):
                if intComposition[ii] is None:
                    print("{} value not found!".format(wantedElements[ii]))
            
            raise Exception()
        
        # Get maximum intershell cmposition
        maxComposition = np.maximum(maxComposition, intComposition)
        
        # Get maximum temperatures
        if tempHyd > maxTempHyd:
            maxTempHyd = tempHyd
        
        if tempBCE > maxTempBCE:
            maxTempBCE = tempBCE
        
        maxTemp = getMaximumTemperature(model)
        if maxTemp > maxTDUTemp:
            maxTDUTemp = maxTemp
        
        if not inTDU and foundTDU:
            break
        
        # Exit if repeating models
        if prevModNum == model[0][0]:
            break
        
        prevModNum = model[0][0]
    
    val = (maxCoreMass, minCoreMass, totMass, maxTDUTemp, maxTempBCE,
            maxTempHyd, maxComposition)
    return val

def main():
    '''Main program'''
    
    if len(sys.argv) < 2:
        print("Usage: python {} <SnuppatOutputBIN>".format(sys.argv[0]), end = " ")
        print("[mass threshold [age threshold]]")
        return 1
    
    # Get mass threshold
    if len(sys.argv) > 2:
        massTh = float(sys.argv[2])
    else:
        massTh = 4e-5
    
    # Get age threshold
    if len(sys.argv) > 3:
        ageThreshold = float(sys.argv[3])
    else:
        ageThreshold = 1e3
    
    # Wanted elements
    wantedElements = ["he4", "c12", "o16", "ne22"]
    
    species = "../../data/species.dat"
    speciesDict = getSpeciesDict(species, 4)
    
    # Create object
    binObj = readBinaryModels(sys.argv[1])
    
    # Print header
    s = "Lambda, Temperature, BCE Temperature, Hyd Temperature, Core Mass, "
    s += "Envelope mass, Dredged Mass"
    for elem in wantedElements:
        s += ", " + elem.upper()
    s += "; Mass Threshold (Msun) = {:.2E}".format(massTh)
    s += " Age Threshold (years) = {:.2E}".format(ageThreshold)
    print(s)
    
    # Look for max and min core mass
    prevMinCore = None
    while True:
        val = getTDUMass(binObj, speciesDict, wantedElements, massTh, ageThreshold)
        maxCoreM, minCoreM, totMass, temp, tempBCE, tempHyd, composition = val
        if maxCoreM == 0 or totMass is None:
            break
        
        if prevMinCore is None:
            prevMinCore = minCoreM
            continue
        
        coreGrowth = maxCoreM - prevMinCore
        dredg = maxCoreM - minCoreM
        lambd = dredg/coreGrowth
        
        prevMinCore = minCoreM
        
        values = [lambd, temp, tempBCE, tempHyd, minCoreM, totMass - minCoreM,
                  dredg]
        for comp in composition:
            values.append(comp)
        
        for ii in range(len(values)):
            if ii == 0:
                s = "{:.2E}".format(values[ii])
            else:
                s += " {:.2E}".format(values[ii])
        print(s)

if __name__ == "__main__":
    main()
