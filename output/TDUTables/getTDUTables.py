import struct, sys

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
    
    ii = 0; spD = {}
    with open(species, "r") as fread:
        for line in fread:
            wgt, nam = line.split()[0:2]
            nam = "{}{}".format(nam, wgt)
            
            spD[nam] = ii + corrPos
            ii += 1
    
    return spD

def getCoreMass(model, speciesDict):
    '''Return the core mass of this model, defined
    as the point where H < 0.1*maxH'''
    
    hindx = speciesDict["p1"]
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
    
    hindx = speciesDict["p1"]
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

def getIntershellC12O16(model, coreMass, speciesDict):
    '''Return the intershell c12 and o16'''
    
    # Store values
    c12Indx = speciesDict["c12"]
    o16Indx = speciesDict["o16"]
    he4Indx = speciesDict["he4"]
    mass = []; c12 = []; o16 = []; he4 = []
    c122 = None; he42 = None; m2 = None
    for ii in range(model[0][3]):
        c121 = model[ii + 1][c12Indx]*12
        o161 = model[ii + 1][o16Indx]*16
        he41 = model[ii + 1][he4Indx]*4
        m1 = model[ii + 1][0]
        if c122 is None:
            c122 = c121; o162 = o161; he42 = he41; m2 = m1
            continue
        
        c12.append((c121 + c122)*0.5)
        o16.append((o161 + o162)*0.5)
        he4.append((he41 + he42)*0.5)
        mass.append((m1 + m2)*0.5)
        c122 = c121; o162 = o161; m2 = m1; he42 = he41
    
    # Now look for intershell value
    maxHe = max(he4); avgC12 = 0; avgO16 = 0; totMass = 0
    for ii in range(len(c12)):
        if abs(he4[ii] - maxHe) < 0.4*maxHe:
            if he4[ii] > c12[ii] and c12[ii] > 0.1:
                avgC12 += c12[ii]*mass[ii]
                avgO16 += o16[ii]*mass[ii]
                totMass += mass[ii]
    
    try:
        avgC12 /= totMass
        avgO16 /= totMass
    except ZeroDivisionError:
        pass
    except:
        raise
    
    return avgC12, avgO16

def getTDUMass(binObj, speciesDict, massTh):
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
    maxC12 = 0
    maxO16 = 0
    prevModNum = 0
    
    # Read model
    while True:
        newModel = binObj.nextModel()
        if not newModel:
            break
        
        model = binObj.model
        
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
            if coreMass < minCoreMass:
                minCoreMass = coreMass
                totMass = model[0][1]
            
        else:
            inTDU = False
        
        # Get intershell c12 and o16
        intC12, intO16 = getIntershellC12O16(model, coreMass, speciesDict)
        if intC12 is None:
            raise Exception("C12 value not found!")
        if intO16 is None:
            raise Exception("O16 value not found!")
        
        # Get maximum intershell c12 and o16
        if intC12 > maxC12:
            maxC12 = intC12
        if intO16 > maxO16:
            maxO16 = intO16
        
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
            maxTempHyd, maxC12, maxO16)
    return val

def main():
    '''Main program'''
    
    if len(sys.argv) < 2:
        print "Usage: python {} <SnuppatOutputBIN>".format(sys.argv[0])
        return 1
    
    # Get mass threshold
    if len(sys.argv) > 2:
        massTh= float(sys.argv[2])
    else:
        massTh= 4e-5
    
    species = "../../data/species.dat"
    speciesDict = getSpeciesDict(species, 4)
    
    # Create object
    binObj = readBinaryModels(sys.argv[1])
    
    prevMinCore = None
    # Look for max and min core mass
    print "Lambda, Temperature, BCE Temperature, Hyd Temperature, Core Mass,",
    print "Envelope mass, Dredged Mass, C12, O16 ;",
    print "Mass Threshold = {}".format(massTh)
    while True:
        val = getTDUMass(binObj, speciesDict, massTh)
        maxCoreM, minCoreM, totMass, temp, tempBCE, tempHyd, c12, o16 = val
        if maxCoreM == 0 or totMass is None:
            break
        
        if prevMinCore is None:
            prevMinCore = minCoreM
            continue
        
        coreGrowth = maxCoreM - prevMinCore
        dredg = maxCoreM - minCoreM
        lambd = dredg/coreGrowth
        
        prevMinCore = minCoreM
        print lambd, temp, tempBCE, tempHyd, minCoreM,
        print totMass - minCoreM, dredg, c12, o16

if __name__ == "__main__":
    main()
