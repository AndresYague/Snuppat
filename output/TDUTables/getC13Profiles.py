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
        
        self.model = []
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

class readConvVel(object):
    '''Class for reading binary models'''
    
    def __init__(self, fil):
        '''Initialize'''
        
        super(readConvVel, self).__init__()
        if fil is not None:
            self.fread = open(fil, "r")
        else:
            self.fread = None
        
        self.modNum = None
        self.velocities = None
    
    def close(self):
        '''Close file'''
        
        self.fread.close()
    
    def nextModel(self):
        '''Read next model'''
        
        if self.fread is None:
            return None
        
        # Read header
        while True:
            line = self.fread.readline()
            if len(line) == 0:
                return False
            
            if "Model" in line:
                break
        
        velocities = []
        modNum = int(line.split()[-1])
        while True:
            line = self.fread.readline()
            if "Model" in line:
                self.fread.seek(-len(line), 1)
                break
            
            lnlst = line.split()
            mass = float(lnlst[-2])
            vel = float(lnlst[-1])
            
            velocities.append((mass, vel))
        
        self.modNum = modNum
        self.velocities = velocities
        
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

def getCoreMass(model, head, speciesDict):
    '''Return the hydrogen depleted core mass of this model, defined
    as the point where H < 0.1*maxH'''
    
    hindx = speciesDict["p1"]
    mass = []; hydro = []
    h2 = None; m2 = None
    for ii in range(head[3]):
        h1 = model[ii][hindx]
        m1 = model[ii][0]
        if h2 is None:
            h2 = h1; m2 = m1
            continue
        
        hydro.append((h1 + h2)*0.5)
        mass.append((m1 + m2)*0.5)
        h2 = h1; m2 = m1
    
    maxHydr = max(hydro)
    ii = len(hydro) - 1
    while ii > -1:
        if hydro[ii] <= 1e-1*maxHydr:
            return mass[ii]*head[1]
        
        ii -= 1
    
    return None

def getCCoreMass(model, head, speciesDict):
    '''Return the carbon core mass of this model, defined
    as the point where He4 < 1e-5'''
    
    he4indx = speciesDict["he4"]
    mass = []; helium = []
    he2 = None; m2 = None
    for ii in range(head[3]):
        he1 = model[ii][he4indx]
        m1 = model[ii][0]
        if he2 is None:
            he2 = he1; m2 = m1
            continue
        
        helium.append((he1 + he2)*2)
        mass.append((m1 + m2)*0.5)
        he2 = he1; m2 = m1
    
    ii = len(helium) - 1
    while ii > -1:
        if helium[ii] <= 1e-4:
            return mass[ii]*head[1]
        
        ii -= 1
    
    return None

def getPocket(model, head, coreMass, speciesDict):
    '''Return the effective c13 pocket'''
    
    # Store values
    c13Indx = speciesDict["c13"]
    n14Indx = speciesDict["n14"]
    p1Indx = speciesDict["p1"]
    mass = []; c13 = []; n14 = []; p1 = []
    c132 = None; n142 = None; p12 = None; m2 = None
    for ii in range(head[3]):
        c131 = model[ii][c13Indx]
        n141 = model[ii][n14Indx]
        p11 = model[ii][p1Indx]
        m1 = model[ii][0]
        if c132 is None:
            c132 = c131; n142 = n141; p12 = p11; m2 = m1
            continue
        
        c13.append((c131 + c132)*0.5)
        n14.append((n141 + n142)*0.5)
        p1.append((p11 + p12)*0.5)
        mass.append((m1 + m2)*0.5)
        c132 = c131; n142 = n141; m2 = m1
    
    # Now look for effective c13Pocket
    m0 = None; m1 = None; maxp1 = 0
    effPock = []; width = 0
    for ii in range(len(mass)):
        effC13 = 13*(c13[ii] - n14[ii])
        
        if effC13 > 1.e-3:
            effPock.append((mass[ii], effC13))
            
            if p1[ii] > maxp1:
                maxp1 = p1[ii]
            
            m1 = mass[ii]
            if m0 is None:
                m0 = mass[ii]
    
    if m0 is not None:
        width = m1 - m0
    
    return effPock, width, maxp1

def getMaximumTemperature(model, head):
    '''Return the maximum temperature'''
    
    temp = []
    t2 = None
    for ii in range(head[3]):
        t1 = model[ii][1]*1e3
        if t2 is None:
            t2 = t1
            continue
        
        temp.append((t1 + t2)*0.5)
        t2 = t1
    
    return max(temp)

def getTDUMass(binObj, convVelReader, speciesDict, massTh):
    '''Finds next core masses'''
    
    # For storage
    foundTDU = False
    inTDU = False
    totMass = None
    maxCoreMass = 0
    minCoreMass = 0
    maxVel = None
    minIntMass = None
    maxWidth = 0
    maxC13Pocket = None
    maxTPTemp = 0
    
    # Read model
    while True:
        isNewModel = binObj.nextModel()
        
        if not isNewModel:
            break
        
        head = binObj.head
        model = binObj.model
        
        # Get core mass
        coreMass = getCoreMass(model, head, speciesDict)
        cCoreMass = getCCoreMass(model, head, speciesDict)
        
        if coreMass is None:
            raise Exception("Core mass not found!")
        
        if coreMass > maxCoreMass and not inTDU:
            maxCoreMass = coreMass
            minCoreMass = coreMass
        elif (maxCoreMass - coreMass) > massTh:
            inTDU = True
            foundTDU = True
            modNum = head[0]
            
            # Look for modNum in velocities
            while True:
                newConvVel = convVelReader.nextModel()
                if not newConvVel:
                    break
                elif convVelReader.modNum == modNum:
                    break
            
            # Check velocities
            if newConvVel:
                vels = convVelReader.velocities
                minDiff = None
                for mass, vel in vels:
                    diff = abs(coreMass - mass)
                    if minDiff is None or minDiff > diff:
                        minDiff = diff
                        minVel = vel
                
                velBCE = minVel
                if maxVel is None or maxVel < velBCE:
                    maxVel = velBCE
            
            if coreMass < minCoreMass:
                minCoreMass = coreMass
                totMass = head[1]
            
            intMass = coreMass - cCoreMass
            if minIntMass is None or minIntMass > intMass:
                minIntMass = intMass
            
        else:
            inTDU = False
        
        # Get maximum temperature
        maxTemp = getMaximumTemperature(model, head)
        if maxTemp > maxTPTemp:
            maxTPTemp = maxTemp
        
        # Get effective c13 pocket if any
        c13Pocket, width, maxp1 = getPocket(model, head, coreMass, speciesDict)
        
        # Select maximum width effective c13 pocket
        if maxp1 < 1e-6 and width > maxWidth:
            maxC13Pocket = c13Pocket
            maxWidth = width
        
        if not inTDU and foundTDU:
            break
    
    return maxCoreMass, minCoreMass, minIntMass, maxC13Pocket, \
            maxTPTemp, maxVel, totMass

def main():
    '''Main program'''
    
    if len(sys.argv) < 2:
        print "Usage: python {} ".format(sys.argv[0]),
        print "<SnuppatOutputBIN> [convVel file]"
        return 1
    
    # Get mass threshold
    convVelFile = None
    if len(sys.argv) > 2:
        convVelFile = sys.argv[2]
    
    massTh= 4e-5
    
    species = "../../data/species.dat"
    speciesDict = getSpeciesDict(species, 4)
    
    # Create objects
    convVelReader = readConvVel(convVelFile)
    binObj = readBinaryModels(sys.argv[1])
    
    prevMinCore = None
    # Look for max and min core mass
    print "# Mass Threshold = {}".format(massTh)
    while True:
        val = getTDUMass(binObj, convVelReader, speciesDict, massTh)
        maxCoreMass, minCoreMass, intMass, maxPocket, maxTPTemp, convVel, totMass = val
        if maxCoreMass == 0 or totMass is None:
            break
        
        if prevMinCore is None:
            prevMinCore = minCoreMass
            continue
        
        coreGrowth = maxCoreMass - prevMinCore
        dredg = maxCoreMass - minCoreMass
        lambd = dredg/coreGrowth
        
        # Get pocket mass
        pockMass = 0
        if maxPocket is not None:
            pockMass = (maxPocket[-1][0] - maxPocket[0][0])*totMass
        
        prevMinCore = minCoreMass
        print "# Lambda = {}; Pocket Size = {}; ".format(lambd, pockMass),
        print "Core mass = {}; Intershell mass = {};".format(maxCoreMass, intMass),
        if convVel is not None:
            print "Convective velocity = {}".format(convVel),
        
        print ""
        if maxPocket is None:
            print "# --"
        else:
            for elem in maxPocket:
                print elem[0]*totMass, elem[1]
        
        print

if __name__ == "__main__":
    main()
