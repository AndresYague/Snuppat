import struct, sys
import matplotlib.pyplot as plt

class readModelsMonash(object):
    '''Class for reading Monash models'''
    
    def __init__(self, fil):
        '''Initialize'''
        
        super(readModelsMonash, self).__init__()
        self.fread = open(fil, "r")
    
    def close(self):
        '''Close file'''
        
        self.fread.close()
    
    def nextModel(self, currPos = None, onlyHead = False):
        '''Return next model'''
        
        if currPos is not None:
            self.fread.seek(currPos)
        
        # Look for the next model
        head = None; model = []; nMod = None
        while True:
            currPos = self.fread.tell()
            line = self.fread.readline()
            lnlst = line.split()
            
            # Exit if end of the file
            if len(line) == 0:
                break
            
            # Add new line to the model or exit if new model
            if nMod is not None:
                if len(lnlst) > 1:
                    model.append(map(float, line.split()))
                else:
                    break
                
            elif len(lnlst) == 1:
                nMod = int(lnlst[0])
                age = float(nMod) # TODO for now I don't have info on this
                head = [nMod, age]
                
                # Return now if only looking for head
                if onlyHead:
                    return head, currPos
        
        if nMod is None:
            return False
        else:
            self.model = [head + [len(model)], model]
            return True

def getDist(kk, storeCurr, storePrev):
    dist = 0; lenStore = len(storeCurr)
    for jj in range(lenStore):
        ii = lenStore - 1 + jj - kk
        if ii >= 0 and ii < lenStore:
            dist += abs(storeCurr[jj] - storePrev[ii])
            dist += abs(storePrev[jj] - storeCurr[ii])
        else:
            dist += abs(storeCurr[jj] + 1)
            dist += abs(storePrev[jj] + 1)
    
    return dist

def plotConve(ax, xxArr, convRegArr):
    '''Plot the convective regions'''
    
    # Make an array with one column per mass limit
    maxLen = 0
    gray = (0.5, 0.5, 0.5)
    for ii in range(len(xxArr)):
        for li in convRegArr[ii]:
            ax.plot([xxArr[ii], xxArr[ii]], li, color = gray, lw = 2)

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

def getCarbonMass(model, speciesDict, cMasses):
    '''Return carbon masses'''
    
    # Store all carbon masses
    cindex = speciesDict["c12"]
    mass = []; cc = []; temp = []
    c2 = None; m2 = None; t2 = None
    for ii in range(model[0][2]):
        c1 = model[1][ii][cindex]
        m1 = model[1][ii][0]
        if c2 is None:
            c2 = c1; m2 = m1
            continue
        
        cc.append((c1 + c2)*0.5)
        mass.append((m1 + m2)*0.5)
        c2 = c1; m2 = m1
    
    # Search for every carbon mass in cMasses from the surface to the core
    carbonMasses = None
    for cM in cMasses:
        ii = len(cc) - 1
        while ii > -1:
            if cc[ii] >= cM:
                if carbonMasses is None:
                    carbonMasses = []
                
                carbonMasses.append(mass[ii])
                break
            
            ii -= 1
        
        if ii == -1:
            carbonMasses.append(0)
    
    return carbonMasses

def getBCEMass(model, speciesDict):
    '''Return the BCE mass of this model, defined
    as the last convective poin where H ~ maxH'''
    
    hindx = speciesDict["p1"]
    mass = []; hydro = []; convV = []
    h2 = None; m2 = None; convV2 = None
    for ii in range(model[0][2]):
        h1 = model[1][ii][hindx]
        m1 = model[1][ii][0]
        convV1 = model[1][ii][3]
        if h2 is None:
            h2 = h1; m2 = m1; convV2 = convV1
            continue
        
        hydro.append((h1 + h2)*0.5)
        mass.append((m1 + m2)*0.5)
        convV.append((convV1 + convV2)*0.5)
        h2 = h1; m2 = m1; convV2 = convV1
    
    maxHydr = max(hydro)
    ii = len(hydro) - 1
    bceMass = None
    while ii > -1:
        if hydro[ii] >= (maxHydr*.9) and convV[ii] > 0:
            bceMass = mass[ii]
        
        ii -= 1
    
    return bceMass

def getConvRegions(model):
    '''Return convective regions'''
    
    mass = []; convV = []
    m2 = None; convV2 = None
    for ii in range(model[0][2]):
        m1 = model[1][ii][0]
        convV1 = model[1][ii][2]
        if m2 is None:
            m2 = m1; convV2 = convV1
            continue
        
        mass.append((m1 + m2)*0.5)
        convV.append((convV1 + convV2)*0.5)
        m2 = m1; convV2 = convV1
    
    conveZones = []
    inConve = False; m0 = None; m1 = None
    for ii in range(len(mass)):
        if convV[ii] > 0:
            if not inConve:
                inConve = True
                m0 = mass[ii]
            
        else:
            if inConve:
                m1 = mass[ii - 1]
                inConve = False
                if m1 - m0 > 1e-3:
                    conveZones.append((m0, m1))
    
    # Append last one if any
    if inConve:
        conveZones.append((m0, mass[-1]))
    
    return conveZones

def getValues(monashObj, speciesDict, cMasses):
    '''Gets core mass, envelope mass per model'''
    
    # Read model
    isNewMod = monashObj.nextModel()
    model = monashObj.model
    
    # Get core mass, BCE mass and convective regions
    carbonMass = getCarbonMass(model, speciesDict, cMasses)
    BCEMass = getBCEMass(model, speciesDict)
    convReg = getConvRegions(model)
    #age = 10**model[0][1] # TODO fix this at some point
    age = model[0][1]
    modNum = model[0][0]
    
    if carbonMass is None:
        raise Exception("Core mass not found!")
    if BCEMass is None:
        raise Exception("BCE mass not found!")
    
    return age, modNum, convReg, carbonMass, BCEMass, isNewMod

def main():
    '''Main program'''
    
    if len(sys.argv) < 2:
        print("Usage: python {} <MonashOutput>".format(sys.argv[0]))
        return 1
    
    species = "speciesMonash.dat"
    speciesDict = getSpeciesDict(species, 3)
    
    # Create object
    monashObj = readModelsMonash(sys.argv[1])
    
    # Define carbon masses
    cMasses = (0.3, 0.4, 0.5)
    
    firstAge = None
    ageArr = []; modNumArr = []; convRegArr = []
    carbonMassArr = []; envMassArr = []
    while True:
        vals = getValues(monashObj, speciesDict, cMasses)
        age, modNum, convReg, carbonMass, envMass, cont = vals
        if firstAge is None:
            firstAge = age
            for mass in carbonMass:
                carbonMassArr.append([])
            
            lenCarbon = len(carbonMass)
        
        ageArr.append(age - firstAge)
        modNumArr.append(modNum)
        convRegArr.append(convReg)
        for ii in range(lenCarbon):
            carbonMassArr[ii].append(carbonMass[ii])
        envMassArr.append(envMass)
        
        if not cont:
            break
        
        print(len(ageArr))
    
    xxArr = modNumArr
    #xxArr = ageArr
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plotConve(ax, xxArr, convRegArr)
    for ii in range(lenCarbon):
        ax.plot(xxArr, carbonMassArr[ii], lw = 2,
                label = "XC12 = {}".format(cMasses[ii]))
    ax.plot(xxArr, envMassArr, lw = 2, label = "H exhausted core")
    
    ax.set_xlabel("Model number")
    #ax.set_xlabel("Time in yrs")
    ax.set_ylabel("M/M$_\odot$")
    
    plt.legend(loc = 0)
    plt.show()

if __name__ == "__main__":
    main()
