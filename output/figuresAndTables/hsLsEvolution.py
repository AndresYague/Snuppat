import sys, math, numpy, struct
import matplotlib.pyplot as plt

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
        if self.head is None:
            return False
        
        # Skip file
        for ii in range(head[3]):
            for jj in range(head[4]):
                self.fread.read(8)
        
        return True

def main():
    '''Get evolution of hs/ls vs s in [X/Fe]'''
    
    # Check arguments
    if len(sys.argv) < 2:
        print("Usage python {} <model>".format(sys.argv[0]))
        return 1
    
    data = "../../data/species.dat"
    archivo = sys.argv[1]
    hsLsElems = ["sr", "y", "zr", "ba", "la", "ce"]
    
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
    modelObj = readBinaryModels(archivo)
    
    # Each line has mass, temperature, rho, radiat
    # and elements in number fraction
    evolXFe = []; jj = -1
    while True:
        isNewModel = modelObj.nextModel()
        if not isNewModel:
            break
        
        header = modelObj.head
        model = modelObj.model
        
        # Report some progress
        print(jj)
        
        # Find the surface for this model
        for ii in range(1, len(model)):
            mass = (model[ii - 1][0] + model[ii][0])*0.5
            
            # If found surface, extract information
            if mass >= 0.85:
                
                prevLine = model[ii - 1]
                newLine = model[ii]
                
                # Take all abundances
                dens = [(x + y)*0.5 for (x, y) in zip(prevLine[4:], newLine[4:])]
                
                xFeVals = {}
                # Add the values for each element
                for ii in range(len(atomicNum)):
                    key = atomicNum[ii]
                    xFeVals[key] = xFeVals.get(key, 0) + dens[ii]*atomicMass[ii]
                
                # Now calculate values of interest
                selectedFe = []
                feVal = xFeVals[namesZ["fe"]]
                sunFeVal = solarValues[namesZ["fe"]]
                for elem in hsLsElems:
                    try:
                        val = xFeVals[namesZ[elem]]/feVal + 1e-100
                    except KeyError:
                        print("{} is not on the list".format(elem))
                    except:
                        raise
                    sunVal = solarValues.get(namesZ[elem], 1e-100)/sunFeVal
                    val = math.log10(val) - math.log10(sunVal)
                    selectedFe.append(val)
                
                break
        
        evolXFe.append(selectedFe)
    
    # Calculate hs/Fe, ls/Fe and s/Fe
    sFe = []; hsLs = []
    for arr in evolXFe:
        lsFe = sum(arr[0:3])/3.
        hsFe = sum(arr[3:])/3.
        
        sFe.append(sum(arr)/6.)
        hsLs.append(hsFe - lsFe)
    
    # Now plot values
    plt.plot(sFe, hsLs, lw = 2)
    
    plt.xlabel("[s/Fe]")
    plt.ylabel("[hs/ls]")
    
    print("# file: {}".format(archivo))
    print("# sFe hsLs")
    for ii in range(len(sFe)):
        print("{} {}".format(sFe[ii], hsLs[ii]))
    
    plt.show()
    
    return 0

if __name__ == "__main__":
    main()
