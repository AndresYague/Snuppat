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
    '''Get evolution of one element in epsilon or [X/Fe]'''
    
    # Check arguments
    if len(sys.argv) < 4:
        print "Usage python {} <(eps|xfe|massf)> <model>".format(sys.argv[0]),
        print "<elem1> [elem2, elem3, ...]"
        return 1
    
    data = "../../data/species.dat"
    mode = sys.argv[1]
    archivo = sys.argv[2]
    elms = sys.argv[3:]
    
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
    ages = []; evolEps = []; evolXFe = []; evolMassF = []
    while True:
        isNewModel = modelObj.nextModel()
        if not isNewModel:
            break
        
        header = modelObj.head
        model = modelObj.model
        
        # Get the age
        age = 10**(header[2] - 3)
        if len(ages) == 0:
            ages.append(age)
        else:
            ages.append(age - ages[0])
        
        # Report some progress
        print len(ages)
        
        # Find the surface for this model
        for ii in range(1, len(model)):
            mass = (model[ii - 1][0] + model[ii][0])*0.5
            
            # If found surface, extract information
            if mass >= 0.85:
                
                prevLine = model[ii - 1]
                newLine = model[ii]
                
                # Take all abundances
                dens = [(x + y)*0.5 for (x, y) in zip(prevLine[4:], newLine[4:])]
                
                epsVals = {}; xFeVals = {}; mFVals = {}
                # Add the values for each element
                for ii in range(len(atomicNum)):
                    key = atomicNum[ii]
                    epsVals[key] = epsVals.get(key, 0) + dens[ii]
                    mFVals[key] = mFVals.get(key, 0) + dens[ii]*atomicMass[ii]
                    xFeVals[key] = mFVals[key]
                
                # Now calculate values of interest
                hydroVal = epsVals[namesZ["h"]]
                feVal = xFeVals[namesZ["fe"]]
                sunFeVal = solarValues[namesZ["fe"]]
                selectedEps = []; selectedFe = []; selectedMassF = []
                for elem in elms:
                    try:
                        val = epsVals[namesZ[elem]]/hydroVal + 1e-100
                    except KeyError:
                        print "{} is not on the list".format(elem)
                    except:
                        raise
                    
                    val = math.log10(val) + 12
                    selectedEps.append(val)
                    
                    try:
                        val = xFeVals[namesZ[elem]]/feVal + 1e-100
                    except KeyError:
                        print "{} is not on the list".format(elem)
                    except:
                        raise
                    sunVal = solarValues.get(namesZ[elem], 1e-100)/sunFeVal
                    val = math.log10(val) - math.log10(sunVal)
                    selectedFe.append(val)
                    
                    try:
                        val = mFVals[namesZ[elem]] + 1e-100
                    except KeyError:
                        print "{} is not on the list".format(elem)
                    except:
                        raise
                    selectedMassF.append(val)
                
                break
        
        evolEps.append(selectedEps)
        evolXFe.append(selectedFe)
        evolMassF.append(selectedMassF)
    
    # Transform age and evol values to something plottable
    ages[0] = 0
    evolEps = numpy.transpose(numpy.array(evolEps))
    evolXFe = numpy.transpose(numpy.array(evolXFe))
    evolMassF = numpy.transpose(numpy.array(evolMassF))
    
    # Now plot values
    if mode == "eps":
        for ii in range(len(elms)):
            evEps = evolEps[ii]
            minLen = min(len(ages), len(evEps))
            plt.plot(ages[0:minLen], evEps[0:minLen], label = elms[ii], lw = 2)
        
        plt.xlabel("TP-AGB time in ky")
        plt.ylabel("Log epsilon")
        plt.ylim([-2, 5])
        
    elif mode == "xfe":
        for ii in range(len(elms)):
            evXFe = evolXFe[ii]
            minLen = min(len(ages), len(evXFe))
            plt.plot(ages[0:minLen], evXFe[0:minLen], label = elms[ii], lw = 2)
        
        plt.xlabel("TP-AGB time in ky")
        plt.ylabel("[X/Fe]")
        plt.ylim([-0.2, 2])
        
    elif mode == "massf":
        for ii in range(len(elms)):
            evMassF = evolMassF[ii]
            minLen = min(len(ages), len(evMassF))
            plt.plot(ages[0:minLen], evMassF[0:minLen], label = elms[ii], lw = 2)
        
        plt.xlabel("TP-AGB time in ky")
        plt.ylabel("Mass fraction")
        plt.yscale("log")
    
    plt.legend(loc = 0)
    plt.show()
    
    return 0

if __name__ == "__main__":
    main()
