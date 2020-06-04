import sys, os, struct

class readBinaryModels(object):
    '''Class for reading binary models'''
    
    def __init__(self, fil):
        '''Initialize'''
        
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

def getMassFracs(abunds, specs):
    '''Multiply number fraction abundances by mass for mass fraction'''
    
    massFracs = []
    for ii in range(len(specs)):
        nam, zz, mass = specs[ii]
        massFracs.append(abunds[ii]*mass)
    
    return massFracs

def main():
    '''Program to extract yields from the SNUPPAT models'''
    
    # Check input
    if len(sys.argv) < 3:
        print("Usage: python {}".format(sys.argv[0]),)
        print("<input file> <net yields? (y/n)>")
        return 1
        
    else:
        outFil = sys.argv[1]
        netYields = True if sys.argv[2] == "y" else False
    
    mode = "isotopes" # Default mode
    if len(sys.argv) >= 4:
        mode = sys.argv[3]
        if mode != "isotopes" and mode != "elements":
            print('The mode must be either "isotopes" or "elements"')
            return 1
    
    # Species dictionary
    specFile = os.path.join("..", "..", "data", "species.dat")
    specs = {}; idNam = {}; idZZ = {}; ii = 0
    with open(specFile) as fread:
        for line in fread:
            lnlst = line.split()
            
            zz = int(lnlst[0]) - int(lnlst[2])
            nam = lnlst[1]; mass = int(lnlst[0])
            
            if zz == 1:
                nam = "h"
            
            nam += lnlst[0]
            idNam[ii] = nam
            idZZ[ii] = zz
            specs[ii] = (nam, ii, mass)
            
            ii += 1
    
    totYields = None
    firstAge = None; prevMass = None
    prevMassFrac = None; initialAbund = None
    fread = readBinaryModels(outFil)
    # Read each model
    while fread.nextModel():
        age = 10**(fread.head[2] - 3) # Age in thousand years
        mass = fread.head[1]          # Mass in solar masses
        model = fread.model
        
        # Get abundances
        for ii in range(1, len(model) - 1):
            line = model[ii]
            nextLine = model[ii + 1]
            if (line[0] + nextLine[0])*0.5 <= 0.85:
                continue
            
            abundances = [(x + y)*0.5 for x, y in zip(line[4:], nextLine[4:])]
            break
        
        # If first model, store everything and skip
        if firstAge is None:
            if netYields:
                initMassFrac = getMassFracs(abundances, specs)
            else:
                initMassFrac = [0 for x in abundances]
            
            prevMassFrac = [x - y for x, y in
                    zip(getMassFracs(abundances, specs), initMassFrac)]
            firstAge = age
            prevMass = mass
            continue
        
        # Now calculate yields
        if totYields is None:
            totYields = [0 for x in abundances]
        
        # Integrate with trapezoidal rule
        massFracs = getMassFracs(abundances, specs)
        massFracs = [x - y for x, y in zip(massFracs, initMassFrac)]
        
        # Apply trapezoidal rule
        dMass = prevMass - mass
        yields = [(x + y)*0.5*dMass for x, y in zip(massFracs, prevMassFrac)]
        totYields = [x + y for x, y in zip(yields, totYields)]
        
        # Update all
        prevMassFrac = massFracs
        prevMass = mass
        
        # Correct age
        age -= firstAge
        
        # Print results
        if age > 0:
            print("# Mass: {} MSun. Time: {} ky".format(mass, age))
            print("Species, Current yield, Accumulated yields")
            for ii in range(len(idNam)):
                if mode == "elements":
                    formStr = "{:2} {:2} {:11.4E} {:11.4E}"
                elif mode == "isotopes":
                    formStr = "{:5} {:2} {:11.4E} {:11.4E}"
                
                s = formStr.format(idNam[ii], idZZ[ii], yields[ii], totYields[ii])
                print(s)
            
            print()

if __name__ == "__main__":
    main()
