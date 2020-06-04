import sys, struct

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
        if self.head is None:
            return False
        
        # Skip file
        for ii in range(head[3]):
            for jj in range(head[4]):
                self.fread.read(8)
        
        return True

def main():
    if len(sys.argv) < 3:
        print("Usage {} <input> <output> [time in ky]".format(sys.argv[0]))
        print("If time is 0, it will fetch the first model.")
        print("If no time is given it will get the last model.")
        return 1
    
    if len(sys.argv) == 4:
        targetAge = float(sys.argv[3])
    else:
        targetAge = None
    
    inpt = sys.argv[1]
    outpt = sys.argv[2]
    
    stillSearching =  True
    firstAge = None; oldStore = None
    binObj = readBinaryModels(inpt)
    while stillSearching:
        
        gotNewModel = binObj.nextModel()
        
        # If in last model, exit
        if not gotNewModel:
            if oldStore is not None:
                binObj.model = oldStore
                stillSearching = False
                continue
            else:
                return Exception("File {} is empty".format(inpt))
        
        # Check age
        if targetAge is not None:
            # Get current age
            head = binObj.model[0]
            thisAge = 10**(head[2] - 3) # in ky
            
            # Store first age or check if bigger
            if firstAge is None:
                firstAge = thisAge
                
                # If age is 0, return first model
                if targetAge == 0:
                    stillSearching = False
                
            elif (thisAge - firstAge) >= targetAge:
                stillSearching = False
        
        # Update oldStore for the event we reach the final model
        oldStore = binObj.model
    
    with open(outpt, "w") as fwrite:
        # Write header first
        head = binObj.model[0]
        fwrite.write("# Model number = {}\n".format(head[0]))
        fwrite.write("# Mass = {}\n".format(head[1]))
        fwrite.write("# Age = {}\n".format(head[2]))
        
        for line in binObj.model[1:]:
            line = [str(x) for x in line]
            fwrite.write(" ".join(line) + "\n")

if __name__ == "__main__":
    main()
