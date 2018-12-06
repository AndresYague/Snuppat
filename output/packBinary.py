import struct, sys

class readAsciiModels(object):
    '''Class for reading ascii models'''
    
    def __init__(self, fil):
        '''initialize'''
        
        super(readAsciiModels, self).__init__()
        self.fread = open(fil, "r")
    
    def close(self):
        '''Close file'''
        
        self.fread.close()
    
    def nextModel(self):
        '''Get next model'''
        
        # First get header
        try:
            nMod = int(self.fread.readline().split()[-1])
        except IndexError:
            return None
        except:
            raise
        
        mass = float(self.fread.readline().split()[-1])
        age = float(self.fread.readline().split()[-1])
        head = [nMod, mass, age]
        
        # Now read the model lines
        ii = 1
        model = [[]]
        
        line = self.fread.readline().split()
        nCols = len(line)
        while len(line) > 0:
            model.append([float(x) for x in line])
            line = self.fread.readline().split()
            ii += 1
        
        ii -= 1
        head += [ii, nCols]
        
        model[0] = head
        return model

def main():
    if len(sys.argv) < 2:
        print "Usage: python {} <input file>".format(sys.argv[0])
        return 1
    
    arch = sys.argv[1]
    models = readAsciiModels(arch)
    
    with open(arch + "BIN", "wb") as fwrite:
        while True:
            model = models.nextModel()
            if model is None:
                break
            
            # Pack head
            s = struct.pack('i', model[0][0])
            fwrite.write(s)
            s = struct.pack('d', model[0][1])
            fwrite.write(s)
            s = struct.pack('d', model[0][2])
            fwrite.write(s)
            s = struct.pack('i', model[0][3])
            fwrite.write(s)
            s = struct.pack('i', model[0][4])
            fwrite.write(s)
            
            # Now pack the rest
            for ii in range(model[0][3]):
                for jj in range(model[0][4]):
                    s = struct.pack('d', model[ii + 1][jj])
                    fwrite.write(s)
            
            print model[0]
    
    models.close()

if __name__ == "__main__":
    main()
