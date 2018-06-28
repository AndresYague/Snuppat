'''Cleans printout.out from half written models. This usually happens
if the integration was forcibly halted.'''

import os, sys, math

class ReadWriteObject(object):
    def __init__(self, readFile, writeFile):
        '''Initalize object'''
        
        self.readFile = readFile
        self.writeFile = writeFile
        self.fread = open(readFile, "r")
        self.fwrite = open(writeFile, "w")
    
    def writeWhenFind(self, string):
        '''
        Search document for <string>. Each time it is found, write
        everything it has found from the last <string> ocurrence or
        the beginning, and start storing the lines again. If, at some
        point, <string> is not found before the end of the document,
        the stored information is discarded.
        '''
        while True:
            # Start storing lines
            lines = ""; line = ""; found = False
            while True:
                lines += line
                line = self.fread.readline()
                
                if string in line:
                    found = True
                    lines += line
                    break
                elif len(line) == 0:
                    break
            
            if found:
                self.fwrite.write(lines)
            else:
                break
    
    def renameAndClose(self):
        '''
        Close the documents and rename writeFile as
        readFile.
        '''
        self.close()
        os.rename(self.writeFile, self.readFile)
    
    def close(self):
        '''Close the documents'''
        self.fread.close()
        self.fwrite.close()

def main():
    # Aton output file
    inFilePath = os.path.join("printout.out")
    
    # Temporal file for operations
    tempFile = os.path.join("temp.dat")
    
    # Size warning operations
    conversion = {0: " B", 1: " KiB", 2: " MiB", 3: " GiB", 4: " TiB"}
    flsiz = float(os.path.getsize(inFilePath))
    powr = int(math.log(flsiz, 1024))
    warn = str(flsiz/(1024**powr)) + conversion[powr]
    warn = "Warning, you need at most {} of free disk space to run this"\
            .format(warn) + " program. Press any key to continue or Ctrl" + \
            " + C to abort.\n"
    
    # Issue warning
    try:
        raw_input(warn)
    except KeyboardInterrupt:
        return 1
    except:
        raise
    
    # Create file object, clean it and rename it
    files = ReadWriteObject(inFilePath, tempFile)
    files.writeWhenFind("MODELLO")
    files.renameAndClose()

if __name__ == "__main__":
    main()
