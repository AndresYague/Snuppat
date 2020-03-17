'''Remove repeated models'''

import os, sys, math

class ReadWriteObject(object):
    def __init__(self, readFile, writeFile):
        '''Initalize object'''
        
        self.readFile = readFile
        self.writeFile = writeFile
        self.fread = open(readFile, "r")
        self.fwrite = open(writeFile, "w")
        self.repeatedFis = {}
        self.repeatedChim = {}
    
    def searchRemoveRepeated(self):
        '''
        Search document for repeated models and store them. Then pass
        a second time and write the last copy of every model.
        '''
        
        # Check if in unix
        inUnix = True if os.name == "posix" else False
        
        if inUnix:
            totModels = os.popen("tail -n 1 {}".format(self.readFile)).read()
            totModels = int(totModels.split()[1])
            iniModel = None
            
        else:
            totModels = 0
        
        # Save useful strings
        pattMod = "MODELLO"
        pattFis = "FISICA"
        pattChim = "CHIMICA"
        
        # Create a dictionary where the keys are all the models found and the
        # value is the number of repetitions of said model.
        oldPrctg = 0
        print("Searching for repeated models")
        print()
        for line in self.fread:
            if pattMod in line:
                modNum = int(line.split()[1])
                
                if pattFis in line:
                    if self.repeatedFis.get(modNum) is None:
                        self.repeatedFis[modNum] = 0
                    else:
                        self.repeatedFis[modNum] += 1
                    
                elif pattChim in line:
                    if self.repeatedChim.get(modNum) is None:
                        self.repeatedChim[modNum] = 0
                    else:
                        self.repeatedChim[modNum] += 1
                
                # Write percentage to completion
                if iniModel is None:
                    iniModel = modNum
                
                prctg = float(modNum - iniModel)/\
                             (totModels - iniModel)
                prctg *= 100
                
                
                if inUnix and ((prctg - oldPrctg) > 1):
                    # Move in shell: up one line, back three columns
                    # and erase line
                    print("[1A", end = " ")
                    print("[30D", end = " ")
                    print("[K", end = " ")
                    
                    # Write precentage.
                    print("Done {}%".format(int(prctg)))
                    oldPrctg = prctg
        
        # Rewind fread and write models to fwrite if not repeated
        self.fread.seek(0)
        
        print("Copying models")
        print()
        oldPrctg = 0; lines = ""
        for line in self.fread:
            lines += line
            
            if pattMod in line:
                modNum = int(line.split()[1])
                
                # Check if repeated any of them. Write if not.
                if pattFis in line:
                    if self.repeatedFis[modNum] > 0:
                        self.repeatedFis[modNum] -= 1
                    else:
                        self.fwrite.write(lines)
                    
                elif pattChim in line:
                    if self.repeatedChim[modNum] > 0:
                        self.repeatedChim[modNum] -= 1
                    else:
                        self.fwrite.write(lines)
                
                # Clean lines
                lines = ""
                
                # Write percentage to completion
                if iniModel is None:
                    iniModel = modNum
                    oldPrctg = 0
                
                prctg = float(modNum - iniModel)/\
                             (totModels - iniModel)
                prctg *= 100
                
                
                if inUnix and ((prctg - oldPrctg) > 1):
                    # Move in shell: up one line, back three columns
                    # and erase line
                    print("[1A", end = " ")
                    print("[30D", end = " ")
                    print("[K", end = " ")
                    
                    # Write precentage.
                    print("Done {}%".format(int(prctg)))
                    oldPrctg = prctg
    
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
    files.searchRemoveRepeated()
    files.renameAndClose()

if __name__ == "__main__":
    main()
