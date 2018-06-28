import sys, os

def getNextModel(openFile):
    store = []
    readModel = False
    while True:
        line = openFile.readline()
        
        if len(line) == 0:
            if readModel:
                break
            else:
                return None
        
        if "Model" in line:
            if not readModel:
                readModel = True
            else:
                openFile.seek(-len(line), 1)
                break
        
        store.append(line)
    
    return store

def main():
    '''Change the 3 column format from pgfortran to the expected one'''
    
    if len(sys.argv) < 2:
        print "Use: python {} <input>".format(sys.argv[0])
        return 1
    
    arch = sys.argv[1]
    lastModel = None
    with open(arch, "r") as fread:
        with open("out.dat", "w") as fwrite:
            while True:
                
                store = getNextModel(fread)
                
                # If in last model, exit
                if store is None:
                    break
                
                # Get model number
                lnlst = store[0].split()
                modelNumb = int(lnlst[-1])
                if lastModel is not None and lastModel > modelNumb:
                    continue
                else:
                    lastModel = modelNumb
                
                # Now fix model
                oneLine = ''
                for line in store:
                    lnlst = line.split()
                    
                    # Join information
                    if '#' not in line and len(lnlst) > 0:
                        oneLine += ' '.join(lnlst) + ' '
                        
                        if len(lnlst) <= 2:
                            # Write line
                            fwrite.write(oneLine.strip() + '\n')
                            oneLine = ''
                    else:
                        # Just write the comments
                        fwrite.write(line)
    
    # Rename
    os.rename("out.dat", arch)

if __name__=="__main__":
    main()
