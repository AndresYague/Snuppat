import sys

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
    if len(sys.argv) < 3:
        print "Usage {} <input> <output> [time in ky]".format(sys.argv[0])
        print "If time is 0, it will fetch the first model."
        print "If no time is given it will get the last model."
        return 1
    
    if len(sys.argv) == 4:
        targetAge = float(sys.argv[3])
    else:
        targetAge = None
    
    inpt = sys.argv[1]
    outpt = sys.argv[2]
    
    stillSearching =  True
    firstAge = None; oldStore = None
    with open(inpt, "r") as fread:
        while stillSearching:
            
            store = getNextModel(fread)
            
            # If in last model, exit
            if store is None:
                if oldStore is not None:
                    store = oldStore
                    stillSearching = False
                    continue
                else:
                    return Exception("File {} is empty".format(inpt))
            
            # Check age
            if targetAge is not None:
                # Get current age
                for line in store:
                    if "Age" in line:
                        lst = line.split()
                        thisAge = 10**(float(lst[-1]) - 3) # in ky
                        break
                
                # Store first age or check if bigger
                if firstAge is None:
                    firstAge = thisAge
                    
                    # If age is 0, return first model
                    if targetAge == 0:
                        stillSearching = False
                    
                elif (thisAge - firstAge) >= targetAge:
                    stillSearching = False
            
            # Update oldStore for the event we reach the final model
            oldStore = store
    
    with open(outpt, "w") as fwrite:
        for line in store:
            fwrite.write(line)

if __name__ == "__main__":
    main()
