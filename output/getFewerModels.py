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
    '''Reduce the number of models in an output'''
    
    if len(sys.argv) < 3:
        print("Use: python {} <input>".format(sys.argv[0]), end = " ")
        print("<model writting frequency>")
        return 1
    
    arch = sys.argv[1]
    skipMods = int(sys.argv[2])
    countModels = -1
    with open(arch, "r") as fread:
        with open("out.dat", "w") as fwrite:
            while True:
                
                store = getNextModel(fread)
                countModels += 1
                
                # If in last model, exit
                if store is None:
                    break
                
                # If avoid writing, skip
                if countModels%skipMods > 0:
                    continue
                
                # Write the model
                for line in store:
                    fwrite.write(line)
    
    # Rename
    os.rename("out.dat", arch)

if __name__=="__main__":
    main()
