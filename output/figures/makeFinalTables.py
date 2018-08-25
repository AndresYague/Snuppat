def getKeyList(indx, lst):
    '''Return adecuate value from list'''
    
    if indx >= 0 and indx < len(lst):
        return lst[indx]
    else:
        return "--"
    
def printTable(storeNamVal):
    '''Print table in order'''
    
    nCol = 3
    keys = storeNamVal.keys(); keys.sort()
    nEls = len(keys); div = nEls/nCol
    
    # Get number of lines for tables
    nlines = div if nEls % nCol == 0 else div + 1
    
    firstOfSecond = None
    for ii in range(nEls):
        zz1 = getKeyList(ii, keys)
        nam1, val1 = storeNamVal.get(zz1, ("--", "--"))
        
        zz2 = getKeyList(ii + nlines, keys)
        nam2, val2 = storeNamVal.get(zz2, ("--", "--"))
        
        zz3 = getKeyList(ii + nlines*2, keys)
        nam3, val3 = storeNamVal.get(zz3, ("--", "--"))
        
        if firstOfSecond is None:
            firstOfSecond = zz2
        elif zz1 == firstOfSecond:
            break
        
        print "{} & {} & {:5.2f} & ".format(nam1, zz1, float(val1)),
        print "{} & {} & {:5.2f} & ".format(nam2, zz2, float(val2)),
        
        if val3 != "--":
            print "{} & {} & {:5.2f}\\\\".format(nam3, zz3, float(val3))
        else:
            print "{} & {} & {}\\\\".format(nam3, zz3, val3)

def main():
    '''Transform plottedValues.dat into .tex tables'''
    
    arch = "plottedValues.dat"
    data = "../../data/species.dat"
    
    # Index zz and names
    zToName = {}
    with open(data, "r") as fread:
        for line in fread:
            lnlst = line.split()
            zz = int(lnlst[0]) - int(lnlst[2])
            name = lnlst[1]
            name = name[0].upper() + name[1:]
            
            zToName[zz] = name
    
    # Create and print tables
    storeNamVal = {}
    with open(arch, "r") as fread:
        for line in fread:
            if "#" in line:
                if len(storeNamVal) > 0:
                    printTable(storeNamVal)
                    print
                    storeNamVal = {}
                
                print line
                continue
            
            lnlst = line.split()
            if len(lnlst) == 0:
                continue
            
            zz = int(lnlst[0])
            name = zToName[zz]
            val = lnlst[1]
            storeNamVal[zz] = (name, val)
    
    if len(storeNamVal) > 0:
        printTable(storeNamVal)

if __name__ == "__main__":
    main()
