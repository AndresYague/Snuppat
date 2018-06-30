import sys, os, copy

class outObj(object):
    '''Object for input'''
    
    def __init__(self, fil):
        self.fread = open(fil, "r")
        self.prevLine = None
    
    def getNextModel(self):
        '''Find next complete model. Return mass, age and the model itself'''
        
        if self.prevLine is not None:
            line = self.prevLine
            if "Model" not in self.prevLine:
                raise Exception("Not found model in prevLine")
            
        else:
            while True:
                line = self.fread.readline()
                if len(line) == 0:
                    break
                
                if "Model" in line:
                    break
        
        if len(line) == 0:
            return (None, None, None)
        
        mass = float(self.fread.readline().split()[-1])
        age = float(self.fread.readline().split()[-1])
        
        prevLine = None
        while True:
            line = self.fread.readline()
            if len(line) == 0:
                break
            
            if "Model" in line:
                self.prevLine = line
                break
            else:
                newLine = [float(x) for x in line.split()]
                
                # Store only if mass is above 0.95
                if newLine[0] > 0.95:
                    self.prevLine = None
                    break
                
                prevLine = newLine
        
        # Get values for store
        if prevLine is not None:
            store = [(x + y)*0.5 for x,y in zip(prevLine, newLine)]
        else:
            store = None
        
        return (mass, age, store)

def getSpecies(spec):
    '''Extract species information from file'''
    
    with open(spec, "r") as fread:
        orderElem = []
        namSpec = {}; ii = 0
        for line in fread:
            lnlst = line.split()
            
            name = lnlst[1]
            # Be careful with metastable names
            if "-" in name or "*" in name:
                name = name[0:-1]
            
            mass = int(lnlst[0])
            zz = mass - int(lnlst[2])
            tup = (name, zz, mass)
            order = (name, zz)
            
            # Add tuples to namSpec and orderElem
            if order not in orderElem:
                orderElem.append(order)
            namSpec[ii] = tup
            ii += 1
    
    return namSpec, orderElem

def transformAbund(li, spec):
    '''Program to change abundances to our preferred format'''
    
    liAbund = li[4:]
    newDi = {}
    for ii in range(len(spec)):
        key = spec[ii]
        nam, zz, mass = key
        newDi[(nam, zz)] = newDi.get((nam, zz), 0) + liAbund[ii]*mass
    
    return newDi

def main():
    '''Obtain the yield tables for selected model'''
    if len(sys.argv) < 2:
        print "Usage: python {} <outFile>".format(sys.argv[0])
        return 1
    
    # Store the indices corresponding to a specific atomic mass
    spec = os.path.join("..", "..", "data", "species.dat")
    namSpec, orderElem = getSpecies(spec)
    
    # Now open the output file and go model by model obtaining the values
    models = outObj(sys.argv[1])
    
    # Get first model the output of getNextModel is (mass, age, [line])
    firstModel = models.getNextModel()
    prevMass = firstModel[0]; firstAge = 10**firstModel[1]
    firstAbund = transformAbund(firstModel[2], namSpec)
    
    yields = copy.copy(firstAbund)
    for key in yields:
        yields[key] = 0
    
    
    # Now write results in the output
    outpt = os.path.split(sys.argv[1])[1]
    outpt += "Yields"
    with open(outpt, "w") as fwrite:
        while True:
            thisModel = models.getNextModel()
            if thisModel[-1] is None:
                break
            
            thisMass = thisModel[0]; thisAge = 10**thisModel[1]
            thisAbund = transformAbund(thisModel[2], namSpec)
            
            # Get dm and time up to now
            dm = prevMass - thisMass
            tt = (thisAge - firstAge)*1e-3
            
            for key in thisAbund:
                changeAbund = thisAbund[key] - firstAbund[key]
                yields[key] += (thisAbund[key] - firstAbund[key])*dm
            
            fwrite.write("# Mass (Msun), time (Ky): {} {}\n".format(thisMass,
                                                                    tt))
            
            for key in orderElem:
                yi = yields[key]
                fwrite.write("{} {} {}\n".format(key[0], key[1], yi))
            fwrite.write("\n")
            
            prevMass = thisMass

if __name__ == "__main__":
    main()
