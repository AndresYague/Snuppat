import sys

'''Program to extract a given model from output for tests or continuation'''

def searchModelChem(inpt, modNum):
    '''Takes an open file and returns a model'''
    
    modLen = -1
    s = ''; foundModel = False
    for line in inpt:
        # Look for the next model
        if "Model" in line:
            # If already found one, we can exit. If the model number is
            # greater than the one we are looking for, we can exit as well
            if foundModel:
                break
            elif int(line.split()[-1]) == modNum:
                foundModel = True
            elif int(line.split()[-1]) > modNum:
                break
        
        if foundModel:
            if "#" not in line:
                s += line
                modLen += 1
            if "Mass" in line:
                modMass = line.split()[-1]
    
    if not foundModel:
        return None, None, None
    
    return modMass, modLen, s

def searchModelPhys(inpt, modNum, leftOver):
    '''Takes an open file and returns a model'''
    
    s = ''; foundModel = False
    for line in inpt:
        # Look for the next model
        lnlst = line.split()
        if len(leftOver) > 0:
            leftLst = leftOver.split()
            
            if int(leftLst[0]) < modNum:
                leftOver = ''
            else:
                lnlst = leftLst
        
        if len(lnlst) == 4:
            # If already found one, we can exit. If the model number is
            # greater than the one we are looking for, we can exit as well
            if foundModel:
                leftOver = line
                break
            elif int(lnlst[0]) == modNum:
                foundModel = True
            elif int(lnlst[0]) > modNum:
                break
        
        if foundModel:
            s += leftOver
            s += line
            leftOver = ''
    
    if not foundModel:
        return None, None
    
    return s, leftOver

def main():
    '''Get input, open the files and look for the models.
    Then write them.'''
    
    if len(sys.argv) < 5:
        print("use: python {} <phys> <chem> model".format(sys.argv[0]), end = " ")
        print("<number of models>".format(sys.argv[0]))
        return 0
    
    phys = open(sys.argv[1], "r")
    chem = open(sys.argv[2], "r")
    modNum = int(sys.argv[3])
    nn = int(sys.argv[4])
    
    # Search for chemical model and exit if not found
    modMass, modLen, modChem = searchModelChem(chem, modNum)
    if modChem is None:
        print("Model unavailable")
        return 0
    
    # Open and write in contChem
    fwrite = open("contChem{}.dat".format(modNum), "w")
    fwrite.write("{} {} {}\n".format(modNum, modMass, modLen))
    for line in modChem.split('\n'):
        # Split line
        lnlst = line.split()
        if len(lnlst) == 0:
            continue
        
        fwrite.write(lnlst[0] + ' ')
        fwrite.write(' '.join(lnlst[4:]) + '\n')
    fwrite.close()
    
    # Get at least n consecutive physical models
    leftOver = ''
    fwrite = open("physics{}plus{}.dat".format(modNum, nn), "w")
    for i in range(nn):
        # Search for next model
        modPhys, leftOver = searchModelPhys(phys, modNum + i, leftOver)
        if modPhys is None:
            print("Model {} not found".format(modNum + i))
            break
        
        # Write it
        fwrite.write(modPhys)
    fwrite.close()
    
    phys.close()
    chem.close()

if __name__ == "__main__":
    main()
