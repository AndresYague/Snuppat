with open("physics3MSunZ018.dat", "r") as fread:
    maxTemp = None
    for line in fread:
        lnlst = line.split()
        if len(lnlst) == 4:
            if maxTemp is not None:
                print(modNum, 10**(maxTemp - 6))
            
            maxTemp = 0
            modNum = lnlst[0]
        else:
            temp = float(lnlst[1])
            if temp > maxTemp:
                maxTemp = temp

print(modNum, 10**(maxTemp - 6))
