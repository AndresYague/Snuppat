import sys

def main():
    '''Get maximum 3A per model from physics file'''
    
    if len(sys.argv) < 2:
        print("Use: python {} <physics file>".format(sys.argv[0]))
        return 1
    
    fil = sys.argv[1]
    with open(fil, "r") as fread:
        
        # Read first line
        firstAge = None
        lnlst = fread.readline().split()
        
        # Start cycle
        inPulse = False
        max3A, maxTemp, TPNum, maxSpecific = 0, 0, 0, 0
        while True:
            nMod = int(lnlst[0])
            mass = float(lnlst[1])
            age = 10**(float(lnlst[2]) - 6)
            
            if firstAge is None:
                firstAge = age
            
            # Read model
            integ3A = 0; integEpsNuc = 0
            prev3A, prevEpsNuc, prevMass = None, None, None
            while True:
                lnlst = fread.readline().split()
                
                if len(lnlst) < 5:
                    break
                
                massCoord = float(lnlst[0])*mass
                thisEpsNuc = float(lnlst[-2])
                this3A = float(lnlst[-1])
                thisTemp = 10**(float(lnlst[1]) - 6)
                
                # Maximum intershell temperature
                if this3A > 0 and thisTemp > maxTemp:
                    maxTemp = thisTemp
                
                # Maximum specific-shell 3A
                if this3A > maxSpecific:
                    maxSpecific = this3A
                
                if prev3A is not None:
                    # Trapezoidal rule
                    deltM = massCoord - prevMass
                    integ3A += deltM*(this3A + prev3A)*0.5
                    integEpsNuc += deltM*(thisEpsNuc + prevEpsNuc)*0.5
                
                prev3A = this3A
                prevEpsNuc = thisEpsNuc
                prevMass = massCoord
            
            if integ3A > max3A:
                max3A = integ3A
            
            if integ3A/integEpsNuc > 0.1:
                inPulse = True
            else:
                if inPulse:
                    TPNum += 1
                    print(TPNum, nMod, age - firstAge, max3A, maxSpecific, maxTemp)
                    max3A, maxTemp, maxSpecific = 0, 0, 0
                
                inPulse = False
            
            # Stop if at the end of the document
            if len(lnlst) == 0:
                break

if __name__ == "__main__":
    main()
