import sys

def main():
    '''Get maximum temperature per model from physics file'''
    
    if len(sys.argv) < 2:
        print("Use: python {} <physics file>".format(sys.argv[0]))
        return 1
    
    fil = sys.argv[1]
    with open(fil, "r") as fread:
        
        # Read first line
        firstAge = None
        lnlst = fread.readline().split()
        
        # Start cycle
        while True:
            nMod = int(lnlst[0])
            mass = float(lnlst[1])
            age = 10**(float(lnlst[2]) - 6)
            
            if firstAge is None:
                firstAge = age
            
            # Read model
            maxTemp = 0
            while True:
                lnlst = fread.readline().split()
                
                if len(lnlst) < 5:
                    break
                
                massCoord = float(lnlst[0])*mass
                tempCoord = 10**(float(lnlst[1]) - 6)
                if massCoord < 0.80:
                    continue
                
                if tempCoord > maxTemp:
                    maxTemp = tempCoord
            
            print(nMod, age - firstAge, maxTemp)
            
            # Stop if at the end of the document
            if len(lnlst) == 0:
                break

if __name__ == "__main__":
    main()
