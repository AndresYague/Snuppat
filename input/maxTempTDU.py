import sys

def main():
    '''Get maximum temperature per model from physics file'''
    
    if len(sys.argv) < 3:
        print("Use: python {} <maxTemps> <table>".format(sys.argv[0]))
        return 1
    
    fil1 = sys.argv[1]
    fil2 = sys.argv[2]
    
    # Open both files and get maximum temp for TDU
    with open(fil1, "r") as fread1:
        with open(fil2, "r") as fread2:
            # Skip first line of fread2
            fread2.readline()
            
            while True:
                # Get age range
                lnlst2 = fread2.readline().split()
                if len(lnlst2) == 0:
                    break
                
                rang = (float(lnlst2[1]), float(lnlst2[3]))
                
                # Read temperatures
                maxTemp = 0
                while True:
                    lnlst1 = fread1.readline().split()
                    if len(lnlst1) == 0:
                        break
                    
                    age = float(lnlst1[1])
                    if age < rang[0]:
                        continue
                    elif age > rang[1]:
                        break
                    
                    temp = float(lnlst1[2])
                    if temp > maxTemp:
                        maxTemp = temp
                
                print(lnlst2[0], lnlst2[1:4], maxTemp)

if __name__ == "__main__":
    main()
