import os, math

def main():
    '''Program to format the solar abundances into relative Nx/NH abundances,
    further dividing them by isotopes'''
    
    outpt = os.path.join("data", "fmt_abund.dat")
    temp = os.path.join("data", "temp.dat")
    
    # Ask for mode
    answer = raw_input("Prepare netsu? (y/n)").lower()
    
    if ("n" in answer):
        print "Preparing nist ratios"
        abundncs = os.path.join("data", "solar_abundances.dat")
        isotrat = os.path.join("data", "isotopic_ratio.dat")
        
        # Open temporal file
        count = 0
        with open(temp, "w") as fwrite:
            with open(abundncs, "r") as fabund:
                # Perform for every element
                for line in fabund:
                    lnlst = line.split()
                    elem = lnlst[0].lower()
                    abund = 10**(float(lnlst[1]) - 12)
                    
                    # Now look for the isotopes
                    with open(isotrat, "r") as fisot:
                        found = False
                        
                        for line in fisot:
                            lnlst = line.split()
                            elem2 = lnlst[0].lower()
                            weight = int(lnlst[1])
                            relativ = float(lnlst[2])
                            
                            if (elem2 == elem):
                                found = True
                                
                                ewrit = elem
                                if (weight == 1) and (elem2 == "h"):
                                    ewrit = "p"
                                    
                                elif (weight == 2) and (elem2 == "h"):
                                    ewrit = "d"
                                
                                st = "{} {} {}".format(weight, ewrit, \
                                                       abund*relativ)
                                
                                fwrite.write(st + "\n")
                                count += 1
                                
                            elif (found):
                                break
        
    else:
        print "Preparing netsu ratios"
        abundncs = os.path.join("data", "initial_comp320AGS2009.dat")
        
        # Open temporal file
        count = 0; factor = None
        with open(temp, "w") as fwrite:
            with open(abundncs, "r") as fabund:
                for line in fabund:
                    lnlst = line.split()
                    
                    # Get factor
                    if factor is None:
                        factor = 0.999885/float(lnlst[1])
                    
                    fwrite.write("{} {} {}\n".format(lnlst[2], lnlst[0],\
                                                     float(lnlst[1])*factor))
                    
                    count += 1
    
    # Now write count and the values to output
    with open(outpt, "w") as fwrite:
        fwrite.write("{}\n".format(count))
        with open(temp, "r") as fread:
            for line in fread:
                fwrite.write(line)
    
    os.remove(temp)

if __name__ == "__main__":
    main()
