import sys

def main():
    if len(sys.argv) < 2:
        print("Give me an input file!")
    
    inpFil = sys.argv[1]
    dicYields = {}; isotopeList = []
    listOfWanted = ["hf182", "pb205", "cs135", "pd107"]
    with open(inpFil, "r") as fread:
        for line in fread:
            lnlst = line.split()
            if "#" in line or len(lnlst) == 0:
                continue
            
            isot = lnlst[0]
            if isot not in listOfWanted:
                continue
            
            zz = int(lnlst[1])
            aa = int(isot[2:]) - zz
            isotName= isot[0:2] + "-" + isot[2:]
            val = float(lnlst[2])
            if isotName not in dicYields:
                isotopeList.append((isotName, zz, aa))
                dicYields[isotName] = val
            else:
                dicYields[isotName] += val
    
    # Write the new file
    with open(inpFil + "Total", "w") as fwrite:
        for isot in isotopeList:
            isotName, zz, aa = isot
            fwrite.write("&{} &{:.2e} &{:d} &{:d}\n".format(isotName, dicYields[isotName], zz, aa))

if __name__ == "__main__":
    main()
