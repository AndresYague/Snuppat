import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main():
    if len(sys.argv) < 2:
        print("Usage: python {} profile1 <profile2 ...>".format(sys.argv[0]))
        return 1
    
    colors = ["b:", "g--", "y-"]
    labels = ["$\omega = 0.10$", "$\omega = 0.12$", "$\omega = 0.14$"]
    
    ii = 0
    for arch in sys.argv[1:]:
        # Read event
        fread = open(arch, "r")
        
        # Get xx, yy
        xxArr = []; yyArr = []
        for line in fread:
            if "#" in line:
                continue
            
            xx, yy = map(float, line.split())
            xxArr.append(xx)
            yyArr.append(yy)
        
        plt.plot(xxArr, yyArr, colors[ii], label = labels[ii], lw = 2)
        ii += 1
    
    # Observations values
    # (AW Cyg, S Sct, SS Vir, SZ Sgr, U Hya, V460 Cyg, Z Psc)
    # Obs values: Rb, Sr, Y, Zr, Ba, La, Ce, Nd, Sm
    xxObs = [37, 38, 39, 40, 56, 57, 58, 60, 62]
    yyErrs = [0.25, 0.20, 0.20, 0.20, 0.3, 0.40, 0.45, 0.40, 0.40]
    abiStars = [
               [0.2, 0.4, 0.3, 0.5, 0, 0.3, "-", 0.5, "-"],
               [0.5, 0.8, 0.7, 0.5, 0.2, 0.1, "-", 0.2, "-"],
               ["-", 0.0, 0.5, 0.4, 0.3, 0.4, "-", 0.3, "-"],
               [0.1, 0.4, 0.9, 0.8, 0.8, 0.9, "-", 0.9, 0.6],
               [0.5, 0.8, 1.3, 1.1, 1.1, 0.9, 0.6, 1.0, 0.7],
               [0.4, 0.5, 0.7, 0.8, 0.8, 0.7, "-", 0.8, 0.4],
               [0.6, 0.9, 1.0, 1.0, 1.0, 1.1, 0.6, 0.9, 0.8]
              ]
    
    # Plot observations with errors
    gray = (0.75, 0.75, 0.75)
    star_ii = 0
    for star in abiStars:
        # Calculate things
        hsFe = sum(star[1:4])/3.; hsErr = sum(yyErrs[1:4])/3.
        lsFe = sum(star[4:6])/2.; lsErr = sum(yyErrs[4:6])/2.
        
        hsLs = hsFe - lsFe; hsLsErr = hsErr + lsErr
        sFe = sum(star[1:6])/5.; sErr = sum(yyErrs[1:6])/5.
        
        # Set style
        if star_ii == 2 or star_ii == 3:
            col = "sr"
        elif star_ii == 0 or star_ii == 5:
            col = "vb"
        else:
            col = "^k"
        
        # Now plot it with error bars
        plt.errorbar(sFe, hsLs, fmt = col, xerr = sErr, yerr = hsLsErr)
        star_ii += 1
    
    plt.xlabel("[s/Fe]", size = 12)
    plt.ylabel("[hs/ls]", size = 12)
    plt.legend(loc = 0, ncol = 2, prop = {"size": 12})
    plt.show()

if __name__ == "__main__":
    main()
