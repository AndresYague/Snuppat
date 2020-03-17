'''Prepares the input for sporcess using printout.out in AtonOutput.
If not initial model is inputed it will take the first model'''

import os, sys, math

def main():
    nsteps = 2
    current = 1
    
    # Aton output file
    infilepath = os.path.join("AtonOutput", "printout.out")
    
    # Temporal file for operations
    tempfile = os.path.join("temp.dat")
    
    # Size warning operations
    conversion = {0: " B", 1: " KiB", 2: " MiB", 3: " GiB", 4: " TiB"}
    flsiz = float(os.path.getsize(infilepath))
    powr = int(math.log(flsiz, 1024))
    warn = str(flsiz/(2*(1024**powr))) + conversion[powr]
    warn = "Warning, you need at most {} of free disk space to run this"\
            .format(warn) + " program. Press any key to continue or Ctrl" + \
            " + C to abort.\n"
    
    # Issue warning
    try:
        raw_input(warn)
    except KeyboardInterrupt:
        return 1
    except:
        raise
    
    # Input management
    if len(sys.argv) > 1:
        try:
            ini_model = int(sys.argv[-1])
            
        except ValueError:
            print("Usage: python {} [model_number]".format(sys.argv[0]))
            return 2
            
        except:
            raise
        
    else:
        ini_model = None
    
    # First step: Extract chemistry from chosen model
    print("Extracting chemistry, step {} of {}".format(current, nsteps))
    current += 1
    
    chemst = os.path.join("input", "chemistry.dat")
    extracted = extract_chemistry(infilepath, chemst, ini_model)
    
    if not extracted:
        print("Model {} doesn't exist".format(ini_model))
        return 3
    
    # Second step: Extract physics from chosen model
    print("Extracting physics, step {} of {}".format(current, nsteps))
    print("---")
    current += 1
    
    phys = os.path.join("input", "physics.dat")
    minMaxDt = os.path.join("input", "minMaxDt.dat")
    extracted = extract_physics(infilepath, phys, minMaxDt, ini_model)
    
    if not extracted:
        print("Unable to extract physics models")
        return 4

def average(lst):
    return sum(lst)/float(len(lst))

def extract_chemistry(infilepath, chemst, ini_model):
    '''Extracts the chemistry from ini_model. If ini_model is None,
    then it will extract the first model'''
    
    mass_indx = 0
    with open(infilepath, "r") as fread:
        lines = list()
        nlines = 0
        for line in fread:
            if (len(line.split()) > 1) and (line[0] != '#'):
                lines.append(line) # Store each line except comments
                nlines += 1
            
            # If found end of chemistry model, check if we want it
            if ("MODELLO" in line) and ("CHIMICA" in line):
                nlines -= 1
                write_now = False
                
                if ini_model is not None:
                    # We found the desired model; write it
                    if int(line.split()[1]) == ini_model:
                        write_now = True
                        
                    # We assume the models are in increasing order
                    elif int(line.split()[1]) > ini_model:
                        return False
                        
                    else:
                        lines = list()
                        nlines = 0
                    
                # Here we want the first chemical model, so we'll just write it
                else:
                    write_now = True
                
                if write_now:
                    change_mass(lines, mass_indx)
                    with open(chemst, "w") as fwrite:
                        # First line must have the important information:
                        # model, mass, shell number
                        lnlst = lines[-1].split()
                        nelems = len(lines[0].split()) - 1
                        header = "{} {} {}\n".format(lnlst[1], lnlst[4], nlines)
                        fwrite.write(header)
                        
                        # Now everything else
                        for element in lines[:-1]:
                            fwrite.write(element)
                    
                    break
                
            elif "FISICA" in line:
                lines = list()
                nlines = 0
    
    return True

def extract_physics(infilepath, phys, minMaxDt, ini_model):
    '''Extracts the physics from ini_model to the end.
    If ini_model is None, then it will extract from the first model'''
    
    # Check if in unix
    in_unix = True if os.name == "posix" else False
    
    # Get total of models
    if in_unix:
        tot_models = os.popen("tail -n 1 {}".format(infilepath)).read()
        tot_models = int(tot_models.split()[1])
        
    else:
        tot_models = 0
    
    old_prctg = -1
    prctg = 0
    
    oldAge = None
    allDt = []
    
    mass_indx = 1
    found = True if ini_model is None else False
    
    with open(infilepath, "r") as fread:
        with open(phys, "w") as fwrite:
            lines = list()
            comments = list()
            nlines = 0
            for line in fread:
                if (len(line.split()) > 1) and (line[0] != '#'):
                    lines.append(line) # Store each line except comments
                    nlines += 1
                    
                elif (line[0] == '#'):
                    comments.append(line) # Store comments
                
                # If found end of physics model, check if we want it
                if ("MODELLO" in line) and ("FISICA" in line):
                    nlines -= 1
                    model = int(line.split()[1])
                    write_now = False
                    
                    if not found:
                        # We found the desired model; write it
                        if model == ini_model:
                            found = True
                            write_now = True
                            
                        # We assume the models are in increasing order
                        elif model > ini_model:
                            return False
                        
                    # Here we want every physical model onwards,
                    # so we'll just write them.
                    else:
                        write_now = True
                    
                    if write_now:
                        change_mass(lines, mass_indx)
                        
                        # Store age
                        for comment in comments:
                            if "Age" in comment:
                                cmnlst = comment.split()
                                age = cmnlst[cmnlst.index("Age=") + 1]
                                
                                # Calculate dt
                                if oldAge is not None:
                                    currDt = 10**float(age) - oldAge
                                    allDt.append(currDt)
                                    
                                    oldAge += currDt
                                else:
                                    oldAge = 10**float(age)
                        
                        # First line must have the important information:
                        # model, mass, age, shell number
                        lnlst = lines[-1].split()
                        fwrite.write("{} {} {} {}\n".format(lnlst[1],\
                                     lnlst[4], age, nlines))
                        
                        for element in lines[:-1]:
                            # Write M/Mtot, log T, log rho, radiat, radius, Hp,
                            # pressure, velocity, nuclear energy, 3alpha energy
                            lnlst = element.split()
                            fwrite.write("{} {} ".format(lnlst[1], lnlst[4]))
                            fwrite.write("{} {} ".format(lnlst[5], lnlst[12]))
                            fwrite.write("{} {} ".format(lnlst[2], lnlst[18]))
                            fwrite.write("{} {} ".format(lnlst[3], lnlst[20]))
                            fwrite.write("{} {}\n".format(lnlst[7], lnlst[8]))
                        
                        if ini_model is None:
                            ini_model = model
                    
                    lines = list()
                    comments = list()
                    nlines = 0
                    
                    if tot_models != ini_model:
                        prctg = float(model - ini_model)/\
                                     (tot_models - ini_model)
                        prctg *= 100
                        
                        if in_unix and ((prctg - old_prctg) > 1):
                            # Move in shell: up one line, back three columns
                            # and erase line
                            print("[1A", end = " ")
                            print("[30D", end = " ")
                            print("[K", end = " ")
                            
                            # Write precentage.
                            print("Done {}%".format(int(prctg)))
                            old_prctg = prctg
                    
                elif "CHIMICA" in line:
                    lines = list()
                    comments = list()
                    nlines = 0
    
    # Calculate min and max dt
    
    # Sort dt list
    allDt.sort()
    
    # Get average of first quarter:
    minDt = average(allDt[:len(allDt)/4])
    
    # Get average of last quarter:
    maxDt = average(allDt[int(len(allDt)*3./4.):])
    
    # Put min and max dt
    with open(minMaxDt, "w") as fwrite:
        fwrite.write("{} {}".format(minDt, maxDt))
    
    return True

def change_mass(lines, mass_indx):
    '''Program to change the mass coordinate to increase monotonically.
    The variable mass_indx holds the position on each line where the mass
    is located, starting at 0.'''
    
    ii = 0
    prev_mass = -1.
    while ii < len(lines):
        lnlst = lines[ii].split()
        
        # Comment, "MODELLO" or blank left alone
        if ("#" in lines[ii]) or ("MODELLO" in lnlst) or (len(lnlst) < 1):
            ii += 1
            continue
        
        if prev_mass <= float(lnlst[mass_indx]):
            prev_mass = float(lnlst[mass_indx])
            
        else:
            # Calculate new mass
            new_mass = 1 - float(lnlst[mass_indx])
            
            # Modify lnlst and copy
            lnlst[mass_indx] = str(new_mass)
            newline = " ".join(lnlst)
            lines[ii] = newline + "\n"
        
        ii += 1

if __name__ == "__main__":
    main()
