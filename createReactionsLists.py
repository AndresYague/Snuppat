'''Creates reactions list using data/species.dat and the data/xay.dat. It then
writes into the data/xay.lst and partition.lst'''

import os

class Element(object):
    '''Element(name[, weight = '' [, count = "no"]])
    
    Class Element(object)
    Defined method x.__eq__(y)
    DEfined attributes self.name and self.weight'''
    
    count = 0
    
    def __init__(self, name, weight = '', count = "no"):
        if count == "yes":
            Element.count += 1
            self.posit = str(Element.count)
        
        self.name = name + weight
        self.weight = weight
    
    def __eq__(self, other):
        if type(self) is type(other):
            if self.name == other.name:
                self.posit = other.posit
                return True
        
        return False
    
    def __str__(self):
        stri = "{} {}".format(self.posit, self.name)
        return stri

def change_name(name):
    '''Changes "n" for "n1", "p" for "p1", "d" for "d2",
    al-6 for al-26, al*26 for al*6, kr-5 for kr-85 and kr*5 for kr*85'''
    
    names = {"n":"n1", "p":"p1", "d":"d2", "al-6":"al-26", "al*6":"al*26",\
             "kr-5":"kr-85", "kr*5":"kr*85"}
    
    if name in names.keys():
        return names[name]
    
    return name
 
def isInList(nucleon, netwrk, nucleon_list = None):
    '''isInList(nucleon, network[, newlist = None])
    
    Take a nucleon and check if in netwrk. If a nucleon
    list is provided, a class Element with the nucleon is
    appended to it'''
    
    nucleon = change_name(nucleon)
    temp = Element(nucleon, '')
    
    if nucleon_list is not None:
        nucleon_list.append(temp)
    
    return (temp in netwrk)

def sortHighTempReactions():
    '''This script breaks the file data/Aton_Sp* into the data/xay.dat lists'''
    database = "Aton_Sprocess_compilation"
    database = os.path.join("data", database)
    
    files = ("1a1", "1a2", "1a3", "2a1", "2a2", "2a3", "2a4", "3a1", "3a2")
    
    nn = 0; file_open = False
    with open(database, "r") as fread:
        for line in fread:
            lnlst = line.split()
            
            if len(lnlst) == 1 and int(lnlst[0]) in range(len(files) + 1):
                nn = int(lnlst[0]) - 1
                if file_open:
                    fwrite.close()
                    file_open = False
                
            elif len(lnlst) > 1:
                if not file_open:
                    fwrite = open(os.path.join("data", files[nn] + ".dat"), "w")
                    file_open = True
                
                if line[0] == "0":
                    fwrite.write(line[1:])
                else:
                    fwrite.write(line)

def writePartition(refpartition, partition, netwrk):
    '''writePartition(refpartition, partition, netwrk)
    
    Function to write partition function. Argument
    refpartition is the file data/partition.lst,
    partition is the file data/partition.dat, and
    netwrk is the list of elements'''
    with open(refpartition, "w") as fwrite:
        for ele in netwrk:
            lines = 0
            with open(partition, "r") as fread:
                for line in fread:
                    lnlst = line.split()
                    lnlst[0] = change_name(lnlst[0])
                    
                    if lines == 0 and lnlst[0] == ele.name:
                        fwrite.write("{}\n".format(ele.posit))
                        lines = 1
                    elif lines > 0:
                        fwrite.write(line)
                        lines += 1
                        if lines > 3:
                            break

def highTemp():
    '''Extract and order the reaction fits'''
    
    lists = ("1a1", "1a2", "1a3", "2a1", "2a2", "2a3", "2a4", "3a1", "3a2")
    species = os.path.join("data", "species.dat")
    partition = os.path.join("data", "partition.dat")
    refpartition = os.path.join("data", "partition.lst")
    
    # Remove files
    for arch in lists:
        path = os.path.join("data", arch)
        
        try:
            os.remove(path + ".lst")
            
        except OSError as e:
            if e.strerror != "No such file or directory":
                raise
            
        except:
            raise
    
    
    # Create .dat lists
    sortHighTempReactions()
    
    # Create list with elements:
    netwrk = list()
    with open(species, "r") as fread:
        for line in fread:
            parts = line.split()
            netwrk.append(Element(parts[1], parts[0], count = "yes"))
    
    # Create partition.lst
    print("Creating partition.lst...")
    writePartition(refpartition, partition, netwrk)
    print("Finished partition functions")
    
    # List of elements that should be written at the beginning of the
    # first file with their index
    check = ("p", "d", "t", "he3", "he4", "n", "p", "he4", "c12", "c13", \
             "n14", "ne22")
    
    # Search every file for reactions that matches our network
    # and write it:
    completed = 0; reactions = 0
    for data in lists:
        if os.path.isfile(os.path.join("data", data + ".dat")):
            # Open files and start line counting
            count = 1
            fread = open(os.path.join("data", data + ".dat"), "r")
            fwrite = open(os.path.join("data", data + ".lst"), "w")
            
            # At the beginning of 1a1.lst write out how many lines more
            # should be skipped when reading the actual cross sections,
            # the number of elements, and some elements name and
            # index (the ones in the tuple "check")
            if data == "1a1":
                fwrite.write("{}\n".format(2*len(check) + 1))
                fwrite.write("{}".format(Element.count))
                
                templst = list()
                for elmnt in check:
                    if isInList(elmnt, netwrk, templst):
                        fwrite.write("\n{}\n{}".format(elmnt, \
                                                  templst[-1].posit))
                    else:
                        fwrite.write("\n{}\n0".format(elmnt))
                
                fwrite.write("\n")
            
            # Now fill the ".lst" files
            for line in fread:
                if count == 1:
                    parts = line.split()
                    nucleons = parts[0:-2]
                    
                    # Check that all reactions are
                    nuclList = list()
                    allin = [isInList(x, netwrk, nuclList) for x in nucleons]
                    reduced = True
                    for elem in allin:
                        if not elem:
                            reduced = False
                            break
                    allin = reduced
                    
                    if allin:
                        reactions += 1
                        for nucleon in nuclList:
                            fwrite.write("{:>4}".format(nucleon.posit))
                        
                        fwrite.write("{:>7}".format(parts[-2]))
                
                if count > 1 and allin:
                    for part in line.split():
                        fwrite.write("  {}".format(part))
                    
                    if count == 3:
                        fwrite.write('\n')
                
                if count == 3:
                    count = 0
                count += 1
            
            fwrite.close()
            fread.close()
            
            # Remove the ".dat" file
            os.remove(os.path.join("data", data + ".dat"))
        
        completed += 1
        if completed > 1:
            # Move in shell: up one line, back three columns
            # and erase line
            print("[1A", end = " ")
            print("[30D", end = " ")
            print("[K", end = " ")
        
        # Write progress
        print("Done {}/{}".format(completed, len(lists)))
    
    print("{} elements in {} reactions".format(Element.count, reactions))
    
    # Restart count
    Element.count = 0

def lowTemp():
    '''Extract and order the reaction tables'''
    
    # Initialize
    files = ("1a1lowTemp", "1a2lowTemp", "1a3lowTemp")
    files += ("2a1lowTemp", "2a2lowTemp", "2a3lowTemp")
    files += ("2a4lowTemp", "3a1lowTemp", "3a2lowTemp")
    species = os.path.join("data", "species.dat")
    database = os.path.join("data", "starlib_v6.dat")
    
    # Check if in unix
    inUnix = True if os.name == "posix" else False
    
    # Get total number of lines
    if inUnix:
        totLines = os.popen("wc -l {}".format(database)).read()
        totLines = int(totLines.split()[0])
    else:
        totLines = None
    
    # Remove files
    for arch in files:
        path = os.path.join("data", arch)
        
        try:
            os.remove(path + ".lst")
            
        except OSError as e:
            if e.strerror != "No such file or directory":
                raise
            
        except:
            raise
    
    # Create list with elements:
    netwrk = list()
    with open(species, "r") as fread:
        for line in fread:
            lnlst = line.split()
            netwrk.append(Element(lnlst[1], lnlst[0], count = "yes"))
    
    # Create the files
    for arch in files:
        fwrite = open(os.path.join("data", arch + ".lst"), "w")
        fwrite.close()
    
    # Print initial non-progress
    print("Done {}%".format(0))
    
    # Search and store
    reactions = 0; currLine = 0; oldPrctg = 0
    with open(database, "r") as fread:
        num = None
        stri1 = ""; stri2 = ""; stri3 = ""
        
        for line in fread:
            currLine += 1
            lnlst = line.split()
            
            if len(lnlst) > 3:
                if num is not None:
                    # Append to file
                    fapnd = open(os.path.join("data", files[num] + ".lst"), "a")
                    
                    # Add to stri1 the length of the table and the newlines
                    stri1 += " {}\n".format(len(stri2.split()))
                    stri2 += "\n"
                    stri3 += "\n"
                    
                    # Write
                    fapnd.write(stri1)
                    fapnd.write(stri2)
                    fapnd.write(stri3)
                    
                    # Close
                    fapnd.close()
                
                num = int(lnlst[0]) - 1
                stri1 = ""; stri2 = ""; stri3 = ""
                
                # If not considering the kind, ignore it
                if num > len(files) - 1:
                    num = None
                else:
                    nucleons = lnlst[1:-2]
                    
                    # Check that all reactions are
                    nuclList = list()
                    allin = [isInList(x, netwrk, nuclList) for x in nucleons]
                    reduced = True
                    for elem in allin:
                        if not elem:
                            reduced = False
                            break
                    allin = reduced
                    
                    # Add to stri1 the isotope indices and the source
                    if allin:
                        reactions += 1
                        for nucleon in nuclList:
                            stri1 += "{:>4}".format(nucleon.posit)
                        
                        stri1 += "{:>7}".format(lnlst[-2])
                        
                    else:
                        num = None
                
            elif len(lnlst) == 3:
                if num is not None:
                    stri2 += lnlst[0] + " "
                    stri3 += lnlst[1] + " "
            
            # Print progression
            if totLines is not None and currLine%1000 == 0:
                prctg = float(currLine)/totLines
                prctg = int(prctg*100)
                
                if prctg > oldPrctg:
                    # Move in shell: up one line, back three columns
                    # and erase line
                    print("[1A", end = " ")
                    print("[30D", end = " ")
                    print("[K", end = " ")
                    
                    # Write precentage.
                    print("Done {}%".format(prctg))
                    oldPrctg = prctg
        
        # Check if write last reaction
        if num is not None:
            # Append to file
            fapnd = open(os.path.join("data", files[num] + ".lst"), "a")
            
            # Add to stri1 the length of the table
            stri1 += " {}\n".format(len(stri2.split()))
            stri2 += "\n"
            stri3 += "\n"
            
            # Write
            fapnd.write(stri1)
            fapnd.write(stri2)
            fapnd.write(stri3)
            
            # Close
            fapnd.close()
    
    print("{} elements in {} reactions".format(Element.count, reactions))
    
    # Restart count
    Element.count = 0

if __name__ == "__main__":
    print("----> High temperature reactions")
    highTemp()
    print("--------------------------")
    print("----> Low temperature reactions")
    lowTemp()
