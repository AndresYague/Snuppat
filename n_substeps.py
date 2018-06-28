import sys

def main():
    '''Program to generate number of substeps for Bader-Deuflhard scheme'''
    # Generate first the number of substeps
    
    if len(sys.argv) < 2:
        print "Usage: {} <n of substeps>".format(sys.argv[0])
        return -1
    else:
        nsubs = int(sys.argv[1])
    
    # We have to generate the list 2*(2j + 1) where also
    # n_i/n_{i+1} <= 5/7
    
    nums = list(); j = 0
    while len(nums) < nsubs:
        newnum = 4*j + 2
        
        if len(nums) > 0:
            if nums[-1]/float(newnum) <= 5/7.:
                nums.append(newnum)
            
        else:
            nums.append(newnum)
        
        j += 1
    
    print "The substeps are: "
    print nums
    
    # Work coefficients for the bader-deuflhard method
    # cf := cost for f-evaluation
    # cj := cost for jacobian evaluation
    # cs := cost of backwards substitution
    # clr := cost of gaussian elimination
    # In this case we have: cf = 1, cj = 0, cs = 1, clr = 0
    # a1 = cj + clr + m_1*(cf + cs)
    # aj = a_(j-1) + m_j*(cf + cs) + clr
    cf, cj, cs, clr = 7, 0, 10, 0
    aj = list()
    for i in range(nsubs):
        if len(aj) > 0:
            aj.append(int(aj[-1] + (nums[i] + 1)*(cf + cs) + clr))
        else:
            aj.append(int(cj + clr + (nums[i] + 1)*(cf + cs)))
    
    print "The aj are: "
    print aj

if __name__ == "__main__":
    main()
