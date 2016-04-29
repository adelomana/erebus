### this script creates the input files for MEGAN5 and runs them in batch mode

import os,sys

# 0. user defined variables
print 'defining user variables...'
numberOfFiles=64 # this should be 64

# 1. creating the runner file
for i in range(numberOfFiles):
    
    count=str(i+1).zfill(2)
    print 'creating runner file %s...'%count

    runnerFile='extractors/runner_%s.txt'%count
    with open(runnerFile,'w') as f:
        
        f.write("open file='/Volumes/omics4tb/alomana/projects/ornl/data/rmaFiles/nr/concatenated_%s.rma';\n"%count)
        f.write("collapse rank='Species';\n")
        f.write("select rank='Species';\n")
        f.write("set drawer=RectangularPhylogram;\n")
        f.write("show window=message;\n")
        f.write("list summary=all;\n")
        f.write("quit;\n")

    # 2. executing
    print 'executing runner file %s...'%count
    cmd='/Users/alomana/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E -c /Users/alomana/gDrive2/projects/ornl/src/%s > /Volumes/omics4tb/alomana/projects/ornl/data/hitsFiles/nr/sample_%s.txt'%(runnerFile,count)

    print
    print cmd
    print

    os.system(cmd)
    print
    print '%s job completed.'%count
    
