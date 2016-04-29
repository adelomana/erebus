### this script creates the input files for MEGAN5 and runs them in batch mode

import os,sys

def suffixesFinder(count):

    '''
    this function finds the suffixes of the reads
    '''

    inputFileName='/Volumes/omics4tb/alomana/projects/ornl/data/fastaFiles/concatenated_%s.fasta'%count
    with open(inputFileName,'r') as f:
        line=f.readline()
        vector=line.split()
        suffix1=vector[-1]

        suffix2=list(suffix1)
        suffix2[0]='2'
        suffix2=''.join(suffix2)

        suffixes=[suffix1,suffix2]

    return suffixes

### MAIN

# 0. user defined variables
print 'defining user variables...'
numberOfFiles=64 # this should be 64

# 1. creating the runner file
for i in range(numberOfFiles):
    
    count=str(i+1).zfill(2)
    print 'creating runner file %s...'%count

    suffixes=suffixesFinder(count)

    runnerFile='runners/runner_%s.txt'%count
    with open(runnerFile,'w') as f:
        
        f.write("load taxGIFile='/Users/alomana/Applications/MEGAN/class/resources/files/gi_taxid-March2015X.bin';\n")
        f.write("set taxUseGIMap=true;\n")
        f.write("load keggGIFile='/Users/alomana/Applications/MEGAN/class/resources/files/gi2kegg-Feb2015X.bin';\n")
        f.write("set keggUseGIMap=true;\n")
        f.write("load seedGIFile='/Users/alomana/Applications/MEGAN/class/resources/files/gi2seed.map';\n")
        f.write("set seedUseGIMap=true;\n")
        f.write("load cogGIFile='/Users/alomana/Applications/MEGAN/class/resources/files/ncbi/ncbi.map';\n")
        f.write("set cogUseGIMap=true;\n")

        f.write("import blastFile='/Volumes/omics4tb/alomana/projects/ornl/data/diamondFiles/nr.85107862/concatenated.strict_%s.m8' "%count)
        f.write("fastaFile='/Volumes/omics4tb/alomana/projects/ornl/data/fastaFiles/concatenated_%s.fasta' "%count)
        f.write("meganFile='/Volumes/omics4tb/alomana/projects/ornl/data/rmaFiles/nr/concatenated_%s.rma' "%count)
        f.write("maxMatches=100 minScore=50.0 maxExpected=0.01 topPercent=10.0 minSupport=1 minComplexity=0.0 useMinimalCoverageHeuristic=false useSeed=true useCOG=true useKegg=true paired=true ")
        f.write("suffix1='%s' suffix2='%s' "%(suffixes[0],suffixes[1]))
        f.write("useIdentityFilter=false textStoragePolicy=InRMAZ blastFormat=BlastTAB mapping='Taxonomy:BUILT_IN=true,Taxonomy:GI_MAP=true,SEED:GI_MAP=true,KEGG:GI_MAP=true,COG:GI_MAP=true,COG:REFSEQ_MAP=true';\n")
        f.write("quit;\n")

    # 2. executing
    print 'executing runner file %s...'%count
    cmd='/Users/alomana/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E -c /Users/alomana/gDrive2/projects/ornl/src/%s'%runnerFile

    print
    print cmd
    print

    os.system(cmd)
    print
    print '%s job completed.'%count

    
