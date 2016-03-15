### this script removes all hits from DIAMOND whose identity is lower than 0.70.
### the script expects a parallel environment, 5 threads are launched

import sys,multiprocessing
from multiprocessing import pool

def selector(label):

    '''
    this function writes a new file with the selected hits
    '''

    inputFile=inputFileDir+'concatenated_%s.m8'%label
    outputFile=inputFile.replace('_','.strict_')
    
    freq=0.
    nS=0
    nT=0

    with open(inputFile,'r') as f:
        with open(outputFile,'w') as g:
            for line in f:
                vector=line.split('\t')
                identity=float(vector[2])
                nT=nT+1
                if identity > threshold:
                    g.write(line)
                    nS=nS+1
                    
    freq=float(nS)/float(nT)

    return freq,nS,nT

### MAIN

# 0. user defined variables
inputFileDir='/Volumes/omics4tb/alomana/projects/ornl/data/diamondFilesRefSeq/'
numberOfFiles=64
threshold=70.

# 1. main
print 'computing...'

all_freq=[]
all_nS=[]
all_nT=[]

labels=[str(i).zfill(2) for i in range(1,numberOfFiles+1)]

# entering a parallel instance
hydra=multiprocessing.pool.Pool(5)
output=hydra.map(selector,labels)
    
# 2. messages printing
print
for i in range(len(output)):

    print 'about file %s ...'%str(i+1)
    print 'total number of hits %s'%(output[i][2])
    print 'number of accepted hits %s'%(output[i][1])
    print 'frequency of acceptance %s'%(output[i][0])
    print

allFreq=output[:][0]
print 'most successful frequencty',max(allFreq)
print 'least successful frequency',min(allFreq)
    
