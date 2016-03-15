### this script removes all hits from DIAMOND whose identity is lower than 0.70

import sys

inputFileDir='/Volumes/omics4tb/alomana/projects/ornl/data/diamondFilesRefSeq/'
numberOfFiles=64
threshold=70.


for i in range(1,numberOfFiles+1):

    print 'working with file',i,'...'

    label=str(i).zfill(2)

    inputFile=inputFileDir+'concatenated_%s.m8'%label
    outputFile=inputFile.replace('_','.strict_')

    nT=0
    nS=0

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
    print 'total number of hits',nT
    print 'number of accepted hits',nS
    print 'frequency of acceptance',freq
    print
    
