### this script generate heatmap figures for domains and phyla for each well and size fraction

import sys,os

def broadAbundancesReader():

    '''
    this function reads the abundances of samples down to domain and phyla
    '''

    flotilla={}

    samples=os.listdir(rootDataDir)
    for sample in samples:

        path=rootDataDir+sample
        label=sample.split('_')[1].split('.')[0]
        noahsArk,taxonomy=sampleReader(path)

        flotilla[label]=[noahsArk,taxonomy]

    return flotilla

def metaDataReader():

    '''
    this function reads builds the metadata
    it associates to each sample a well, a filter size and a sorted population
    '''

    metaData={}
    allWells=[]
    allSizes=[]

    with open(metaDataFile,'r') as f:
        allData=f.readline()
        lines=allData.split('\r')
        for i in range(1,len(lines)):
            vector=lines[i].split('\t')

            # sample name
            sampleName=vector[0].split('_')[0]

            # well name
            wellName=vector[3]

            # filter size
            if vector[4] == 'NA':
                filterSize=0
            else:
                filterSize=int(vector[4])

            # sorted population id
            tag=vector[6]
            if 'P' in tag:
                counter=int(tag.split('P')[1])
                tag='P%02d'%counter
            sortedPopulation=tag

            # populating variable
            metaData[sampleName]={}
            metaData[sampleName]['well']=wellName
            metaData[sampleName]['size']=filterSize
            metaData[sampleName]['population']=sortedPopulation

            if wellName not in allWells and wellName != 'FCMbackground':
                allWells.append(wellName)
            
            if filterSize not in allSizes and filterSize != 0:
                allSizes.append(filterSize)

    return metaData,allWells,allSizes

def sampleReader(inputFileName):

    noahsArk={}
    taxonomy=[]

    reading=False

    # f.1. retrieving the data
    with open(inputFileName,'r') as f:
        for line in f:
            vector=line.split('  ')
            cleanVector=[]
            for element in vector:
                element=element.replace('\n','')
                elements=element.split(': ')
                for item in elements:
                    cleanVector.append(item)

            # choosing which lines to read
            if cleanVector[0] == 'root':
                reading=True
                totalReads=float(cleanVector[1])
            if cleanVector[0][:3] == '+++':
                reading=False
                
            # filling up Noah's ark. Building a dictionary with unique names and number of reads and a list for the structure of the taxonomy.
            if reading == True:
                # dealing with the taxonomy
                taxonomyDepth=len(cleanVector)-1
                if taxonomyDepth < taxonomyResolution:
                    for i in range(taxonomyDepth):
                        putative=cleanVector[i]
                        # creating a list if it is a new taxon level
                        if putative != '':
                            if len(taxonomy) < i+1:
                                taxonomy.append([])
                            # adding the taxon if it does not exists
                            if putative not in taxonomy[i]:
                                taxonomy[i].append(putative)

                    # adding the number of reads to the taxon
                    numberOfReads=int(cleanVector[-1])
                    verboseTaxon=[taxonomy[taxon][-1] for taxon in range(taxonomyDepth)]
                    tag='/'.join(verboseTaxon)
                    print tag
                    noahsArk[tag]=numberOfReads

    sys.exit()

    return noahsArk,taxonomy

### MAIN

# 0. preliminaries
metaDataFile='/Volumes/omics4tb/alomana/projects/ornl/data/metadata/ORNL_Sept2015_sequence_metadata.txt'
rootDataDir='/Volumes/omics4tb/alomana/projects/ornl/data/hitsFilesRefSeq/'
taxonomyResolution=5 # a value of 5 reads taxonomy down to phylum
abundanceThreshold=5e-3

# 0.1. metadata reader
print 'reading metadata...'
metaData,allWells,allSizes=metaDataReader()

# 0.2. abundances reader
print 'reading data...'
flotilla=broadAbundancesReader()

# 1. defining subsets of samples
print 'treating the data...'

for well in allWells:
    for size in allSizes:
        
        print 'working with well %s, size %s...'%(well,str(size))

        # ordering the columns
        selectedSamples={}
        for sampleID in metaData.keys():
            if metaData[sampleID]['well'] == well and metaData[sampleID]['size'] == size:
                selectedSamples[sampleID]=metaData[sampleID]['population']

        sortedColumns=sorted(selectedSamples,key=selectedSamples.__getitem__)

        # merging the total number of reads for each sample into an ensemble taxonomy
        ensembleSet={}
        for sampleID in sortedColumns:
            
            noahsArk=flotilla[sampleID][0]
            taxonomy=flotilla[sampleID][1]

            taxa=noahsArk.keys()
            for taxon in taxa:
                if taxon in ensembleSet.keys():
                    ensembleSet[taxon]=ensembleSet[taxon]+noahsArk[taxon]
                else:
                    ensembleSet[taxon]=noahsArk[taxon]

        # computing the relative frequencies
        totalN=float(ensembleSet['root'])
        for taxon in ensembleSet.keys():
            print taxon
            ensembleSet[taxon]=float(ensembleSet[taxon])/totalN

            if ensembleSet[taxon] < abundanceThreshold:
                del ensembleSet[taxon]

        for element in ensembleSet.keys():
            print element,ensembleSet[element]
            

            
        sys.exit()
