### this script selects the samples to be compared and plots heatmaps of different samples with the relevant abundances

import sys,os,numpy,matplotlib
from matplotlib import pyplot

def abundancesReader():

    '''
    this function reads the relative abundances for each sample
    '''

    flotilla={}

    samples=os.listdir(rootDataDir)
    for sample in samples:

        path=rootDataDir+sample
        label=sample.split('_')[1].split('.')[0]
        noahsArk,taxonomy=sampleReader(path)

        flotilla[label]=[noahsArk,taxonomy]

    return flotilla

def matrixDiscreteMaker(t):

    '''
    this function converts the frequencies matrix t into an integer valued-matrix
    rules used:

    less than 1% or nan --> nan
    between 1% and 5% --> 1
    between 5% and 10% --> 2
    between 10% and 20% --> 3
    more than 20% --> 4
    '''

    # less than 1% or nan --> 1
    numpy.place(t,numpy.isnan(t),0)
    numpy.place(t,t<0.01,float('nan'))

    # between 1% and 5% --> 2
    t[numpy.where(numpy.logical_and(t>=0.01, t<0.05))]=1

    # between 5% and 10% --> 3
    t[numpy.where(numpy.logical_and(t>=0.05, t<0.1))]=2

    # between 10% and 20% --> 4
    t[numpy.where(numpy.logical_and(t>=0.1, t<0.2))]=3

    # more than 20% --> 5
    t[numpy.where(numpy.logical_and(t>=0.2, t<1-1e-10))]=4

    return t

def discrete_cmap(N, base_cmap=None):
    
    '''
    Create an N-bin discrete colormap from the specified input map
    Modified from https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    '''

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = matplotlib.pyplot.cm.get_cmap(base_cmap)
    color_list = base(numpy.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    
    return base.from_list(cmap_name, color_list, N)

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
            sortedPopulation=vector[6]

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
                relativeAbundance=float(numberOfReads)/totalReads
                verboseTaxon=[taxonomy[taxon][-1] for taxon in range(taxonomyDepth)]
                tag='/'.join(verboseTaxon)
                noahsArk[tag]=relativeAbundance

    return noahsArk,taxonomy

### MAIN

# 0. preliminaries
metaDataFile='/Volumes/omics4tb/alomana/projects/ornl/data/metadata/ORNL_Sept2015_sequence_metadata.txt'
rootDataDir='/Volumes/omics4tb/alomana/projects/ornl/data/hitsFiles/'
abundanceThreshold=5e-2
nMax=50

# 0.1. metadata reader
print 'reading metadata...'
metaData,allWells,allSizes=metaDataReader()

# 0.2. abundances reader
print 'reading data...'
flotilla=abundancesReader()

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

        # selecting the most abundant taxa over each of the samples
        abundantTaxa=[]
        for sampleID in sortedColumns:
            
            noahsArk=flotilla[sampleID][0]
            taxonomy=flotilla[sampleID][1]

            taxa=noahsArk.keys()
            for taxon in taxa:
                if noahsArk[taxon] > abundanceThreshold and taxon not in abundantTaxa:
                    abundantTaxa.append(taxon)
        print len(abundantTaxa),'taxa found above the abundance threshold.'

        # selecting the unique taxa along all samples
        deepTaxa=[]
        for taxon in abundantTaxa:
            # is there a deeper taxon, then it is not deep enough
            deep=0
            for element in abundantTaxa:
                if element.find(taxon) != -1:
                    deep=deep+1
            if deep == 1:
                deepTaxa.append(taxon)
        print len(deepTaxa),'final unique taxa found with abundance higher than threshold.'

        # computing accumulated abundances, used as a proxy for ranking taxa
        accumulatedAbundances={}
        sortedAccumulatedAbundances=[]
        for taxon in deepTaxa:
            f=0.
            for sampleID in sortedColumns:
                noahsArk=flotilla[sampleID][0]
                if taxon in noahsArk:
                    fi=noahsArk[taxon]
                else:
                    fi=0.
                f=f+fi
            accumulatedAbundances[taxon]=f
            sortedAccumulatedAbundances.append(f)
        sortedAccumulatedAbundances.sort(reverse=True)

        # selecting the 50 or less most abundant taxa
        rowNames=[]
        if len(deepTaxa) > nMax:
            threshold=sortedAccumulatedAbundances[nMax-1]
        else:
            threshold=sortedAccumulatedAbundances[-1]
        for taxon in deepTaxa:
            if accumulatedAbundances[taxon] >= threshold:
                rowNames.append(taxon)
        rowNames.sort()
        print len(rowNames),'taxa selected for plotting...'

        # building the variables required for the figure
        columnNames=[selectedSamples[sampleID] for sampleID in sortedColumns]
        finalRowNames=[rowName.split('/')[-1] for rowName in rowNames]

        m=[]
        for sampleID in sortedColumns:
            noahsArk=flotilla[sampleID][0]
            v=[]
            for rowName in rowNames:
                if rowName in noahsArk:
                    value=noahsArk[rowName]
                else:
                    value=float('nan')
                v.append(value)
            m.append(v)
        m=numpy.array(m)
        t=numpy.transpose(m)
        discreteT=matrixDiscreteMaker(t)
        
        # plotting the figure
        colorResolution=4
        discrete_cmap = matplotlib.pyplot.get_cmap('Blues',colorResolution)
        matplotlib.pyplot.imshow(t,interpolation='none',cmap=discrete_cmap,vmin=0.5,vmax=4.5)
        theBar=matplotlib.pyplot.colorbar(label='frequency',orientation='vertical',fraction=0.1,ticks=[0.5,1.5,2.5,3.5,4]) # original fraction was 0.01
        theBar.ax.set_yticklabels(['1%','5%','10%','20%','>20%'])
        
        matplotlib.pyplot.grid(False)
        matplotlib.pyplot.xticks(range(len(columnNames)),columnNames,size=8,rotation=90)
        matplotlib.pyplot.yticks(range(len(finalRowNames)),finalRowNames,size=8)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.tick_params(axis='x',which='both',bottom='off',top='off')
        matplotlib.pyplot.tick_params(axis='y',which='both',right='off',left='off')
        figureName='well_%s_size_%s.png'%(well,size)
        matplotlib.pyplot.savefig('figures/%s'%figureName)
        matplotlib.pyplot.clf()
        
        print
        
        
