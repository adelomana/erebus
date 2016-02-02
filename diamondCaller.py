### this script reads the input fasta files and calls DIAMOND. The output files are in .m8 format. Script to be run in aegir under SGE environment.

import os,sys

def runnerCreator(tag):

    '''
    this function creates the runner files of aegir SGE system
    '''

    minimalTag='c%s'%(tag.split('_')[1])
    
    inputFile='sgeRunners/%s.sh'%tag
    with open(inputFile,'w') as g:
        g.write('#!/bin/bash\n\n')
        g.write('#$ -N %s\n'%minimalTag)
        g.write('#$ -o %s/messagesDIAMOND_%s.o.txt\n'%(scratchDir,minimalTag))
        g.write('#$ -e %s/messagesDIAMOND_%s.e.txt\n'%(scratchDir,minimalTag))
        g.write('#$ -P Bal_alomana\n')
        g.write('#$ -pe serial 20\n')
        g.write('#$ -q baliga\n')
        g.write('#$ -S /bin/bash\n\n')
        g.write('cd /users/alomana\n')
        g.write('source .bash_profile\n\n')

        cmd=diamondPath+' blastx -d '+nrPath+' -q '+fastaFilesDir+tag+'.fasta -a '+diamondOutputDir+tag+' -t '+scratchDir+' --threads '+str(threads)+' --sensitive'
        g.write('%s\n\n'%cmd)

        cmd=diamondPath+' view -a '+diamondOutputDir+tag+'.daa -o '+diamondOutputDir+tag+'.m8'
        g.write('%s\n\n'%cmd)


    g.close()

    return None

# 0. user defined variables
threads=20
diamondPath='/proj/omics4tb/alomana/software/diamond_linux64binary_v0.7.9/diamond'
nrPath='/proj/omics4tb/alomana/projects/rossSea/data/metagenomics/db/nr'
fastaFilesDir='/proj/omics4tb/alomana/projects/ornl/data/fastaFiles/'
diamondOutputDir='/proj/omics4tb/alomana/projects/ornl/data/diamondFiles/'
scratchDir='/proj/omics4tb/alomana/scratch/diamond/'

# 1. define the inputs
inputFiles=os.listdir(fastaFilesDir)

# 2. create launching the SGE calling files
for inputFile in inputFiles:
    tag=inputFile.split('.')[0]
    runnerCreator(tag)

    # 2.1. launching
    cmd='qsub sgeRunners/%s.sh'%tag
    os.system(cmd)
    sys.exit()

print '... all done.'
