### this file reads the structure of directories of the trimmed sequence files and creates a single FASTA file for each sample

import os,sys

# 0. user defined variables
path2OriginalFiles='/Volumes/omics4tb/alomana/scratch/tempo/ORNL_Seqs/'
path2fastaFiles='/Volumes/omics4tb/alomana/projects/ornl/data/fastaFiles/'

# 1. defining the data paths
print 'defining paths to original data...'
elements=os.listdir(path2OriginalFiles)
paths=[element for element in elements if 'sort' in element]

# 2. converting trimmed files into
for path in paths:
    label=path.split('_')[0]
    print 'working with sample',label,'...'

    # 2.1. uncompressing files
    print '\t uncompressing files...'
    elements=os.listdir(path2OriginalFiles+path)
    gunzipFiles=[element for element in elements if '_paired_trimmed.fastq.gz' in element]
    for file in gunzipFiles:
        os.system()

    # 2.2. concatenating files
    print '\t concatenating files...'

    # 2.3. converting to fasta file
    print '\t converting FASTQ into FASTA file...'

    sys.exit()
