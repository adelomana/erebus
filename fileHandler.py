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
    gunzipFiles=[path2OriginalFiles+path+'/'+element for element in elements if '_paired_trimmed.fastq.gz' in element]
    for file in gunzipFiles:
        os.system('gunzip -k -f %s'%file)

    # 2.2. concatenating files
    print '\t concatenating files...'
    concatenatedLabel=path2OriginalFiles+path
    fastaFiles=[element.replace('.gz','') for element in gunzipFiles]
    cmd='cat '
    for element in fastaFiles:
        cmd=cmd+' '+element
    cmd=cmd+' > %s/concatenated.fastq'%(concatenatedLabel)

    os.system(cmd)

    # 2.3. converting to fasta file
    print '\t converting FASTQ into FASTA file...'
    cmd="awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,\">\");print}; if(P==4)P=0; P++}' %s/concatenated.fastq > %s/concatenated_%s.fasta"%(concatenatedLabel,concatenatedLabel,label)
    os.system(cmd)
