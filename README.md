# erebus
Tools for quantitative metagenomic analysis of underground samples.  
    
fileHandler.py; script to uncompress, merge and convert FASTQ into FATA files.  
diamondCaller.py; script to call DIAMOND and convert reads into BLAST hits.  
blastCleaner.py; script to remove hits below a certain identity threshold.  
rmaGenerator.py; script to call MEGAN5 in batch mode and convert BLAST files into rma files.  
rmaExtractor.py; script to call MEGAN5 in batch mode to extract hits per species from rma files.  
abundancesGrapher.py; script to generate the heatmap figures of taxa abundances.  
