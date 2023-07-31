
kunmixer example
================

This example uses the excellent, public, and largely unrestricted dataset of matched tumor-normal pairs part of the [Texas Cancer Research Biobank (TCRB)](http://stegg.hgsc.bcm.edu/open.html). Registration and agreement not to attempt to re-identify the participants is required to access these data, but no additional DUA.

Steps include cloning the kunmixer repository, downloading the hg38 reference genome and the TCRB data (FASTQ files), running kunmixer (ktype) to generate sample fingerprints/profiles, and comparing these fingerprints to identify sample similarities and matches. You will need to register for access to TCRB via the link above and replace "USERNAME" with your username - note that the sftp command will ask for your password on the command line.


    git clone https://github.com/jwanglab/kunmixer
    cd kumnmixer/example
    mkdir -p data
    
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    
    sftp -c aes128-cbc -oPort=9910 USERNAME@hgsc-sftp1.hgsc.bcm.tmc.edu:case_*/*.fastq.bz2 ./data/
    
    for f in TCRBOA1-N-WEX TCRBOA1-T-WEX TCRBOA3-N-WEX TCRBOA3-T-WEX TCRBOA4-N-WEX TCRBOA4-T-WEX TCRBOA5-N-WEX TCRBOA5-T-WEX TCRBOA6-N-WEX TCRBOA6-N-WGS.lane1 TCRBOA6-N-WGS.lane2 TCRBOA6-T-WEX TCRBOA6-T-WGS.lane1 TCRBOA7-N-WEX TCRBOA7-N-WGS.lane1 TCRBOA7-N-WGS.lane2 TCRBOA7-T-RNA TCRBOA7-T-WEX TCRBOA7-T-WGS.lane1 TCRBOA6-T-WGS.lane2 TCRBOA7-T-WGS.lane2 TCRBOA6-T-WGS.lane3 TCRBOA7-T-WGS.lane3 TCRBOA7-T-WGS.lane4 TCRBOA6-T-WGS.lane4 TCRBOA2-N-WEX TCRBOA2-T-WEX; do
      ../ktype 21 ../data/1kgenomes_0.5af_exome_snps.bed hg38.fa.gz data/${f}.read1.fastq.bz2 data/${f}.read2.fastq.bz2 > ${f}.ktype
    done
    
    python ../compareByRatioV4.py TCRBOA1-N-WEX.ktype TCRBOA1-T-WEX.ktype TCRBOA2-N-WEX.ktype TCRBOA2-T-WEX.ktype TCRBOA3-N-WEX.ktype TCRBOA3-T-WEX.ktype TCRBOA4-N-WEX.ktype TCRBOA4-T-WEX.ktype TCRBOA5-N-WEX.ktype TCRBOA5-T-WEX.ktype TCRBOA6-N-WEX.ktype TCRBOA6-T-WEX.ktype TCRBOA6-N-WGS.lane1.ktype TCRBOA6-N-WGS.lane2.ktype TCRBOA6-T-WGS.lane1.ktype TCRBOA6-T-WGS.lane2.ktype TCRBOA6-T-WGS.lane3.ktype TCRBOA6-T-WGS.lane4.ktype TCRBOA7-N-WEX.ktype TCRBOA7-T-WEX.ktype TCRBOA7-N-WGS.lane1.ktype TCRBOA7-N-WGS.lane2.ktype TCRBOA7-T-WGS.lane1.ktype TCRBOA7-T-WGS.lane2.ktype TCRBOA7-T-WGS.lane3.ktype TCRBOA7-T-WGS.lane4.ktype TCRBOA7-T-RNA.ktype


Several output files will be created:

- SimpleSampleFingerPrint.csv
    - Comma-separated matrix including all sample fingerprints
- kunmixer_SNP_check_result.csv
    - Comma-separated list of sample pairings and scores
- compare_matrix.png
    - Heatmap of pairwise similarities
- kunmixer_matchRate_55_distribution.png
    - Histogram distribution of pairwise similarities

Kunmixer results for the TCRB dataset are shown below, clearly identifying matched tumor and normal samples and illustrating its utility across and within genome, exome, and transcriptome sequencing datasets.

![Kunmixer pairwise distance heatmap](/example/compare_matrix.png)

![Kunmixer pairwise distance histogram](/example/kunmixer_matchRate_55_distribution.png)
