fastvar
=======

Some quick 'n dirty variant calling tools

Installation
------------

    git clone http://github.com/txje/fastvar
    cd fastvar
    mkdir incl
    cd incl
    
    # install htslib
    git clone http://github.com/samtools/htslib
    cd htslib
    autoheader
    autoconf
    ./configure
    make
    make install
    
    # get klib headers
    git clone http://github.com/attractivechaos/klib
    
    cd ..
    make

Usage
-----

    ./fasttype <BAM> <BED> <REF FASTA>

BAM: reads aligned to &lt;REF FASTA&gt; or '-' to read SAM format from the input stream

BED: list of variant positions to capture in BED format (like dbSNP's [discontinued] batch query results)

    <chrom>  <st> <en> <snp>  <score> <strand>
    chr7  6026829 6026830 rs2639  0 +

We actually only care about the first two fields - chrom and start position.

IMPORTANT: in this format, for whatever reason, the canonical reference position of the SNP is the *end*, and the start is (pos - 1). In practice, we just use the "start" as an index into the 0-indexed reference, so if you're off by one, most of the time all of your loci will appear homozygous and there will be little or no differences between samples. If you see this, *double-check your BED file*. See data/exomeChip\_fingerprint\_snps.bed as an example.

REF FASTA: fasta file used for alignments


Or, directly from one or more unaligned read files, using a k-mer approach:

    ./ktype <BED> <REF FASTA> <READ FASTA/Q> [<READ FASTA/Q> ...]

READ FASTA/Q: may be gzipped


Output
------

A tab-separated file in the following format:

    chrom	start	end	ref_allele	A	C	G	T	N	Unknown

Where the first three fields match the input BED file, followed by the allele from the reference fasta, and counts for each read allele observed at that position, which - in rare cases - may include an "Unknown" that accounts for all manner of sins, including deletions.


Notes
-----

For the cleanest and fastest rapid genotyping, here's what I recommend:

    git clone https://github.com/lh3/minimap2
    cd minimap2
    make
    cd ..

    
    ref=h38.fa                                             # or something
    r1=SAMPLE_A_R1.fastq.gz                                # or something
    r2=SAMPLE_A_R2.fastq.gz                                # or something
    snp_bed=fastvar/data/exomeChip_fingerprint_snps.bed    # or something
    
    ./minimap2/minimap2 -t $(nproc) -ax sr $ref $r1 $r2 | ./fastvar/fasttype - $snp_bed $ref > SAMPLE_A.fast_variants

This will be limited in time and memory by the alignment in most cases (except for *very* large BED SNP files)

I've included a small set of common exome genotyping SNPs (data/exomeChip\_fingerprint\_snps.bed) (see https://genome.sph.umich.edu/wiki/Exome_Chip_Design) that are typically sufficient to match and distinguish samples by individual for fast sanity checks (read barcode mixups).

There is also another, much bigger set of exome SNPs with AF 0.45 - 0.55 from 1000 Genomes data: data/1kgenomes\_0.5af\_exome\_snps.bed

Feel free to submit issues, pull requests, and hate mail.


