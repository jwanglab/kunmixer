kunmixer
=======

Kunmixer uses a set of maximally informative *k*-mers as in silico probes to estimate genomic identity across samples.
The intent of this project is not to create *accurate* genotypes, per se, as there are a great many good pieces of software for doing this, but to create the fastest possible genotypes that are passable to assess sample integrity, match samples from the same person, and approximate relatedness among individuals.

Installation
------------

Dependencies:

[HTSLib v1.10+] (https://github.com/samtools/htslib/releases/tag/1.13) must be [installed] (https://github.com/samtools/htslib/blob/develop/INSTALL) globally on your machine, or if you don't have root access to your machine, tell GCC where a local htslib is.

For example, if you have built htslib in `/home/user/htslib`:

    export LIBRARY_PATH=/home/usr/htslib/htslib:$LIBRARY_PATH
    export CPATH=/home/usr/htslib/htslib:$CPATH

Install:

    git clone http://github.com/jwanglab/kunmixer
    cd kunmixer
    
    # get klib headers
    mkdir incl
    cd incl
    git clone http://github.com/attractivechaos/klib
    cd ..
    
    make


Usage
-----

    ./ktype <BED> <REF FASTA> <READ FASTA/Q> [<READ FASTA/Q> ...]

BED: list of variant positions to capture in BED format (like dbSNP's [discontinued] batch query results)

    <chrom>  <st> <en> <snp>  <score> <strand>
    chr7  6026829 6026830 rs2639  0 +

We actually only care about the first two fields - chrom and start position.

IMPORTANT: in this format, for whatever reason, the canonical reference position of the SNP is the *end*, and the start is (pos - 1). In practice, we just use the "start" as an index into the 0-indexed reference, so if you're off by one, most of the time all of your loci will appear homozygous and there will be little or no differences between samples. If you see this, *double-check your BED file*. See data/exomeChip\_fingerprint\_snps.bed as an example.

REF FASTA: fasta file used for alignments

READ FASTA/Q: may be gzipped or bzip2'd

Or, a slower but potentially more precise using aligned reads:

    ./fasttype <BAM> <BED> <REF FASTA>

BAM: reads aligned to &lt;REF FASTA&gt; or '-' to read SAM format from the input stream

Here's an example run script:

    ref=h38.fa
    r1=SAMPLE_A_R1.fastq.gz
    r2=SAMPLE_A_R2.fastq.gz
    snp_bed=kunmixer/data/exomeChip_fingerprint_snps.bed    # or something

    ./kunmixer/ktype $snp_bed $ref $r1 $r2 > SAMPLE_A.k_variants

    # -- OR -- (using minimap2)
    
    git clone https://github.com/lh3/minimap2
    cd minimap2
    make
    cd ..

    ./minimap2/minimap2 -t $(nproc) -ax sr $ref $r1 $r2 | ./kunmixer/fasttype - $snp_bed $ref > SAMPLE_A.fast_variants

Included is a small set of common exome genotyping SNPs (data/exomeChip\_fingerprint\_snps.bed) (see https://genome.sph.umich.edu/wiki/Exome_Chip_Design) that are typically sufficient to match and distinguish samples by individual for fast sanity checks (read barcode mixups).

There is also another, much bigger set of exome SNPs with AF 0.45 - 0.55 from 1000 Genomes data: data/1kgenomes\_0.5af\_exome\_snps.bed

Both are based on hg38 and will only work given the appropriate reference.


Output
------

A tab-separated file in the following format:

    chrom	start	end	ref_allele	A	C	G	T	N	Unknown

Where the first three fields match the input BED file, followed by the allele from the reference fasta, and counts for each read allele observed at that position, which - in rare cases - may include an "Unknown" that accounts for all manner of sins, including deletions.


Performance
-----------

For both fasttype and ktype, memory usage is typically dominated the size of the reference (~3GB for humans).
If you have *very* many SNPs, the variant bookkeeping will begin to contribute.

Runtime is typically dominated by the BAM or read FASTQ size. In our hands, for high-depth (>100x) exome sequencing, fasttype takes ~2 min per sample and ktype takes ~7 min per sample from gzipped reads (on one thread).

Alignment-based fasttype results correspond very closely to the ktype results in most cases. For the limited trials, with relatively few SNPs, described below, I've tried k = 32 and k = 64, and there's little difference. There are still some undiagnosed differences where, in some cases, ktype finds anywhere from ~1/2 - 2x as many of one or both nucleotides as fasttype. For simple sample matching, this isn't a huge problem because at least it's consistent across samples, but is certainly not the expected behavior.


Example
-------

See a [full example](https://github.com/jwanglab/kunmixer/tree/master/example) including run script and expected output.

This example uses the excellent, public, and largely unrestricted dataset of matched tumor-normal pairs part of the [Texas Cancer Research Biobank (TCRB)](http://stegg.hgsc.bcm.edu/open.html). Registration and agreement not to attempt to re-identify the participants is required to access these data, but no additional DUA.

Kunmixer results for this dataset are shown below, clearly identifying matched tumor and normal samples and illustrating its utility across and within genome, exome, and transcriptome sequencing datasets.

![Kunmixer pairwise distance heatmap](/example/compare_matrix.png)

![Kunmixer pairwise distance histogram](/example/kunmixer_matchRate_55_distribution.png)



Contributions
-------------

Feel free to submit issues and pull requests.
