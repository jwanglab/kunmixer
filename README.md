fastvar
=======

Some quick 'n dirty variant calling tools

To compile:

    git clone http://github.com/txje/fastvar
    cd fastvar
    mkdir incl
    cd incl
    
    // install htslib
    git clone http://github.com/samtools/htslib
    cd htslib
    autoheader
    autoconf
    ./configure
    make
    make install
    
    // get klib headers
    git clone http://github.com/attractivechaos/klib
    
    cd ..
    make

Usage:

    ./fasttype <BAM> <BED> <REF FASTA>

BAM: reads aligned to <REF FASTA> or '-' to read SAM format from the input stream

BED: list of variant positions to capture in BED format (like dbSNP's [discontinued] batch query results)

    <chrom>  <st> <en> <snp>  <score> <strand>
    chr7  6026829 6026830 rs2639  0 +

We actually only care about the first two fields - chrom and start position.

IMPORTANT: in this format, for whatever reason, the canonical reference position of the SNP is the *end*, and the start is (pos - 1). In practice, we just use the "start" as an index into the 0-indexed reference, so if you're off by one, most of the time all of your loci will appear homozygous and there will be little or no differences between samples. If you see this, *double-check your BED file*.

REF FASTA: fasta file used for alignments

Output: a tab-separated file in the following format:

    chrom	start	end	ref_allele	A	C	G	T	N	Unknown

Where the first three fields match the input BED file, followed by the allele from the reference fasta, and counts for each read allele observed at that position, which - in rare cases - may include an "Unknown" that accounts for all manner of sins, including deletions.

Feel free to submit issues, pull requests, and hate mail.
