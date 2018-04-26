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
