# get the hg38 reference genome for use with 1kgenomes BED file
if [ ! -s hg38.fa.gz ]; then
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
fi

mkdir -p data

for c in {001..007}; do
  # have to enter password each time...
  sftp -c aes128-cbc -oPort=9910 USERNAME@hgsc-sftp1.hgsc.bcm.tmc.edu:case_${c}/*.fastq.bz2 ./data/
done

for f in TCRBOA1-N-WEX TCRBOA1-T-WEX TCRBOA3-N-WEX TCRBOA3-T-WEX TCRBOA4-N-WEX TCRBOA4-T-WEX TCRBOA5-N-WEX TCRBOA5-T-WEX TCRBOA6-N-WEX TCRBOA6-N-WGS.lane1 TCRBOA6-N-WGS.lane2 TCRBOA6-T-WEX TCRBOA6-T-WGS.lane1 TCRBOA7-N-WEX TCRBOA7-N-WGS.lane1 TCRBOA7-N-WGS.lane2 TCRBOA7-T-RNA TCRBOA7-T-WEX TCRBOA7-T-WGS.lane1 TCRBOA6-T-WGS.lane2 TCRBOA7-T-WGS.lane2 TCRBOA6-T-WGS.lane3 TCRBOA7-T-WGS.lane3 TCRBOA7-T-WGS.lane4 TCRBOA6-T-WGS.lane4 TCRBOA2-N-WEX TCRBOA2-T-WEX; do
  if [ ! -s ${f}.ktype ]; then
    ../ktype 21 ../data/1kgenomes_0.5af_exome_snps.bed hg38.fa.gz data/${f}.read1.fastq.bz2 data/${f}.read2.fastq.bz2 > ${f}.ktype
  fi
done

python ../compareByRatioV4.py \
  TCRBOA1-N-WEX.ktype \
  TCRBOA1-T-WEX.ktype \
  TCRBOA2-N-WEX.ktype \
  TCRBOA2-T-WEX.ktype \
  TCRBOA3-N-WEX.ktype \
  TCRBOA3-T-WEX.ktype \
  TCRBOA4-N-WEX.ktype \
  TCRBOA4-T-WEX.ktype \
  TCRBOA5-N-WEX.ktype \
  TCRBOA5-T-WEX.ktype \
  TCRBOA6-N-WEX.ktype \
  TCRBOA6-T-WEX.ktype \
  TCRBOA6-N-WGS.lane1.ktype \
  TCRBOA6-N-WGS.lane2.ktype \
  TCRBOA6-T-WGS.lane1.ktype \
  TCRBOA6-T-WGS.lane2.ktype \
  TCRBOA6-T-WGS.lane3.ktype \
  TCRBOA6-T-WGS.lane4.ktype \
  TCRBOA7-N-WEX.ktype \
  TCRBOA7-T-WEX.ktype \
  TCRBOA7-N-WGS.lane1.ktype \
  TCRBOA7-N-WGS.lane2.ktype \
  TCRBOA7-T-WGS.lane1.ktype \
  TCRBOA7-T-WGS.lane2.ktype \
  TCRBOA7-T-WGS.lane3.ktype \
  TCRBOA7-T-WGS.lane4.ktype \
  TCRBOA7-T-RNA.ktype
