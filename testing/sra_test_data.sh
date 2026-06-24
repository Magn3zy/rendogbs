# Barley - citation: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032253#s4
# https://www.ncbi.nlm.nih.gov/sra/SRX112337[accn]
# combination: PstI,MspI
# Genotyping-by-sequencing of Oregon Wolfe Barley DH

fasterq-dump \
SRR391278 \
-O /home/raidfolder/rendogbs/data \
--temp /home/raidfolder/rendogbs/data/tmp \
--threads 8 \
--split-files \
--progress

#Wheat - citation: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032253#s4
# https://www.ncbi.nlm.nih.gov/sra/SRX112339[accn]
# combination: PstI,MspI
# Genotyping-by-sequencing of Synthetic W9784 x Opata M85 DH population

fasterq-dump \
SRR391280 \
-O /home/raidfolder/rendogbs/data \
--temp /home/raidfolder/rendogbs/data/tmp \
--threads 8 \
--split-files \
--progress

# reference: 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz  
# annotation: 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.gff.gz 

sh ./rendogbs.sh \
  --ref /home/raidfolder/rendogbs/pipeline_references/zmays.fna \
  --workdir /home/raidfolder/rendogbs/pipeline_references/zm_rendogbs/ \
  --parallel 11 \
  --size 300-600 \
  --combinations custom \
  --combinations-file /home/raidfolder/rendogbs/pipeline_references/ta_rendogbs/combs_maize.csv \
  --chroms 10 \
  --annotation /home/raidfolder/rendogbs/pipeline_references/zmays.gff

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz 
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz 

cat combs_maize.csv 
enzyme_a,enzyme_b
PstI,MspI
PstI,HaeIII
EcoRI,MspI
HindIII,MspI
HindIII,BfaI
HindIII,MseI
EcoRI,HaeIII
PstI,MseI
HindIII,HaeIII
PstI,BfaI
EcoRI,BfaI

# citation: https://link.springer.com/article/10.1186/s12864-017-3765-8?





# unique fragments
awk -F',' '
NR>1 {
    key=$1","$2
    cov[key] += $5
}
END {
    print "accession,start_pos,coverage"
    for (k in cov)
        print k","cov[k]
}
' SRR5360101_1.csv > SRR5360101_1.unique.csv

./ddrad_test.sh /home/raidfolder/rendogbs/pipeline_references/qr_rendogbs/results/EcoRI_MseI/filtered.csv SRR391278.bam 100
Step 1: loading reference...
Step 2: extracting BAM loci...
Step 3: interval overlap (bedtools)...

Results:
Predicted fragments : 66157
Matched fragments   : 522
Recovery rate       : 0.0079
(base) magn3zy@cipotac2:/home/raidfolder/rendogbs/alignment$ ./ddrad_test.sh /home/raidfolder/rendogbs/pipeline_references/qr_rendogbs/results/EcoRI_MseI/filtered.csv SRR5360103_1.bam 100
Step 1: loading reference...
Step 2: extracting BAM loci...
Step 3: interval overlap (bedtools)...

Results:
Predicted fragments : 66157
Matched fragments   : 36813
Recovery rate       : 0.5564
(base) magn3zy@cipotac2:/home/raidfolder/rendogbs/alignment$ ./ddrad_test.sh /home/raidfolder/rendogbs/pipeline_references/qr_rendogbs/results/EcoRI_MseI/filtered.csv SRR5360104_1.bam 100
Step 1: loading reference...
Step 2: extracting BAM loci...
Step 3: interval overlap (bedtools)...

Results:
Predicted fragments : 66157
Matched fragments   : 36960
Recovery rate       : 0.5587
(base) magn3zy@cipotac2:/home/raidfolder/rendogbs/alignment$ ./ddrad_test.sh /home/raidfolder/rendogbs/pipeline_references/qr_rendogbs/results/EcoRI_MseI/filtered.csv SRR5360101_1.bam 100
Step 1: loading reference...
Step 2: extracting BAM loci...
Step 3: interval overlap (bedtools)...

Results:
Predicted fragments : 66157
Matched fragments   : 32076
Recovery rate       : 0.4848

