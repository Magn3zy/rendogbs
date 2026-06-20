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