#!/bin/bash
#AJOUTER VERSION MMSEQS en fonction de ce qui est utilisé (env conda, img singularity, docker etc)

UNIPROT50="/dsimb/cafe/DATABASES/uniref50/uniref50.fasta"
DB="/home/chili/chalvin/projet_long/uniprot50to30/uniprot50.mmseqsdb"
DB_clu="/home/chili/chalvin/projet_long/uniprot50to30/uniprot30.mmseqsdb"
DB_clu_rep="/home/chili/chalvin/projet_long/uniprot50to30/uniprot30.rep"
DATE=$(date '+%Y-%m-%d %H:%M:%S')
LOG="/home/chili/chalvin/projet_long/uniprot50to30/uniprot50to30.log"
TMP_DIR="tmp"


echo $DATE >> $LOG
#mmseqs createdb $UNIPROT50 $DB
#cov-mode 0 : > 80% sequence of both query and target aligned
#mmseqs linclust --min-seq-id 0.3 --cov-mode 0 $DB $DB_clu $TMP_DIR 2>&1>>$LOG
mmseqs createsubdb $DB_clu $DB $DB_clu_rep >> $LOG


