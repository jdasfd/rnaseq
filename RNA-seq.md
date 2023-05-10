# RNA-seq pipeline

Using RNA-seq of _medicago_truncatula_ to identify RLK genes expressed in the rhizobial symbiosis.

## Prepare

```bash
cargo install --git https://github.com/wang-q/anchr --branch main
anchr help
```

```bash
sudo apt install trim-galore
brew install brewsci/bio/subread
brew install brewsci/bio/hisat2
brew install brewsci/bio/gffread
brew install samtools
```

## RNA-seq pipeline

An automatically running RNA-seq scripts named `rnaseq_auto.pl` is written for the goal.

- Acquiring the Mtru genome and annotation

```bash
mkdir ~/data/rnaseq/GENOMES
cd ~/data/rnaseq/GENOMES

# mannually download them from IE
# here old genomes are used directly
cp -r ~/data/RLK_family/GENOMES/medicago_truncatula/ ./
```

- Selecting RNA-seq files

All RNA-seq `PRJNA*.txt` were downloaded from the NCBI. Go and check dir info.

```bash
mkdir -p ~/data/rnaseq/info
# all PRJNA*.txt were saved here

mkdir -p ~/data/rnaseq/sra
cd ~/data/rnaseq/info
cat *.txt | wc -l
#441

cat << EOF > ../sra/source.csv
#Experiment,Sample_Name,Bases
EOF

ls *.txt |
    parallel -j 1 '
        cat {} |
            mlr --icsv --otsv cat |
            tsv-select -H -f Experiment,"Sample\ Name",Bases |
            tr " " "_" |
            sed '\''1 s/^/#/'\'' |
            keep-header -- sort -k2,2 -k3,3nr |
            tsv-uniq -H -f Sample_Name --max 1 |
            sed 1d |
            mlr --itsv --ocsv cat \
            >> ../sra/source.csv
    '

cd ~/data/rnaseq/sra

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

# check them on the screen
#mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt
```

- Getting the count result from raw reads

```bash
cd ~/data/rnaseq/sra

ls *.fastq.gz | wc -l
#591

ls *_1.fastq.gz *_2.fastq.gz |
    perl -pe 's/_\d\.fastq\.gz$//' |
    tsv-uniq \
    > ../info/PE_SRR.lst
cat ../info/PE_SRR.lst | wc -l
#155

ls *.fastq.gz |
    perl -nle 'next if /.+_\d\.fastq\.gz$/; print;' |
    perl -pe 's/\.fastq\.gz$//' \
    > ../info/SE_SRR.lst
cat ../info/SE_SRR.lst | wc -l
#281

# echo '155*2+281' | bc
#591

cd ~/data/rnaseq

# single end
cat info/SE_SRR.lst |
    parallel -j 1 -k '
        perl scripts/rnaseq_auto.pl --thread 20 -i sra/{}.fastq.gz \
            -g GENOMES/medicago_truncatula/genome.fa \
            -a GENOMES/medicago_truncatula/genome.gtf
    '

# paired end
cat info/PE_SRR.lst |
    parallel -j 1 -k '
        perl scripts/rnaseq_auto.pl --thread 20 -t PE \
            -i sra/{}_1.fastq.gz -i sra/{}_2.fastq.gz \
            -g GENOMES/medicago_truncatula/genome.fa \
            -a GENOMES/medicago_truncatula/genome.gtf
    '
```
