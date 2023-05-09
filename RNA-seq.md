```bash
sudo apt install trim-galore
```

```bash
brew install brewsci/bio/subread
brew install brewsci/bio/hisat2
brew install brewsci/bio/gffread
brew install samtools
```

```bash
mkdir ~/data/rnaseq/GENOMES
```

```bash
cd ~/data/rnaseq
cat *.txt | wc -l
#441

cat << EOF > source.csv
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
            >> source.csv
    '

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

# check them on the screen
#mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt
```
