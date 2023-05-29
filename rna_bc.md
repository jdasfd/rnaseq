# Breast cancer RNA-seq data

Extract KEGG pathway genelist expression matrix from a breast cancer article. Article link: [Cell, 2016, 164: 293-309](https://www.cell.com/cell/fulltext/S0092-8674(15)01624-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415016244%3Fshowall%3Dtrue).

## Preparation

- Software

```bash
cd ~/Scripts
git clone https://github.com/wang-q/fig_table.git
# check how to use the script
perl ~/Scripts/fig_table/xlsx2csv.pl -h
```

- [Convert CR 2 LF](https://toolslick.com/conversion/text/new-line). The web-tool for converting CR to LF EOF.

- R packages

```bash
Rscript -e '
    BiocManager::install("KEGGREST")
    BiocManager::install("EnrichmentBrowser")
    '
```

## Data derived from the origin analyzation

### Extract all necessary data

- Data orginated from [neellab/bfg](https://github.com/neellab/bfg/tree/gh-pages).

```bash
mkdir -p ~/data/rna_bc
cd ~/data/rna_bc

# bfg
git clone https://github.com/neellab/bfg.git
```

- Info derived from article supplementary.

The data was collected from the supplementary material directly. All `.xlsx` were collected into `article_SM`.

```bash
mkdir -p ~/data/rna_bc/info
cd ~/data/rna_bc/info

cp ../bfg/data/annotations/cell_line_subtypes.txt.zip .
cp ../bfg/data/rnaseq/breast_rnaseq_fpkm_nonnormalized.txt.zip .
cp ../bfg/data/rnaseq/breast_rnaseq_qn.txt.zip .

rm -rf __MACOSX/
rm *.zip

# converting cr2lf
for file in cell_line_subtypes breast_rnaseq_qn
do
cat ${file}.txt | tr "\r" "\n" > tmp && mv tmp ${file}.txt
echo >> ${file}.txt
done

# extract cell and its subtype info
cat cell_line_subtypes.txt |
    tsv-select -H -f cell_line,subtype_neve |
    perl -nlae '
        print if /^cell_line/;
        print "$F[0]\tbasal_A" if $F[1] =~ /basala/;
        print "$F[0]\tbasal_B" if $F[1] =~ /basalb/;
        print "$F[0]\tHER2" if $F[1] =~ /her2/;
        print "$F[0]\tLuminal" if $F[1] =~ /luminal/;
    '\
    > cell_anno.tsv

cat > tmp.tsv << EOF 
$(echo -e "H\tbasal_A")
$(echo -e "I\tbasal_B")
$(echo -e "J\tHER2")
$(echo -e "K\tLuminal")
EOF

# extract each gene pathway
rm pathway.tsv gene.tsv
cat tmp.tsv |
    parallel -j 1 --colsep "\t" '
        perl ~/Scripts/fig_table/xlsx2csv.pl \
            -f ../article_SM/1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3{1} |
            sed "1,2d" |
            mlr --icsv --otsv cat |
            tsv-select -H -f "Pathway\ id" |
            sed 1d |
            awk -v GROUP={2} '\''{print ($0"\t"GROUP)}'\'' \
            >> pathway.tsv
        perl ~/Scripts/fig_table/xlsx2csv.pl \
            -f ../article_SM/1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3{1} |
            sed "1,2d" |
            mlr --icsv --otsv cat |
            tsv-select -H -f "Genes" |
            sed 1d |
            tr "," "\n" |
            sort |
            uniq |
            awk -v GROUP={2} '\''{print ($0"\t"GROUP)}'\'' \
            >> gene.tsv
    '
rm tmp.tsv

cat pathway.tsv | tsv-summarize -g 2 --count
#basal_A 580
#basal_B 90
#HER2    264
#Luminal 69

cat gene.tsv | tsv-summarize -g 2 --count
#basal_A 238
#basal_B 104
#HER2    139
#Luminal 83
```

### Extract KEGG pathway via R

```bash
cd ~/data/rna_bc/info

# using KEGGREST
Rscript ../../rnaseq/scripts/kegg_extract.r hsa05224
```

### Extract expression matrices

```bash
cd ~/data/rna_bc/info

for group in fpkm_nonnormalized qn
do
    cat breast_rnaseq_${group}.txt |
        tsv-join -H -f <(
            cat hsa05224.lst |
            sed '1isymbol'
            ) -k symbol \
        > hsa05224_${group}.tsv
done
```
