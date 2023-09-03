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
mkdir -p ~/data/rnaseq/rna_bc/info
cd ~/data/rnaseq/rna_bc/info

# bfg
#git clone https://github.com/neellab/bfg.git

# GSE96058
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
gzip -d GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/soft/GSE96058_family.soft.gz
gzip -d
```

- Info derived from article supplementary.

The data was collected from the supplementary material directly. All `.xlsx` were collected into `article_SM`.

```bash
cd ~/data/rnaseq/rna_bc/info

# converting cr2lf
#for file in cell_line_subtypes breast_rnaseq_qn
#do
#cat ${file}.txt | tr "\r" "\n" > tmp && mv tmp ${file}.txt
#echo >> ${file}.txt
#done

# extract sample info and scan-b id
cat GSE96058_family.soft |
    perl -nle '
        print if /^!Sample_title.*/;
        print if /^!Sample_characteristics_ch1\s=\sscan-b.*/;
    ' \
    > tmp.txt

cat tmp.txt | perl -e 'while (<>) {chomp; if ($.%2 != 0) { $_ =~ /^!Sample_title\s=\s(.+)$/; print "$1\t";} else { $_ =~ /^.*id:\sQ\d+\.C\d+\.(.+)$/; print "$1\n";}}' > id_trans.tsv

perl ~/Scripts/fig_table/xlsx2csv.pl -f 41598_2019_48570_MOESM2_ESM.xlsx | sed '1,2d' | mlr --icsv --otsv cat | tsv-select -H -f Assay,AIMS_PAM50 | sed 1d | tsv-select -f 2,1 | tsv-join -f id_trans.tsv -k 2 -a 1 | tsv-select -f 3,1 > id_trans_subtype.tsv

rm tmp.txt

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
cd ~/data/rnaseq/rna_bc/info

# using KEGGREST
Rscript ../../scripts/kegg_extract.r hsa05224
```

### Extract expression matrices

```bash
cd ~/data/rnaseq/rna_bc/info

for group in fpkm_nonnormalized qn
do
    cat breast_rnaseq_${group}.txt |
        tsv-join -H -f <(
            cat hsa05224.lst |
            sed '1isymbol'
            ) -k symbol |
        datamash transpose |
        sed 2,3d |
        perl -pe 's/symbol/cell_line/' \
        > hsa05224_${group}.tsv
done
```
