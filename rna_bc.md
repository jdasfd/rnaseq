# Breast cancer RNA-seq data

```bash
git clone https://github.com/wang-q/fig_table.git
brew install dos2unix
```

Article link: [Cell, 2016, 164: 293-309](https://www.cell.com/cell/fulltext/S0092-8674(15)01624-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415016244%3Fshowall%3Dtrue)

[Convert CR 2 LF](https://toolslick.com/conversion/text/new-line).

```bash
mkdir -p ~/data/rna_bc/article_SM
cd ~/data/rna_bc/article_SM

# download supplementary tables into the dir

# extract cell and its subtype info
cat cell_line_subtypes.tsv |
    tsv-select -H -f cell_line,subtype_neve \
    > cell_anno.tsv

#perl ~/Scripts/fig_table/xlsx2csv.pl -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3B | sed '1,4d' | tsv-select -d ',' -f 1,5,9,13 | sed 1d | sed '1ibasal_A,basal_B,her,luminal' | mlr --icsv --otsv cat > gene_neve.tsv

# extract each gene pathway
# basal_A
perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3H |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Pathway\ id" |
    sed 1d \
    > basal_A_pathway.lst

perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3H |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Genes" |
    sed 1d |
    tr "," "\n" |
    sort |
    uniq \
    > basal_A_gene.lst

# basal_B
perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3I |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Pathway\ id" |
    sed 1d \
    > basal_B_pathway.lst

perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3I |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Genes" |
    sed 1d |
    tr "," "\n" |
    sort |
    uniq \
    > basal_B_gene.lst

# HER2
perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3J |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Pathway\ id" |
    sed 1d \
    > HER2_pathway.lst

perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3J |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Genes" |
    sed 1d |
    tr "," "\n" |
    sort |
    uniq \
    > HER2_gene.lst

# Luminal
perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3K |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Pathway\ id" |
    sed 1d \
    > Luminal_pathway.lst

perl ~/Scripts/fig_table/xlsx2csv.pl \
    -f 1-s2.0-S0092867415016244-mmc4.xlsx --sheet S3K |
    sed '1,2d' |
    mlr --icsv --otsv cat |
    tsv-select -H -f "Genes" |
    sed 1d |
    tr "," "\n" |
    sort |
    uniq \
    > Luminal_gene.lst

for group in basal_A basal_B HER2 Luminal
do
    # combine into pathway
    cat ${group}_pathway.lst |
        awk -v GR=$group '{print (GR"\t"$0)}' \
        >> pathway.tsv
    
    # combine into gene pathway
    cat ${group}_gene.lst |
        awk -v GR=$group '{print (GR"\t"$0)}' \
        >> gene.tsv
done

wc -l *_pathway.lst
#  264 HER2_pathway.lst
#   69 Luminal_pathway.lst
#  580 basal_A_pathway.lst
#   90 basal_B_pathway.lst
# 1003 total

cat pathway.tsv | cut -f 2 | sort | uniq | wc -l
#825
# pathway info has the repeated part

wc -l *_gene.lst
# 139 HER2_gene.lst
#  83 Luminal_gene.lst
# 238 basal_A_gene.lst
# 104 basal_B_gene.lst
# 564 total

cat gene.tsv | cut -f 2 | sort | uniq | wc -l
#544
# about 20 genes repeated into each other

# extract related gene expression matrices
cat breast_rnaseq_fpkm_nonnormalized.txt |
    tsv-select -H -e gene_id,ensembl_id |
    tsv-join -H -f <(
        cat gene.tsv |
        cut -f 2 |
        sort |
        uniq |
        sed '1isymbol'
        ) -k symbol \
    > matrix_exp.tsv

Rscript -e '
    library("KEGGREST")
    library("EnrichmentBrowser")
    hsapathway <- downloadPathways("hsa")
    hsa <- getGenesets(org = "hsa", db = "kegg", gene.id.type = "SYMBOL",cache = TRUE, return.type="list")
    genelist <- hsa$hsa05224
    write.table(genelist, sep = "\t", file = "hsa05224.tsv")
'

cat hsa05224.tsv |
    cut -f 2 |
    sed 1d |
    perl -nle 'print "$1" if /^"(.+)"$/;' \
    > hsa05224_gene.lst

cat breast_rnaseq_fpkm_nonnormalized.txt |
    tsv-join -H -f <(
        cat hsa05224_gene.lst |
        sed '1isymbol'
        ) -k symbol \
    > hsa05224_exp.tsv
```