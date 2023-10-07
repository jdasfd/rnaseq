# The study of AM RNA-seq

## Raw RNA-seq data

A means the AM (treatment). B means the blank.

```bash
mkdir -p ~/data/rnaseq/rna_AM/SEQ
# manually upload all files into SEQ dir
```

## Get AM genomes

Two version of AM genomes.

One is from the [MycoCosm](https://mycocosm.jgi.doe.gov/Rhiir2_1/Rhiir2_1.home.html). Another is from a recently [research](https://doi.org/10.1093/g3journal/jkad077).

```bash
mkdir -p ~/data/rnaseq/rna_AM/GENOME
cd ~/data/rnaseq/rna_AM/GENOME

# mycocosm
mkdir Rhiir_v2_1
# download from the jgi mycocosm
mv Rhiir2_1_GeneCatalog_proteins_20160502.aa.fasta Rhiir2_1.pep
mv Rhiir2_1_GeneCatalog_CDS_20160502.fasta Rhiir2_1.cds
mv Rhiir2_1_GeneCatalog_20160502.gff3 Rhiir2_1.gff3
mv Rhiir2_1_AssemblyScaffolds_Repeatmasked.fasta Rhiir2_1.fa

# a recent article
mkdir Rir_DAOM197198_2023
# if blocked, please download them on your own
# wget https://zenodo.org/record/7713976/files/Rhizophagus_irregularis_DAOM197198.fa?download=1
# wget https://zenodo.org/record/7713976/files/Rhizophagus_irregularis_DAOM197198_Illumina%2BONT_curated.gff3?download=1
# wget https://zenodo.org/record/7713976/files/Rhizophagus_irregularis_DAOM197198_proteins_Illumina%2BONT_curated.fa?download=1
mv Rhizophagus_irregularis_DAOM197198.fa Rir_DAOM197198_2023.fa
mv Rhizophagus_irregularis_DAOM197198.gff3 Rir_DAOM197198_2023.gff3
mv Rhizophagus_irregularis_DAOM197198_proteins_Illumina+ONT_curated.fa Rir_DAOM197198_2023.pep
```

## Deal with RNA-seq

Using `rnaseq_auto.pl` to analyse different raw seq data.

```bash
cd ~/data/rnaseq/rna_AM

for group in A B
do
    for i in {1..3}
    do
        perl ../scripts/rnaseq_auto.pl -t PE -i SEQ/${group}${i}_1.fq.gz -i SEQ/${group}${i}_2.fq.gz \
            -g GENOME/Rir_DAOM197198_2023/Rir_DAOM197198_2023.fa \
            -a GENOME/Rir_DAOM197198_2023/Rir_DAOM197198_2023.gff3 \
            -w result/Rir_DAOM197198_${group} \
            --thread 12
    done
done
```
