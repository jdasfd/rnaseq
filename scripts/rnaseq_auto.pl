#!/usr/bin/perl -w

#   rnaseq_auto.pl - automatically extract reads count from raw RNA-seq files
#
#
#   Author: Yuqian Jiang
#   Created: 2023-05-06
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 23/05/06: The initial version.


use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Path::Tiny;
use List::Util;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

rnaseq_auto.pl - automatically extract reads count from raw RNA-seq files

=head1 SYNOPSIS

    rnaseq_auto.pl (v1.0.0)
    Automatically extracting gene counts from raw RNA-seq files.

    Usage:
    perl rnaseq_auto.pl -i <filename> -g <genome> -a <annotation> [options]
    perl rnaseq_auto.pl -t PE -i <filename_1> -i <filename_2> -g <genome> -a <annotation> [options]


    Options:
    -i,--in             input file (also accept .gz format), required
    -g,--genome         reference genome for hisat2 alignment, required
    -a,--annotation     annotation files - both gtf and gff acceptable, required
    -w,--workdir        working directories, default is the current path
    -t,--type           RNA-seq in single end (SE) or paired end (PE) mode, default: SE
    --thread            threads for hisat2 and featureCount, default: 1
    -h,--help           help information

=cut

GetOptions(
    "i|in=s@"           => \(my $input),
    "g|genome=s"        => \(my $genome),
    "a|annotation=s"    => \(my $annotation),
    "w|workdir"         => \(my $workdir),
    "t|type=s"          => \(my $type = 'SE'),
    "thread=s"          => \(my $thread = '1'),
    "h|help"            => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $input) {
    print STDERR "Error: cannot find input files\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! @{$input} ) {
    print STDERR "Error: cannot find input files\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $workdir ) {
    $workdir = Path::Tiny -> cwd;
}

if ( $type ne "SE" && $type ne "PE") {
    print STDERR "Error: please correctly entered --type parameter\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $genome ) {
    print STDERR "Error: cannot find genome\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($genome) -> is_file ) {
    die "Error: cannot open file [$genome]";
}

if ( !defined $annotation ) {
    print STDERR "Error: cannot find annotation\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($annotation) -> is_file ) {
    die "Error: cannot open file [$annotation]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

# all global variables
my ($name, $inpath, $suffix, $out_name);
my @suffixlist = (".fastq.gz", ".fastq.bz2", ".fastq", ".fq.gz", ".fq.bz2", ".fq");
my @genomsuffix = ("fasta", "fa");
my @feature;

# mode selection
my $in_num = @{$input};
if ( $in_num == 1 ) {
    if ( $type eq "SE" ) {
        print STDERR "==> SE mode detected\n";
    }
    else {
        die "Error: PE mode with only an input file\n";
    }
}
elsif ( $in_num == 2 ) {
    if ( $type eq "PE" ) {
        print STDERR "==> PE mode detected\n";
    }
    else {
        die "Error: SE mode with two input files, error\n";
    }
}
else {
    die "Error: input exceed the range\n";
}

# check index
my ($genom_name, $genom_path, $genom_suffix) = fileparse ($genome, @genomsuffix);
my $genom_ref = $genom_path.$genom_name;
$genom_ref =~ s/\.$//;
if ( glob ("$genom_path"."$genom_name"."*.ht2") ) {
    print STDERR "==> Hisat2 index exists\n";
}
else {
    my $result_idx = &ht2_index ($genom_name, $genom_path, $genom_suffix, $thread);
    if ( $result_idx != 0 ) {
        die "Error: Hisat-index wrong, aborted\n";
    }
}

# check gtf
if ( $annotation =~ /\.gff$/ ) {
    print STDERR "==> Converting gff to gtf via gffread\n";
    my $gffpath = dirname ($annotation);
    my $gtf = "$gffpath"."/"."$genom_name".".gtf";
    system "gffread $annotation -T -o $gtf";
}
elsif ( $annotation =~ /\.gtf/ ) {
    print STDERR "==> Annotation gtf already exists\n";
}

# dealing with filepath and name
if ( $in_num == 1 ) {
    my $path_in = shift (@{$input});
    ($name, $inpath, $suffix) = fileparse ($path_in, @suffixlist);
    if ( $suffix eq "" ) {
        die "Error: input suffix wrong\n";
    }
}
elsif ( $in_num == 2 ) {
    my $path_in_1 = shift (@{$input});
    my $path_in_2 = shift (@{$input});
    my ($name_1, $inpath_1, $suffix_1) = fileparse ($path_in_1, @suffixlist);
    my ($name_2, $inpath_2, $suffix_2) = fileparse ($path_in_2, @suffixlist);
    if ( $suffix_1 eq "" || $suffix_2 eq "" ) {
        die "Error: input suffix wrong\n";
    }
    elsif ( $inpath_1 ne $inpath_2) {
        die "Error: input path not consistent\n";
    }
    else {
        my $n1 = $1 if $name_1 =~ /^(.+?)_\d/;
        my $n2 = $1 if $name_2 =~ /^(.+?)_\d/;
        if ( $n1 eq $n2 ) {
            $name = $n1;
            $inpath = $inpath_1;
            $suffix = $suffix_1;
        }
        else {
            die "Error: input name not consistent\n";
        }
    }
}
else {
    die "Error: input exceed the range\n";
}

# outdir
my $outdir = $workdir."/result";
if ( ! -d $outdir ) {
    mkdir $outdir;
}

$out_name = $outdir."/".$name;

#----------------------------------------------------------#
# main program
#----------------------------------------------------------#

# trim_galore
if ( $type eq "SE" ) {
    my $trim_in = $inpath.$name.$suffix;
    if ( ! glob $out_name."*_trimmed*" ) {
        my $result_trim = &trim_galore_SE($trim_in, $outdir);
        if ( $result_trim != 0 ) {
            die "Error: trim_galore wrong\n";
        }
    }
}
elsif ( $type eq "PE" ) {
    my $trim_in_1 = $inpath.$name."_1".$suffix;
    my $trim_in_2 = $inpath.$name."_2".$suffix;
    if ( ! glob $out_name."*_val*" ) {
        my $result_trim = &trim_galre_PE($trim_in_1,$trim_in_2, $outdir);
        if ( $result_trim != 0 ) {
            die "Error: trim_galore wrong\n";
        }
    }
}

# ht2-align
if ( $type eq "SE" ) {
    my $ht_in = $out_name."_trimmed.fq.gz";
    my $ht_out = $out_name.".tmp.sam";
    if ( ! glob $ht_out || ! glob $out_name.".sort.bam" ) {
        my $result_ht = &ht2_align_SE($ht_in, $ht_out, $genom_ref, $thread);
        if ( $result_ht != 0 ) {
            die "Error: hisat align wrong\n";
        }
    }
    else {
        print STDERR "==> Hisat2 align result exists\n";
    }
}
elsif ( $type eq "PE" ) {
    my $ht_in_1 = $out_name."_1_val_1.fq.gz";
    my $ht_in_2 = $out_name."_2_val_2.fq.gz";
    my $ht_out = $out_name.".tmp.sam";
    if ( ! glob $ht_out || ! glob $out_name.".sort.bam" ) {
        my $result_ht = &ht2_align_PE($ht_in_1, $ht_in_2, $ht_out, $genom_ref, $thread);
        if ( $result_ht != 0 ) {
            die "Error: hisat align wrong\n";
        }
    }
    else {
        print STDERR "==> Hisat2 align result exists\n";
    }
}

# sam2bam
my $bamfile = $out_name.".sort.bam";
if ( ! glob $bamfile ) {
    print STDERR "==> Sam to bam converting\n";
    my $samfile = $out_name.".tmp.sam";
    system "samtools sort -@ $thread $samfile > $bamfile";
    system "samtools index $bamfile";
    path($samfile) -> remove if glob $out_name.".tmp.sam";
}
else {
    print STDERR "==> Bam file exists\n";
}

# feature counts
my $countfile = $out_name.".tmp.tsv";
my $featurefile = $out_name.".count.tsv";
print STDERR "==> Counting reads via featureCounts\n";
system "featureCounts -T $thread -a $annotation -o $countfile $bamfile";

open my $TSV_IN, "<", $countfile;
while ( <$TSV_IN> ) {
    chomp;
    next if /^#/;
    my @array = split/\t/, $_;
    my $print = "$array[0]\t$array[6]\n";
    push (@feature, $print);
}
close $TSV_IN;

path("$featurefile") -> spew(@feature);

#----------------------------------------------------------#
# sub-program
#----------------------------------------------------------#

# trim_galore 4 cores is the best
sub trim_galore_SE {
    my ($file, $out) = @_;
    print STDERR "==> Trimming adapter for SE\n";
    my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 $file -o $out";
    return $result;
}

sub trim_galre_PE {
    my ($file_5, $file_3, $out) = @_;
    print STDERR "==> Trimming adapter for PE\n";
    my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 --paired $file_5 $file_3 -o $out";
    return $result;
}

sub ht2_index {
    my ($GENOM, $PATH, $SUFFIX, $THREAD) = @_;
    my $NAME = $1 if $GENOM =~ /^(.+?)\./;
    my $OLD = "$PATH"."$GENOM"."$SUFFIX";
    my $NEW = "$PATH"."$NAME";
    print STDERR "==> Hisat2 genome indexing\n";
    my $result = system "hisat2-build -p $THREAD $OLD $NEW";
    return $result;
}

sub ht2_align_SE {
    my ($HT_IN, $HT_OUT, $REF, $THREAD) = @_;
    print STDERR "==> Hisat2 alignment SE mode\n";
    my $result = system "hisat2 -p $THREAD -x $REF -U $HT_IN -S $HT_OUT";
    return $result;
}

sub ht2_align_PE {
    my ($HT_IN_1, $HT_IN_2, $HT_OUT, $REF, $THREAD) = @_;
    print STDERR "==> Hisat2 alignment PE mode\n";
    my $result = system "hisat2 -p $THREAD -x $REF -1 $HT_IN_1 -2 $HT_IN_2 -S $HT_OUT";
    return $result;
}
