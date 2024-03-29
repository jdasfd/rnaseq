#!/usr/bin/perl -w

#   rnaseq_auto.pl - automatically extract reads count from raw RNA-seq files
#
#
#   Author: Yuqian Jiang
#   Created: 2023-05-06
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 23-05-06: The initial version.
#   Version 1.0.1 23-07-25: Bug fixes: gff2gtf name error;
#                                      mkdir error reported by trim-galore.
#                                      genome suffix wrong for ht2-index error.
#                                      glob error so realize by Path::Tiny.
#   Version 1.1.0 23-07-25: Add new parameter: -attribute for featureCounts.
#   Version 1.1.1 23-09-25: Bug fixes: mkdir not work; gff3 not recognized.
#                           Add: log files to record all processes.
#   Version 1.1.2 23-09-26: Bug fixes: log to avoid covering old log.
#                           Add: local time to record each task.
#   Version 1.1.3 23-09-28: Bug fixes: gtf and gff/gff3 conversion but featureCount failed.
#   Version 1.2.1 24-03-03: Bug fixes: the whole name before _1/_2 PE will be treated as the name.
#                           Add: judge the different mode for featureCounts.

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Path::Tiny;
use List::Util;
use IO::Tee;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

rnaseq_auto.pl - automatically extract reads count from raw RNA-seq files

=head1 SYNOPSIS

    rnaseq_auto.pl (v1.2.1)
    Automatically extracting gene counts from raw RNA-seq files.

    Usage:
    SE mode:
    perl rnaseq_auto.pl -i <filename> -g <genome> -a <annotation> [options]
    PE mode:
    perl rnaseq_auto.pl -t PE -i <filename_1> -i <filename_2> -g <genome> -a <annotation> [options]

    Required:
    -i,--in          STR       input file (also accept .gz format), required
    -g,--genome      STR       reference genome for hisat2 alignment, required
    -a,--annotation  STR       annotation files - both gtf and gff acceptable, required

    Options:
    -w,--workdir     STR       working directories, default is the current path
    -t,--type        [SE|PE]   RNA-seq in single end (SE) or paired end (PE) mode, default: SE
    --thread         INT       threads for hisat2 and featureCount, default: 1
    --attribute      STR       attribute for featureCounts, default: gene_id
    -h,--help                  help information

=cut

GetOptions(
    "i|in=s@"           => \(my $input),
    "g|genome=s"        => \(my $genome),
    "a|annotation=s"    => \(my $annotation),
    "w|workdir=s"       => \(my $workdir),
    "attribute=s"       => \(my $attribute = 'gene_id'),
    "t|type=s"          => \(my $type = 'SE'),
    "thread=i"          => \(my $thread = '1'),
    "h|help"            => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: cannot find input file.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! @{$input} ) {
    print STDERR "Error: cannot find input files.\n";
    die Getopt::Long::HelpMessage(1);
}

if ( $type ne "SE" && $type ne "PE") {
    print STDERR "Error: please choose correct mode.\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $genome ) {
    print STDERR "Error: cannot find genome.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($genome) -> is_file ) {
    die "Error: cannot open file [$genome].";
}

if ( !defined $annotation ) {
    print STDERR "Error: cannot find annotation.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($annotation) -> is_file ) {
    die "Error: cannot open file [$annotation].";
}

if ( !defined $workdir ) {
    my $current = Path::Tiny -> cwd;
    $workdir = $current."/result";
}
else {
    if ( $workdir =~ /\/$/ ) {
        $workdir =~ s/\/$//;
    }
    path ($workdir) -> mkdir;
}

unless ( path("$workdir/log.txt") -> is_file ) {
    path("$workdir/log.txt") -> touch;
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

# file output handle
my $gtf;
my $tee_add = IO::Tee -> new ( ">> $workdir/log.txt", \*STDERR );

# all global variables
my ($name, $inpath, $suffix, $out_name);
my @suffixlist = (".fastq.gz", ".fastq.bz2", ".fastq", ".fq.gz", ".fq.bz2", ".fq");
my @genomsuffix = (".fasta", ".fa");
my @feature;

# mode selection
my $in_num = @{$input};
if ( $in_num == 1 ) {
    if ( $type eq "SE" ) {
        my $time = localtime;
        print $tee_add "RNA-seq start at $time, SE mode detected.\n";
        print $tee_add "Genome:$genome.\n";
        my $in_file = shift (@{$input});
        print $tee_add "Input:$in_file.\n";
        print $tee_add "\n";
    }
    else {
        print $tee_add "Error: PE mode with only an input file.\n";
        die;
    }
}
elsif ( $in_num == 2 ) {
    if ( $type eq "PE" ) {
        my $time = localtime;
        print $tee_add "RNA-seq start at $time, PE mode detected.\n";
        print $tee_add "Genome:$genome.\n";
        for (@{$input}) {
            print $tee_add "Input: $_.\n";
        }
        print $tee_add "\n";
    }
    else {
        print $tee_add "Error: SE mode with two input files, error.\n";
        die;
    }
}
else {
    print $tee_add "Error: input exceed the range.\n";
    die;
}

# dealing with filepath and name
if ( $in_num == 1 ) {
    my $path_in = shift (@{$input});
    ($name, $inpath, $suffix) = fileparse ($path_in, @suffixlist);
    if ( $suffix eq "" ) {
        print $tee_add "Error: input suffix wrong.\n";
        die;
    }
}
elsif ( $in_num == 2 ) {
    my $path_in_1 = shift (@{$input});
    my $path_in_2 = shift (@{$input});
    my ($name_1, $inpath_1, $suffix_1) = fileparse ($path_in_1, @suffixlist);
    my ($name_2, $inpath_2, $suffix_2) = fileparse ($path_in_2, @suffixlist);
    if ( $suffix_1 eq "" || $suffix_2 eq "" ) {
        print $tee_add "Error: input suffix wrong.\n";
        die;
    }
    elsif ( $inpath_1 ne $inpath_2) {
        print $tee_add "Error: input path not consistent.\n";
        die;
    }
    else {
        my $n1 = $1 if $name_1 =~ /^(.+)_[1|2]$/;
        my $n2 = $1 if $name_2 =~ /^(.+)_[1|2]$/;
        if ( $n1 eq $n2 ) {
            $name = $n1;
            $inpath = $inpath_1;
            $suffix = $suffix_1;
        }
        else {
            print $tee_add "Error: input name not consistent.\n";
            die;
        }
    }
}
else {
    print $tee_add "Error: input exceed the range.\n";
    die;
}

$out_name = $workdir."/".$name;
my $samtmp = $out_name.".tmp.sam";
my $bamfile = $out_name.".sort.bam";
my $counttmp = $out_name.".tmp.tsv";
my $featurefile = $out_name.".count.tsv";

# check index
my ($genom_name, $genom_path, $genom_suffix) = fileparse ($genome, @genomsuffix);
my $genom_ref = $genom_path.$genom_name;
if ( glob ("$genom_path"."$genom_name"."*.ht2") ) {
    print $tee_add "==> Hisat2 index exists\n";
}
else {
    print $tee_add "==> Hisat2 genome indexing\n";
    my $result_idx = &ht2_index ($genom_name, $genom_path, $genom_suffix, $thread);
    if ( $result_idx != 0 ) {
        print $tee_add "Error: Hisat-index wrong, aborted.\n";
        die;
    }
}

# check gtf
if ( $annotation =~ /\.gff3?$/ ) {
    print $tee_add "==> Converting gff to gtf via gffread\n";
    my $gffpath = dirname ($annotation);
    $gtf = "$gffpath"."/"."$genom_name".".gtf";
    system "gffread $annotation -T -o $gtf";
}
elsif ( $annotation =~ /\.gtf/ ) {
    print $tee_add "==> Annotation gtf already exists\n";
    $gtf = $annotation;
}

#----------------------------------------------------------#
# main program
#----------------------------------------------------------#

# trim_galore
if ( $type eq "SE" ) {
    my $trim_in = $inpath.$name.$suffix;
    if ( ! glob $out_name."*_trimmed*" ) {
        print $tee_add "==> Trimming adapter for SE\n";
        my $result_trim = &trim_galore_SE($trim_in, $workdir);
        if ( $result_trim != 0 ) {
            print $tee_add "Error: trim_galore wrong.\n";
            die;
        }
    }
    else {
        print $tee_add "Already trimmed.\n";
    }
}
elsif ( $type eq "PE" ) {
    my $trim_in_1 = $inpath.$name."_1".$suffix;
    my $trim_in_2 = $inpath.$name."_2".$suffix;
    if ( ! glob $out_name."*_val*" ) {
        print $tee_add "==> Trimming adapter for PE\n";
        my $result_trim = &trim_galre_PE($trim_in_1,$trim_in_2, $workdir);
        if ( $result_trim != 0 ) {
            print $tee_add "Error: trim_galore wrong\n";
            die;
        }
    }
    else {
        print $tee_add "Already trimmed.\n";
    }
}

# ht2-align
if ( $type eq "SE" ) {
    my $ht_in = $out_name."_trimmed.fq.gz";
    if ( ! path($ht_in) -> is_file ) {
        print $tee_add "Error occured in trim_galore part.\n";
        die;
    }
    else {
        if ( path($bamfile) -> is_file || path($samtmp) -> is_file ) {
            print $tee_add "Hisat2 align result exists.\n";
        }
        else {
            print $tee_add "==> Hisat2 alignment SE mode\n";
            my $result_ht = &ht2_align_SE($ht_in, $samtmp, $genom_ref, $thread);
            if ( $result_ht != 0 ) {
                print $tee_add "Error: hisat align wrong.\n";
                die;
            }
        }
    }
}
elsif ( $type eq "PE" ) {
    my $ht_in_1 = $out_name."_1_val_1.fq.gz";
    my $ht_in_2 = $out_name."_2_val_2.fq.gz";
    if ( path($ht_in_1) -> is_file && path($ht_in_2) -> is_file ) {
        if ( path($samtmp) -> is_file || path($bamfile) -> is_file ) {
            print $tee_add "Hisat2 align result exists.\n";
        }
        else {
            print $tee_add "==> Hisat2 alignment PE mode\n";
            my $result_ht = &ht2_align_PE($ht_in_1, $ht_in_2, $samtmp, $genom_ref, $thread);
            if ( $result_ht != 0 ) {
                print $tee_add "Error: hisat align wrong.\n";
                die;
            }
        }
    }
    else {
        print $tee_add "Error occured in trim_galore part.\n";
        die;
    }
}

# sam2bam
if ( !path($bamfile) -> is_file ) {
    print $tee_add "==> Sam to bam converting\n";

    system "samtools sort -@ $thread $samtmp > $bamfile";
    system "samtools index $bamfile";
    path($samtmp) -> remove if path($samtmp) -> is_file;
}
else {
    print $tee_add "==> Bam file exists\n";
}

# feature counts
print $tee_add "==> Counting reads via featureCounts\n";
if ( $type eq "SE" ) {
    system "featureCounts -T $thread -a $gtf -g $attribute -o $counttmp $bamfile";
}
elsif ( $type eq "PE" ) {
    system "featureCounts -p -T $thread -a $gtf -g $attribute -o $counttmp $bamfile";
}

if ( !path($counttmp) -> is_file ) {
    print $tee_add "featureCounts Error!\n";
    die;
}
else {
    open my $TSV_IN, "<", $counttmp;
    while ( <$TSV_IN> ) {
        chomp;
        next if /^#/;
        my @array = split/\t/, $_;
        my $print = "$array[0]\t$array[6]\n";
        push (@feature, $print);
    }
    close $TSV_IN;
}

path("$featurefile") -> spew(@feature);

# final
my $time = localtime;
print $tee_add "All processes completed at $time!\n";

#----------------------------------------------------------#
# sub-program
#----------------------------------------------------------#

# trim_galore 4 cores is the best
sub trim_galore_SE {
    my ($file, $out) = @_;

    my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 $file -o $out";
    return $result;
}

sub trim_galre_PE {
    my ($file_5, $file_3, $out) = @_;

    my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 --paired $file_5 $file_3 -o $out";
    return $result;
}

sub ht2_index {
    my ($GENOM, $PATH, $SUFFIX, $THREAD) = @_;

    my $OLD = "$PATH"."$GENOM"."$SUFFIX";
    my $NEW = "$PATH"."$GENOM";

    my $result = system "hisat2-build -p $THREAD $OLD $NEW";
    return $result;
}

sub ht2_align_SE {
    my ($HT_IN, $HT_OUT, $REF, $THREAD) = @_;

    my $result = system "hisat2 -p $THREAD -x $REF -U $HT_IN -S $HT_OUT";
    return $result;
}

sub ht2_align_PE {
    my ($HT_IN_1, $HT_IN_2, $HT_OUT, $REF, $THREAD) = @_;

    my $result = system "hisat2 -p $THREAD -x $REF -1 $HT_IN_1 -2 $HT_IN_2 -S $HT_OUT";
    return $result;
}

=head1 VERSION

1.2.1

=head1 AUTHORS

Yuqian Jiang, yuqian_j@outlook.com

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2023 by Yuqian Jiang.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut
