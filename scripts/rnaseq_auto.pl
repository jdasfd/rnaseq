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
    perl rnaseq_auto.pl -i <filename>
    perl rnaseq_auto.pl -i <filename_1> -i <filename_2> -t PE


    Options:
    -i,--in         input file (also accept .gz format), required
    -w,--workdir    working directories, default is the current path
    -t,--type       RNA-seq in single end (SE) or paired end (PE) mode, default: SE
    -h,--help       help information

=cut

GetOptions(
    "i|in=s@"           => \(my $input),
    "w|workdir"         => \(my $workdir),
    "t|type=s"          => \(my $type = 'SE'),
    "thread=s"          => \(my $thread = '1'),
    "g|genome=s"        => \(my $genome),
    "a|annotation=s"    => \(my $annotation),
    "h|help"            => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( ! @{$input} ) {
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
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($genome) -> is_file ) {
    die "Error: can't find file [$genome]";
}

if ( !defined $annotation ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( ! path($annotation) -> is_file ) {
    die "Error: can't find file [$annotation]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

# all global variables
my $result;
my ($name, $inpath, $suffix);
my @suffixlist = (".fastq.gz", ".fastq.bz2", ".fastq", ".fq.gz", ".fq.bz2", ".fq");
my @genomsuffix = ("fasta", "fa");

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
}

#----------------------------------------------------------#
# main program
#----------------------------------------------------------#

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

# trim_galore
if ( $type eq "SE" ) {
    my $trim_in = $inpath.$name.$suffix;
    $result = &trim_galore_SE($trim_in);
}
elsif ( $type eq "PE" ) {
    my $trim_in_1 = $inpath.$name."_1".$suffix;
    my $trim_in_2 = $inpath.$name."_2".$suffix;
    $result = &trim_galre_PE($trim_in_1,$trim_in_2);
}

if ( $result != 0 ) {
    die "Error: trim_galore wrong\n";
}


#----------------------------------------------------------#
# sub-program
#----------------------------------------------------------#

# trim_galore 4 cores is the best
sub trim_galore_SE {
    my $file = shift;
    print STDERR "==> Trimming adapter for SE\n";
    my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 $file -o $outdir";
    return $result;
}

sub trim_galre_PE {
    my ($file_5, $file_3) = @_;
    print STDERR "==> Trimming adapter for PE\n";
    print "$file_5, $file_3\n";
    #my $result = system "trim_galore -j 4 -q 30 --fastqc --length 20 --paired $file_5 $file_3 -o $outdir";
    return $result;
}

sub featureCounts {
    my ($annotation, $thread, $name) = @_;
#    featureCounts -T 10 -t gene -G Name
}
