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

# all global variables
my $result;
my ($name, $inpath);
my @suffixlist = (".fastq.gz", ".fastq.bz2", ".fastq", ".fq.gz", ".fq.bz2", ".fq");

GetOptions(
    "i|in=s@"       => \(my $input),
    "w|workdir"     => \(my $workdir),
    "t|type=s"      => \(my $type = 'SE'),
    "thread=s"      => \(my $thread = '1'),
    "h|help"        => sub { Getopt::Long::HelpMessage(0) },
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

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

# mode selection
my $in_num = @{$input};
if ( $in_num == 1 ) {
    if ( $type eq "SE" ) {
        print STDERR "==> SE mode detected\n";
    }
    else {
        print STDERR "PE mode with only an input file, error\n";
        exit(1);
    }
}
elsif ( $in_num == 2 ) {
    if ( $type eq "PE" ) {
        print STDERR "==> PE mode detected\n";
    }
    else {
        print STDERR "SE mode with two input files, error\n";
    }
}
else {
    print STDERR "Error: input exceed the range\n";
    exit(1);
}

# dealing with filepath and name
if ( $in_num == 1 ) {
    my $path_in = shift (@{$input});
    ($name, $inpath, my $suffix) = fileparse ($path_in, @suffixlist);
    if ( $suffix eq "" ) {
        print STDERR "Error: input suffix wrong\n";
        exit(1);
    }
}
elsif ( $in_num == 2 ) {
    my $path_in_1 = shift (@{$input});
    my $path_in_2 = shift (@{$input});
    my ($name_1, $inpath_1, $suffix) = fileparse ($path_in_1, @suffixlist);
}
else {
    print STDERR "Error: input exceed the range\n";
    exit(1);
}

# outdir
my $outdir = $workdir."/result";
if ( ! -d $outdir ) {
    mkdir $outdir;
}

#----------------------------------------------------------#
# main program
#----------------------------------------------------------#

=pod
# trim_galore
if ( $type eq "SE" ) {
    $result = &trim_galore_SE($filepath);
}
elsif ( $type eq "PE" ) {
    $result = &trim_galre_PE($filepath_1,$filepath_2);
}

if ( $result != 0 ) {
    print STDERR "Error: trim_galore wrong\n";
    exit(1);
}
=cut

=pod
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
=cut
