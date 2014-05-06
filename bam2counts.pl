#!/usr/bin/env perl
use strict;
use warnings;
use v5.10.0;
use autodie;
use File::Basename;

# TODO: Output mapped/unmapped counts summary

my $out_file  = shift;
my @bam_files = @ARGV;

my $usage = "USAGE: $0 output.txt sample1.bam sample2.bam sample3.bam\n";
die $usage if defined $out_file && $out_file =~ /\.bam$/i;
die $usage if scalar @bam_files < 1;

my @header;
open my $header_fh, "-|", "samtools view -H $bam_files[0]";
while (<$header_fh>) {
    next unless /^\@SQ/;
    my ($seqid) = /SN:(\S+)/;
    push @header, $seqid;
}

my %counts;
my @samples;
for my $file (@bam_files) {
    say "Processing $file";
    my $sample_id = fileparse $file, ".bam";
    push @samples, $sample_id;
    open my $bam_fh, "-|", "samtools view $file";
    while (my $read = <$bam_fh>) {
        my $seqid = (split /\t/, $read)[2];
        $counts{$seqid}{$sample_id}++;
    }
}

open my $out_fh, ">", $out_file;
say $out_fh join "\t", @samples;
for my $seqid (sort @header) {
    my @seqid_counts = @{ $counts{$seqid} }{@samples};
    $_ //= 0 for @seqid_counts;
    say $out_fh join "\t", $seqid, @seqid_counts;
}
close $out_fh;
