#!/usr/bin/perl
use strict;
use warnings;

# Check if input file is provided
die "Usage: $0 <reference.fa>\n" unless @ARGV == 1;

# 1. Read FASTA file and store in hash
my %chromosomes;
my $current_chr;
my $sequence = '';
my %hash;  # Will store 20bp sequences and their locations

open my $fh, '<', $ARGV[0] or die "Cannot open file $ARGV[0]: $!";
while (my $line = <$fh>) {
    chomp $line;
    
    # Handle chromosome headers
    if ($line =~ /^>(\S+)/) {
        if ($current_chr) {
            $chromosomes{$current_chr} = $sequence;
            $sequence = '';
        }
        $current_chr = $1;
    } else {
        # Store sequence in uppercase, including Ns
        $sequence .= uc($line);
    }
}
close $fh;

# Store the last chromosome
if ($current_chr) {
    $chromosomes{$current_chr} = $sequence;
}

# 2. Process each chromosome to find 23bp sequences
foreach my $chr (sort keys %chromosomes) {
    my $seq = $chromosomes{$chr};
    my $len = length($seq);
    
    # Slide through sequence with 23bp window
    for (my $pos = 0; $pos <= $len - 23; $pos++) {
        my $subseq = substr($seq, $pos, 23);
        
        # Skip if sequence contains N
        next if $subseq =~ /N/;
        
        # Case 1: Starts with CC but doesn't end with GG
        if ($subseq =~ /^CC/ && $subseq !~ /GG$/) {
            my $true_seq = substr($subseq, 3, 20);
            my $location = "$chr:".($pos+4)."-".($pos+23).":-";
            update_hash(\%hash, $true_seq, $location);
        }
        # Case 2: Doesn't start with CC but ends with GG
        elsif ($subseq !~ /^CC/ && $subseq =~ /GG$/) {
            my $true_seq = substr($subseq, 0, 20);
            my $location = "$chr:".($pos+1)."-".($pos+20).":+";
            update_hash(\%hash, $true_seq, $location);
        }
        # Case 3: Both starts with CC and ends with GG
        elsif ($subseq =~ /^CC/ && $subseq =~ /GG$/) {
            my $true_seq1 = substr($subseq, 0, 20);
            my $true_seq2 = substr($subseq, 3, 20);
            
            my $location1 = "$chr:".($pos+1)."-".($pos+20).":+";
            my $location2 = "$chr:".($pos+4)."-".($pos+23).":+";
            
            update_hash(\%hash, $true_seq1, $location1);
            update_hash(\%hash, $true_seq2, $location2);
        }
    }
}

# Helper subroutine to update the hash
sub update_hash {
    my ($hash_ref, $seq, $loc) = @_;
    
    if (exists $hash_ref->{$seq}) {
        $hash_ref->{$seq}[0] .= ";$loc";
        $hash_ref->{$seq}[1]++;
    } else {
        $hash_ref->{$seq}[0] = $loc;
        $hash_ref->{$seq}[1] = 1;
    }
}

# Print results sorted by sequence
foreach my $key (sort keys %hash) {
    print "$key\t$hash{$key}[0]\t$hash{$key}[1]\n";
}
