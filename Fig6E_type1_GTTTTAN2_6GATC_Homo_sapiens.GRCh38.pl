#!/usr/bin/perl
use strict;
use warnings;

##########################################
# Configuration
##########################################
my $gtf_file     = "Homo_sapiens.GRCh38.113.chr.gtf"; # Input GTF file
my $ref_genome   = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"; # Reference genome file
my $output_file  = "Results_Fig6E_type1_GTTTTAN2_6GATC_Homo_sapiens.GRCh38.txt"; # Output file

##########################################
# Processing
##########################################
open my $gtf_fh, '<', $gtf_file or die "Cannot open $gtf_file: $!";

# Hash to store results; key = MatchedSequence + Location
my %results;
my $target_pattern = 'GTTTTA'; # Fixed portion of the search pattern

while (<$gtf_fh>) {
    chomp;
    next if /^#/;

    # Parse GTF fields
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attr_str) = split /\t/;
    next unless $feature eq 'gene'; # Only process entries with the "gene" feature

    # Parse attributes from the GTF record
    my %attr = parse_attributes($attr_str);
    my $gene_id    = $attr{gene_id} // 'NA';
    my $gene_name  = $attr{gene_name} // $gene_id;
    my $biotype    = $attr{gene_biotype} // 'NA';

    ##########################################
    # Extract sequence with 20 bp flanks
    ##########################################
    my $extract_start = $start - 20;
    $extract_start = 1 if $extract_start < 1; # Ensure the start position is not less than 1
    my $extract_end = $end + 20;
    my $region = "$chr:$extract_start-$extract_end";

    # Retrieve the sequence for the region using samtools
    my $seq = `samtools faidx $ref_genome $region 2>/dev/null`;
    $seq =~ s/^>.*\n//; # Remove FASTA header
    $seq =~ s/\n//g;    # Remove line breaks

    # Reverse complement if the strand is negative
    $seq = reverse_complement($seq) if $strand eq '-';

    ##########################################
    # Search for GTTTTA(N2-6)GCTA pattern
    ##########################################
    my $trimmed_seq = substr($seq, 20, length($seq) - 40); # Trim flanking sequences

    while ($trimmed_seq =~ /(GTTTTA([ACGT]{2,6})GCTA)/g) {
        my $matched_seq = $1;       # Full matched sequence
        my $variable_part = $2;     # Variable portion (N2-6)

        # Match start and end positions within the trimmed sequence (1-based)
        my $match_start = pos($trimmed_seq) - length($matched_seq) + 1;
        my $match_end   = pos($trimmed_seq);

        # Convert to genomic coordinates
        my ($genomic_start, $genomic_end);
        if ($strand eq '+') {
            $genomic_start = $start + ($match_start - 1);
            $genomic_end   = $start + ($match_end - 1);
        } else {
            $genomic_end   = $end - ($match_start - 1);
            $genomic_start = $end - ($match_end - 1);
        }

        # Retrieve upstream 20 bp sequence
        my $full_match_start = 20 + $match_start;
        my $upstream_start   = $full_match_start - 20;
        $upstream_start = 0 if $upstream_start < 0;
        my $upstream20bp = substr($seq, $upstream_start-1, 20);

        # Use MatchedSequence + Location as a unique key
        my $location = "$chr:$genomic_start-$genomic_end";
        my $key = "$matched_seq|$location";

        # Initialize result entry if the key does not exist
        if (!exists $results{$key}) {
            $results{$key} = {
                Gene_location => {},
                Strand        => $strand,
                GeneID        => {},
                GeneName      => {},
                GeneBiotype   => {},
                MatchedSeq    => $matched_seq,
                Location      => $location,
                VariablePart  => $variable_part,
                Upstream20bp  => $upstream20bp,
            };
        }

        # Store information in a de-duplicated manner
        $results{$key}->{Gene_location}{"$chr:$start-$end"} = 1;
        $results{$key}->{GeneID}{$gene_id} = 1;
        $results{$key}->{GeneName}{$gene_name} = 1;
        $results{$key}->{GeneBiotype}{$biotype} = 1;
    }
}
close $gtf_fh;

##########################################
# Write results to output file
##########################################
open my $out_fh, '>', $output_file or die "Cannot write to $output_file: $!";

# Print header
print $out_fh join("\t", qw(
    Gene_location Strand GeneID GeneName GeneBiotype MatchedSequence Location VariablePart Upstream20bp
)), "\n";

# Write each result entry
foreach my $key (sort keys %results) {
    my $r = $results{$key};
    print $out_fh join("\t",
        join(";", sort keys %{$r->{Gene_location}}),
        $r->{Strand},
        join(";", sort keys %{$r->{GeneID}}),
        join(";", sort keys %{$r->{GeneName}}),
        join(";", sort keys %{$r->{GeneBiotype}}),
        $r->{MatchedSeq},
        $r->{Location},
        $r->{VariablePart},
        $r->{Upstream20bp},
    ), "\n";
}
close $out_fh;

print "Pipeline completed successfully! Results saved to $output_file\n";

##########################################
# Subroutines
##########################################

# Parse attribute column from GTF file
sub parse_attributes {
    my $attr_str = shift;
    my %attrs;
    foreach my $pair (split /;/, $attr_str) {
        $pair =~ s/^\s+|\s+$//g;           # Remove leading/trailing whitespace
        next unless $pair;                 # Skip empty pairs
        my ($key, $value) = split /\s+/, $pair, 2;
        $value =~ s/^"|"$//g;              # Remove surrounding quotes
        $attrs{$key} = $value;
    }
    return %attrs;
}

# Generate reverse complement of a DNA sequence
sub reverse_complement {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/; # Complement the bases
    return $seq;
}
