#!/usr/bin/perl
use strict;
use warnings;

# Configuration parameters
my $gtf_file       = 'Homo_sapiens.GRCh38.113.chr.gtf'; # Input GTF file
my $ref_genome     = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'; # Reference genome file
my $target_pattern = 'GTTTTAGAGCTA'; # Target sequence
my $flanking_size  = 20; # Flanking region size (upstream and downstream)
my $max_mismatches = 2; # Maximum allowed mismatches
my $output_file    = 'Results_Fig6E_type2_GTTTTAGAGCTA_GA_fixed_02matches_Homo_sapiens.GRCh38.txt'; # Output file

# Check if samtools is available
my $samtools_path = `which samtools` or die "samtools is not found in PATH\n";
chomp $samtools_path;

# Step 1: Extract gene information from the GTF file
print "Step 1: Extracting and parsing gene information from the GTF file...\n";

# Store gene information
my @genes;

open my $gtf_fh, '<', $gtf_file or die "Cannot open $gtf_file: $!";
while (<$gtf_fh>) {
    next if /^#/; # Skip comment lines
    chomp;
    my @fields = split /\t/;
    next unless $fields[2] eq 'gene'; # Only process gene entries

    # Extract chromosome, feature type, start, end, and strand
    my ($chr, $type, $start, $end, $strand) = @fields[0,2,3,4,6];

    # Parse attributes from the ninth column
    my $attributes = $fields[8];
    my %attrs;
    while ($attributes =~ /(\w+)\s+"([^"]+)";/g) {
        $attrs{$1} = $2;
    }

    # Retrieve gene_id, gene_name, and gene_biotype
    my $gene_id     = $attrs{gene_id} || 'NA';
    my $gene_name   = $attrs{gene_name} || $gene_id; # Use gene_id if gene_name is missing
    my $gene_biotype = $attrs{gene_biotype} || 'NA';

    push @genes, {
        chr          => $chr,
        type         => $type,
        start        => $start,
        end          => $end,
        strand       => $strand,
        gene_id      => $gene_id,
        gene_name    => $gene_name,
        gene_biotype => $gene_biotype
    };
}
close $gtf_fh;
print "Loaded ", scalar(@genes), " genes from the GTF file.\n";

# Step 2: Process gene sequences and search for target pattern matches
print "\nStep 2: Processing gene sequences and identifying matches...\n";

# Open output file for writing results
open my $out_fh, '>', $output_file or die "Cannot open $output_file: $!";
print $out_fh join("\t", qw(
    Transcript_location
    Feature
    Strand
    GeneID
    GeneName
    GeneBiotype
    Target_location
    MatchedSequence
    MismatchCount
    Upstream20bp
)), "\n";

# Store merged results for handling overlapping matches
my %merged_data;

foreach my $gene (@genes) {
    # Extract gene sequence including flanking regions
    my $extract_start = $gene->{start} - $flanking_size;
    my $extract_end   = $gene->{end} + $flanking_size;

    my $region = "$gene->{chr}:$extract_start-$extract_end";
    my $cmd = "samtools faidx $ref_genome $region 2>/dev/null";
    my $sequence = `$cmd`;
    $sequence =~ s/^>.*\n//;  # Remove FASTA header
    $sequence =~ s/\n//g;     # Remove newlines

    # Reverse complement if the strand is negative
    if ($gene->{strand} eq '-') {
        $sequence = reverse_complement($sequence);
    }

    # Trim the flanking regions to isolate the main gene sequence
    my $trimmed_seq = substr($sequence, $flanking_size, length($sequence) - 2 * $flanking_size);

    # Iterate through the trimmed sequence to find matches
    for (my $i = 0; $i <= length($trimmed_seq) - length($target_pattern); $i++) {
        my $candidate = substr($trimmed_seq, $i, length($target_pattern));

        # Count mismatches and check if positions 7-8 match 'GA'
        my $mismatch_count = 0;
        my $ga_mismatch = 0;

        for (my $j = 0; $j < length($target_pattern); $j++) {
            my $target_base    = substr($target_pattern, $j, 1);
            my $candidate_base = substr($candidate, $j, 1);

            if ($target_base ne $candidate_base) {
                $mismatch_count++;
                # Check if mismatch occurs at positions 7-8 (index 6-7)
                if (($j == 6 && $target_base ne 'G') || 
                    ($j == 7 && $target_base ne 'A')) {
                    $ga_mismatch = 1;
                    last;
                }
            }
        }

        # Process valid matches (mismatch count <= 2 and positions 7-8 match 'GA')
        if ($mismatch_count <= $max_mismatches && !$ga_mismatch) {
            # Calculate genomic coordinates
            my ($genomic_start, $genomic_end);

            if ($gene->{strand} eq '+') {
                $genomic_start = $gene->{start} + $i;
                $genomic_end   = $genomic_start + length($target_pattern) - 1;
            } else {
                my $trimmed_length = length($trimmed_seq);
                $genomic_end   = $gene->{end} - $i;
                $genomic_start = $genomic_end - length($target_pattern) + 1;
            }

            # Extract the upstream 20bp
            my $upstream_start = $flanking_size + $i - $flanking_size;
            $upstream_start = 0 if $upstream_start < 0;
            my $upstream20bp = substr($sequence, $upstream_start, $flanking_size);

            # Create a transcript location identifier
            my $transcript_location = "$gene->{chr}:$gene->{start}-$gene->{end}";
            my $tar_location = "$gene->{chr}:$genomic_start-$genomic_end";

            # Merge results for overlapping matches
            if (exists $merged_data{$tar_location}) {
                # Append unique attributes to avoid duplication
                unless (index($merged_data{$tar_location}{Feature}, $gene->{type}) != -1) {
                    $merged_data{$tar_location}{Feature} .= ";$gene->{type}";
                }
                unless (index($merged_data{$tar_location}{Transcript_location}, $transcript_location) != -1) {
                    $merged_data{$tar_location}{Transcript_location} .= ";$transcript_location";
                }
                unless (index($merged_data{$tar_location}{GeneID}, $gene->{gene_id}) != -1) {
                    $merged_data{$tar_location}{GeneID} .= ";$gene->{gene_id}";
                }
                unless (index($merged_data{$tar_location}{GeneName}, $gene->{gene_name}) != -1) {
                    $merged_data{$tar_location}{GeneName} .= ";$gene->{gene_name}";
                }
                unless (index($merged_data{$tar_location}{GeneBiotype}, $gene->{gene_biotype}) != -1) {
                    $merged_data{$tar_location}{GeneBiotype} .= ";$gene->{gene_biotype}";
                }
            } else {
                # Create a new entry for the match
                $merged_data{$tar_location} = {
                    Transcript_location => $transcript_location,
                    Feature             => $gene->{type},
                    Strand              => $gene->{strand},
                    GeneID              => $gene->{gene_id},
                    GeneName            => $gene->{gene_name},
                    GeneBiotype         => $gene->{gene_biotype},
                    Target_location     => $tar_location,
                    MatchedSequence     => $candidate,
                    MismatchCount       => $mismatch_count,
                    Upstream20bp        => $upstream20bp
                };
            }
        }
    }
}

# Write merged results to the output file
foreach my $tar_location (sort keys %merged_data) {
    my $data = $merged_data{$tar_location};
    print $out_fh join("\t",
        $data->{Transcript_location},
        #$data->{Feature},
        $data->{Strand},
        $data->{GeneID},
        $data->{GeneName},
        $data->{GeneBiotype},
        $data->{Target_location},
        $data->{MatchedSequence},
        $data->{MismatchCount},
        $data->{Upstream20bp}
    ), "\n";
}

close $out_fh;
print "Final results saved to $output_file\n";

# Subroutine: Reverse complement of a DNA sequence
sub reverse_complement {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

print "\nAll steps completed successfully!\n";
