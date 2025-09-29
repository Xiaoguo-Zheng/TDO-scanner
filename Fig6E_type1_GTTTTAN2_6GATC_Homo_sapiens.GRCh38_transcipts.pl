#!/usr/bin/perl
use strict;
use warnings;

##########################################
# Configuration
##########################################
my $gtf_file     = "Homo_sapiens.GRCh38.113.chr.gtf"; # Input GTF file
my $ref_genome   = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"; # Reference genome file
my $output_file  = "Results_Fig6E_type1_GTTTTAN2_6GATC_Homo_sapiens.GRCh38_transcripts.txt"; # Output file

##########################################
# Processing
##########################################
open my $gtf_fh, '<', $gtf_file or die "Cannot open $gtf_file: $!";

# Hash to store results; key = MatchedSequence + Location
my %results;

while (<$gtf_fh>) {
    chomp;
    next if /^#/;

    # Parse GTF fields
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attr_str) = split /\t/;
    next unless $feature eq 'transcript'; # Only process entries with the "transcript" feature

    # Parse attributes from the GTF record
    my %attr = parse_attributes($attr_str);
    my $gene_id          = $attr{gene_id} // 'NA';
    my $transcript_id    = $attr{transcript_id} // 'NA';
    my $gene_name        = $attr{gene_name} // $gene_id;
    my $transcript_name  = $attr{transcript_name} // $transcript_id;
    my $gene_biotype     = $attr{gene_biotype} // 'NA';
    my $transcript_biotype = $attr{transcript_biotype} // 'NA';

    ##########################################
    # Extract sequence with 20 bp flanks
    ##########################################
    my $region = "$chr:$start-$end";

    # Retrieve the sequence for the region using samtools
    my $seq = `samtools faidx $ref_genome $region 2>/dev/null`;
    $seq =~ s/^>.*\n//; # Remove FASTA header
    $seq =~ s/\n//g;    # Remove line breaks

    # Reverse complement if the strand is negative
    $seq = reverse_complement($seq) if $strand eq '-';

    ##########################################
    # Search for GTTTTA(N2-6)GCTA pattern
    ##########################################
    while ($seq =~ /(GTTTTA([ACGT]{2,6})GCTA)/g) {
        my $matched_seq = $1;       # Full matched sequence
        my $variable_part = $2;     # Variable portion (N2-6)

        # Match start and end positions within the trimmed sequence (1-based)
        my $match_start = pos($seq) - length($matched_seq) + 1;
        my $match_end   = pos($seq);

        my $upstream_start;
        my $upstream_end;

        # Convert to genomic coordinates
        my ($genomic_start, $genomic_end);
        if ($strand eq '+') {
            $genomic_start = $start + ($match_start - 1);
            $genomic_end   = $start + ($match_end - 1);
            $upstream_start  = $genomic_start - 20;
            $upstream_end    = $genomic_start - 1;
        } else {
            $genomic_end   = $end - ($match_start - 1);
            $genomic_start = $end - ($match_end - 1);
            $upstream_start  = $genomic_end + 1;
            $upstream_end    = $genomic_end + 20;
        }

        # Retrieve upstream 20 bp sequence
        my $upstream20 = "$chr:$upstream_start-$upstream_end";
        my $upstream20bp = `samtools faidx $ref_genome $upstream20 2>/dev/null`;
        $upstream20bp =~ s/^>.*\n//; # Remove FASTA header
        $upstream20bp =~ s/\n//g;    # Remove line breaks

        # Use MatchedSequence + Location as a unique key
        my $location = "$chr:$genomic_start-$genomic_end";
        my $key = "$matched_seq|$location";

        # Initialize result entry if the key does not exist
        if (!exists $results{$key}) {
            $results{$key} = {
                Transcript_location => '',
                Strand              => '',
                TranscriptID        => '',
                GeneID              => {},
                GeneName            => {},
                TranscriptName      => '',
                GeneBiotype         => {},
                TranscriptBiotype   => '',
                MatchedSeq          => $matched_seq,
                Location            => $location,
                VariablePart        => $variable_part,
                Upstream20bp        => $upstream20bp,
            };
        }

        # Append information using semicolon to separate values
        $results{$key}->{Transcript_location} .= "$chr:$start-$end; ";
        $results{$key}->{Strand} .= "$strand; ";
        $results{$key}->{GeneID}{$gene_id} = 1;
        $results{$key}->{TranscriptID} .= "$transcript_id; ";
        $results{$key}->{GeneName}{$gene_name} = 1;
        $results{$key}->{TranscriptName} .= "$transcript_name; ";
        $results{$key}->{GeneBiotype}{$gene_biotype} = 1;
        $results{$key}->{TranscriptBiotype} .= "$transcript_biotype; ";
    }
}
close $gtf_fh;

##########################################
# Write results to output file
##########################################
open my $out_fh, '>', $output_file or die "Cannot write to $output_file: $!";

# Print header
print $out_fh join("\t", qw(
    Transcript_location Strand GeneID TranscriptID GeneName TranscriptName GeneBiotype TranscriptBiotype MatchedSequence Location VariablePart Upstream20bp
)), "\n";

# Write each result entry
# Write each result entry
foreach my $key (sort keys %results) {
    my $r = $results{$key};

    # Remove trailing semicolons and spaces
    foreach my $field (qw(Transcript_location Strand TranscriptID TranscriptName TranscriptBiotype)) {
        $r->{$field} =~ s/;\s+$//;  # 去掉尾部的分号和空格
    }

    # 将哈希转换为以分号分隔的字符串
    my $gene_ids = join('; ', keys %{$r->{GeneID}});
    my $gene_names = join('; ', keys %{$r->{GeneName}});
    my $gene_biotypes = join('; ', keys %{$r->{GeneBiotype}});

    print $out_fh join("\t",
        $r->{Transcript_location},
        $r->{Strand},
        $gene_ids,          # 输出 GeneID 字符串
        $r->{TranscriptID}, # 已经是字符串
        $gene_names,        # 输出 GeneName 字符串
        $r->{TranscriptName}, # 已经是字符串
        $gene_biotypes,     # 输出 GeneBiotype 字符串
        $r->{TranscriptBiotype}, # 已经是字符串
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
