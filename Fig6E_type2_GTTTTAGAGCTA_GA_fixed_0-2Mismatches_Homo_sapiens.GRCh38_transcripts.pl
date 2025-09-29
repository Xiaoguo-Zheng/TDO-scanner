#!/usr/bin/perl
use strict;
use warnings;

##########################################
# Configuration
##########################################
my $gtf_file     = "Homo_sapiens.GRCh38.113.chr.gtf"; # Input GTF file
my $ref_genome   = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"; # Reference genome file
my $output_file    = 'Results_Fig6E_type2_GTTTTAGAGCTA_GA_fixed_02matches_Homo_sapiens.GRCh38_transcripts.txt'; # Output file
my $target_pattern = 'GTTTTAGAGCTA';

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
    # Search for GTTTTAGAGCTA mismatch<=2 pattern
    ##########################################
	 for (my $i = 0; $i <= length($seq) - length($target_pattern); $i++) {
        my $candidate = substr($seq, $i, length($target_pattern));

        # Count mismatches and check if positions 7-8 match 'GA'
		 my $upstream_start;
		 my $upstream_end;
        if(substr($candidate,6,2) eq "GA"){
			my $mismatch=($candidate^$target_pattern)=~tr/\0//c;

			# Process valid matches (mismatch count <= 2 and positions 7-8 match 'GA')
			if ($mismatch<= 2) {
				# Calculate genomic coordinates
				my ($genomic_start, $genomic_end);

				if ($strand eq '+') {
					$genomic_start = $start + $i;
					$genomic_end   = $genomic_start + length($target_pattern) - 1;
					$upstream_start=$genomic_start-20;
					$upstream_end=$genomic_start-1;
				} else {			
					$genomic_end   = $end - $i;
					$genomic_start = $genomic_end - length($target_pattern) + 1;
					$upstream_start=$genomic_end+1;
					$upstream_end=$genomic_end+20;
				}

				# Extract the upstream 20bp
				my $upstream_location="$chr:$upstream_start-$upstream_end";
				my $upstream20bp = `samtools faidx $ref_genome $upstream_location 2>/dev/null`;
				$upstream20bp =~ s/^>.*\n//; # Remove FASTA header
				$upstream20bp=~ s/\n//g;    # Remove line breaks

				# Use MatchedSequence + Location as a unique key
				my $location = "$chr:$genomic_start-$genomic_end";
				my $key = "$candidate|$location";

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
					MatchedSeq          => $candidate,
					Location            => $location,
					Mismatch 			=> $mismatch,
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
	 }
}
close $gtf_fh;

#########################################
# Write results to output file
##########################################
open my $out_fh, '>', $output_file or die "Cannot write to $output_file: $!";

# Print header
print $out_fh join("\t", qw(
    Transcript_location Strand GeneID TranscriptID GeneName TranscriptName GeneBiotype TranscriptBiotype MatchedSequence Location Mismatch Upstream20bp
)), "\n";

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
        $r->{Mismatch},
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

	
	
