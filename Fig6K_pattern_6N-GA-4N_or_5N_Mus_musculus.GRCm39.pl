#!/usr/bin/perl
use strict;
use warnings;

##########################################
# Configuration
##########################################
my $gtf_file     = "Mus_musculus.GRCm39.113.chr.gtf";                     # Input GTF file
my $ref_genome   = "Mus_musculus.GRCm39.dna.primary_assembly.fa";         # Reference genome file
my $output_file1 = "Results_Fig6K_pattern_6N-GA-4N_Mus_musculus.GRCm39.txt"; # Output file for 6N-GA-4N
my $output_file2 = "Results_Fig6K_pattern_6N-GA-5N_Mus_musculus.GRCm39.txt"; # Output file for 6N-GA-5N

##########################################
# Initialize
##########################################
open my $gtf_fh, '<', $gtf_file or die "Cannot open $gtf_file: $!";
open my $out_fh1, '>', $output_file1 or die "Cannot open $output_file1: $!";
open my $out_fh2, '>', $output_file2 or die "Cannot open $output_file2: $!";

# Hash to store results
my %ga_sites_6N4N;  # For pattern NNNNNNGANNNN
my %ga_sites_6N5N;  # For pattern NNNNNNGANNNNN

##########################################
# Process GTF file
##########################################
while (<$gtf_fh>) {
    chomp;
    next if /^#/;  # Skip comments

    # Parse GTF fields
    my ($chr, $source, $feature, $start, $end, $score, $strand, $frame, $attr_str) = split /\t/;
    next unless $feature eq 'gene';  # Process only "gene" features


    ##########################################
    # Extract sequence from reference genome
    ##########################################
    my $extract_start = $start;
    my $extract_end   = $end;
    my $region        = "$chr:$extract_start-$extract_end";

    # Retrieve the sequence for the region using samtools
    my $seq = `samtools faidx $ref_genome $region 2>/dev/null`;
    $seq =~ s/^>.*\n//; # Remove FASTA header
    $seq =~ s/\n//g;    # Remove line breaks

    # Reverse complement if the strand is negative
    $seq = reverse_complement($seq) if $strand eq '-';

    ##########################################
    # Search for patterns
    ##########################################
    # NNNNNNGANNNN
	while ($seq =~ /(?=(.{6})GA(.{4}))/g) {
		my $upstream  = $1;
		my $downstream = $2;
		my $context = $upstream . 'GA' . $downstream;  # Full matched context
		
		if($strand eq '+'){
			# Calculate the position of the match based on the current index
			my $match_start = pos($seq)+$extract_start; # Adjust the position
			my $match_end=pos($seq)+$extract_start+11;
			
			if(exists $ga_sites_6N4N{$context}){
				if($ga_sites_6N4N{$context}[0]!~/$chr:$match_start-$match_end/){
					$ga_sites_6N4N{$context}[0]=$ga_sites_6N4N{$context}[0].";"."$chr:$match_start-$match_end";
					$ga_sites_6N4N{$context}[1]+=1;
				}
			}else{
				$ga_sites_6N4N{$context}[0]="$chr:$match_start-$match_end";
				$ga_sites_6N4N{$context}[1]=1;
			}
		}elsif($strand eq "-"){
			
			my $match_start = $extract_end-pos($seq)-length($context)+1; # Adjust the position
			my $match_end=$extract_end-pos($seq);
			
			if(exists $ga_sites_6N4N{$context}){
				if($ga_sites_6N4N{$context}[0]!~/$chr:$match_start-$match_end/){
					$ga_sites_6N4N{$context}[0]=$ga_sites_6N4N{$context}[0].";"."$chr:$match_start-$match_end";
					$ga_sites_6N4N{$context}[1]+=1;
				}
			}else{
				$ga_sites_6N4N{$context}[0]="$chr:$match_start-$match_end";
				$ga_sites_6N4N{$context}[1]=1;
			}

		}
	}

	while ($seq =~ /(?=(.{6})GA(.{5}))/g) {
		my $upstream  = $1;
		my $downstream = $2;
		my $context = $upstream . 'GA' . $downstream;  # Full matched context
		
		if($strand eq '+'){
			# Calculate the position of the match based on the current index
			my $match_start = pos($seq)+$extract_start; # Adjust the position
			my $match_end=pos($seq)+$extract_start+12;
			
			if(exists $ga_sites_6N5N{$context}){
				if($ga_sites_6N5N{$context}[0]!~/$chr:$match_start-$match_end/){
					$ga_sites_6N5N{$context}[0]=$ga_sites_6N5N{$context}[0].";"."$chr:$match_start-$match_end";
					$ga_sites_6N5N{$context}[1]+=1;
				}
			}else{
				$ga_sites_6N5N{$context}[0]="$chr:$match_start-$match_end";
				$ga_sites_6N5N{$context}[1]=1;
			}
		}elsif($strand eq "-"){
			
			my $match_start = $extract_end-pos($seq)-length($context)+1; # Adjust the position
			my $match_end=$extract_end-pos($seq);
			
			if(exists $ga_sites_6N5N{$context}){
				if($ga_sites_6N5N{$context}[0]!~/$chr:$match_start-$match_end/){
					$ga_sites_6N5N{$context}[0]=$ga_sites_6N5N{$context}[0].";"."$chr:$match_start-$match_end";
					$ga_sites_6N5N{$context}[1]+=1;
				}
			}else{
				$ga_sites_6N5N{$context}[0]="$chr:$match_start-$match_end";
				$ga_sites_6N5N{$context}[1]=1;
			}

		}
	}


}

##########################################
# Output Results
##########################################

# Add headers for output files
print $out_fh1 "Sequence\tFrequency\n";
print $out_fh2 "Sequence\tFrequency\n";

# Output for NNNNNNGANNNN
foreach my $site (sort keys %ga_sites_6N4N) {
    #print $out_fh1 "$site\t$ga_sites_6N4N{$site}[0]\t$ga_sites_6N4N{$site}[1]\n";
	 print $out_fh1 "$site\t$ga_sites_6N4N{$site}[1]\n";
}

# Output for NNNNNNGANNNNN
foreach my $site (sort keys %ga_sites_6N5N) {
    #print $out_fh2 "$site\t$ga_sites_6N5N{$site}[0]\t$ga_sites_6N5N{$site}[1]\n";
	print $out_fh2 "$site\t$ga_sites_6N5N{$site}[1]\n";
}

# Close file handles
close $gtf_fh;
close $out_fh1;
close $out_fh2;

print "Pipeline completed successfully! Results saved to $output_file1 and $output_file2\n";


# Generate reverse complement of a DNA sequence
sub reverse_complement {
    my $seq = shift;
    $seq = reverse $seq;                  # Reverse the sequence
    $seq =~ tr/ACGTacgt/TGCAtgca/;        # Complement the bases
    return $seq;
}

