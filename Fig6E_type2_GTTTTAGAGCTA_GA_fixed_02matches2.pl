#!/usr/bin/perl
use strict;
use warnings;

# 配置参数
my $gtf_file = '../Homo_sapiens.GRCh38.113.chr.gtf';
my $ref_genome = '../Homo_sapiens.GRCh38.dna.primary_assembly.fa';
my $target_pattern = 'GTTTTAGAGCTA';
my $flanking_size = 20;
my $max_mismatches = 2;
my $output_file = 'Fig6E_type2_GTTTTAGAGCTA_GA_fixed_02matches.txt';

# 检查samtools是否可用
my $samtools_path = `which samtools` or die "samtools not found in PATH\n";
chomp $samtools_path;

# 步骤1: 从GTF文件中提取基因信息并直接处理
print "Step 1: Extracting and processing gene information from GTF file...\n";

# 存储基因信息
my @genes;

open my $gtf_fh, '<', $gtf_file or die "Cannot open $gtf_file: $!";
while (<$gtf_fh>) {
    next if /^#/;
    chomp;
    my @fields = split /\t/;
    next unless $fields[2] eq 'gene';
    
    my ($chr, $type, $start, $end, $strand) = @fields[0,2,3,4,6];
    
    my $attributes = $fields[8];
    my %attrs;
    while ($attributes =~ /(\w+)\s+"([^"]+)";/g) {
        $attrs{$1} = $2;
    }
    
    my $gene_id = $attrs{gene_id} || 'NA';
    my $gene_name = $attrs{gene_name} || $gene_id;
    my $gene_biotype = $attrs{gene_biotype} || 'NA';
    
    push @genes, {
        chr => $chr,
        type => $type,
        start => $start,
        end => $end,
        strand => $strand,
        gene_id => $gene_id,
        gene_name => $gene_name,
        gene_biotype => $gene_biotype
    };
}
close $gtf_fh;
print "Loaded ", scalar(@genes), " genes from GTF file\n";

# 步骤2: 直接处理基因序列并寻找匹配
print "\nStep 2: Processing gene sequences and searching for target pattern matches...\n";

# 打开输出文件
open my $out_fh, '>', $output_file or die "Cannot open $output_file: $!";
print $out_fh join("\t", qw(
    transcript_location
    Feature
    Strand
    GeneID
    GeneName
    GeneBiotype
    Tar_location
    MatchedSequence
    MismatchCount
    Upstream20bp
)), "\n";

# 用于存储合并的数据
my %merged_data;

foreach my $gene (@genes) {
    # 提取基因序列及侧翼区域
    my $extract_start = $gene->{start} - $flanking_size;
    my $extract_end = $gene->{end} + $flanking_size;
    
    my $region = "$gene->{chr}:$extract_start-$extract_end";
    my $cmd = "samtools faidx $ref_genome $region 2>/dev/null";
    my $sequence = `$cmd`;
    $sequence =~ s/^>.*\n//;
    $sequence =~ s/\n//g;
    
    if ($gene->{strand} eq '-') {
        $sequence = reverse_complement($sequence);
    }
    
    # 去掉头部和尾部各20bp
    my $trimmed_seq = substr($sequence, $flanking_size, length($sequence)-2*$flanking_size);

    # 遍历trimmed_seq，寻找匹配
    for (my $i = 0; $i <= length($trimmed_seq) - length($target_pattern); $i++) {
        my $candidate = substr($trimmed_seq, $i, length($target_pattern));

        # 计算错配数量，并检查第7-8位是否匹配GA
        my $mismatch_count = 0;
        my $ga_mismatch = 0;

        for (my $j = 0; $j < length($target_pattern); $j++) {
            my $target_base = substr($target_pattern, $j, 1);
            my $candidate_base = substr($candidate, $j, 1);
            
            if ($target_base ne $candidate_base) {
                $mismatch_count++;
                # 检查是否在第7-8位（索引6-7）
                if (($j == 6 && $target_base ne 'G') || 
                    ($j == 7 && $target_base ne 'A')) {
                    $ga_mismatch = 1;
                    last;
                }
            }
        }

        # 如果错配数量<=2且GA位置正确，则处理
        if ($mismatch_count <= $max_mismatches && !$ga_mismatch) {
            # 计算基因组位置
            my ($genomic_start, $genomic_end);
            
            if ($gene->{strand} eq '+') {
                $genomic_start = $gene->{start} + $i;
                $genomic_end = $genomic_start + length($target_pattern) - 1;
            } else {
                my $trimmed_length = length($trimmed_seq);
                $genomic_end = $gene->{end} - $i;
                $genomic_start = $genomic_end - length($target_pattern) + 1;
            }

            # 提取上游20bp
            my $upstream_start = $flanking_size + $i - $flanking_size;
            $upstream_start = 0 if $upstream_start < 0;
            my $upstream20bp = substr($sequence, $upstream_start, $flanking_size);

            # 创建位置标识符
            my $transcript_location = "$gene->{chr}:$gene->{start}-$gene->{end}";
            my $tar_location = "$gene->{chr}:$genomic_start-$genomic_end";
            
            # 合并数据
            if (exists $merged_data{$tar_location}) {
                # 合并Feature
                unless (index($merged_data{$tar_location}{Feature}, $gene->{type}) != -1) {
                    $merged_data{$tar_location}{Feature} .= ";$gene->{type}";
                }
                
                # 合并transcript_location
                unless (index($merged_data{$tar_location}{transcript_location}, $transcript_location) != -1) {
                    $merged_data{$tar_location}{transcript_location} .= ";$transcript_location";
                }
                
                # 合并GeneID
                unless (index($merged_data{$tar_location}{GeneID}, $gene->{gene_id}) != -1) {
                    $merged_data{$tar_location}{GeneID} .= ";$gene->{gene_id}";
                }
                
                # 合并GeneName
                unless (index($merged_data{$tar_location}{GeneName}, $gene->{gene_name}) != -1) {
                    $merged_data{$tar_location}{GeneName} .= ";$gene->{gene_name}";
                }
                
                # 合并GeneBiotype
                unless (index($merged_data{$tar_location}{GeneBiotype}, $gene->{gene_biotype}) != -1) {
                    $merged_data{$tar_location}{GeneBiotype} .= ";$gene->{gene_biotype}";
                }
            } else {
                # 创建新条目
                $merged_data{$tar_location} = {
                    transcript_location => $transcript_location,
                    Feature => $gene->{type},
                    Strand => $gene->{strand},
                    GeneID => $gene->{gene_id},
                    GeneName => $gene->{gene_name},
                    GeneBiotype => $gene->{gene_biotype},
                    Tar_location => $tar_location,
                    MatchedSequence => $candidate,
                    MismatchCount => $mismatch_count,
                    Upstream20bp => $upstream20bp
                };
            }
        }
    }
}

# 输出合并后的数据
foreach my $tar_location (sort keys %merged_data) {
    my $data = $merged_data{$tar_location};
    print $out_fh join("\t",
        $data->{transcript_location},
        $data->{Feature},
        $data->{Strand},
        $data->{GeneID},
        $data->{GeneName},
        $data->{GeneBiotype},
        $data->{Tar_location},
        $data->{MatchedSequence},
        $data->{MismatchCount},
        $data->{Upstream20bp}
    ), "\n";
}

close $out_fh;
print "Final results saved to $output_file\n";

# 反向互补子程序
sub reverse_complement {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

print "\nAll steps completed successfully!\n";
