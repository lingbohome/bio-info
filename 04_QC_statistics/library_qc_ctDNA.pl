#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options(
'
Description

	This program is used for library QC process for ctDNA with input list file of bam and fastqc data.

Version

	Name:             library_qc_ctDNA.pl

Requirements

	SAMtools-v1.8
	bedtools-v2.25.0

Options

	library_qc_ctDNA.pl %o <some-arg>
',
	[ 'list|l=s',         "list file: sample_name, sort_bam (contains duplications), rmdup_bam and fastqc_data.txt of read1 and read2 created by FastQC-v0.11.5; seperated by \"\\t\"; multi-samples are avalible (required)",                                                                             ],
	[ 'bed|b=s',          "target region bed file (required)",                                                                        ],
	[ 'prefix|p=s',       "output prefix (required)",                                                                                 ],
	[ 'outdir|o=s',       "output directory (default: ./)",                                                       { default => "./" } ],
	[ 'samtools_dir|S=s', "directory of samtools-v1.8 (default: /data/workdir/software/samtools-1.8)",     { default => "/data/workdir/software/samtools-1.8" }                                                                                                                              ],
	[ 'bedtools_dir|B=s', "directory of bedtools-v2.25.0 (default: /data/workdir/software/bedtools2/bin)", { default => "/data/workdir/software/bedtools2/bin" }                                                                                                                             ],
	[ 'help|h',           "print usage message and exit",                                                                             ],
);
my ($list, $bed, $prefix, $outdir, $samtools_dir, $bedtools_dir, $help) = ($opt->list, $opt->bed, $opt->prefix, $opt->outdir, $opt->samtools_dir, $opt->bedtools_dir, $opt->help);
print $usage->text . "\n" and exit if $help;

print $usage->text . "\n" and exit unless ($list && $bed && $prefix && $outdir && $samtools_dir && $bedtools_dir);

###################################### software ######################################
# Add path of samtools-v1.8
$ENV{PATH} = $samtools_dir . ":$ENV{PATH}";
# Add path of bedtools-v2.25.0
$ENV{PATH} = $bedtools_dir . ":$ENV{PATH}";

#################################### main program ####################################
my $panel_size;
&panel_count;
open LIST, '<', $list or die $!;
open OUT, '>', "$outdir/$prefix.library_QC.tmp" or die $!;
print OUT "#Sample Name\tQ20 of Read1s\tQ20 of Read2s\tQ20 of Paired Reads\tQ30 of Read1s\tQ30 of Read2s\tQ30 of Paired Reads\tTotal Reads\tTotal Duplicate Reads\tDuplicate Rate\tFraction of Reads Mapped to Genome\tTotal Bases Mapped to Genome\tTotal Bases Mapped to Target Region\tFraction of Bases Mapped to Target Region\tTotal Bases Mapped to Target Region after Duplicates Removed\tFraction of Bases Mapped to Target Region after Duplicates Removed\t1X Base Coverage of Target Region\t1X Base Coverage of Target Region after Duplicates Removed\t20X Base Coverage of Target Region\t20X Base Coverage of Target Region after Duplicates Removed\t100X Base Coverage of Target Region\t100X Base Coverage of Target Region after Duplicates Removed\t500X Base Coverage of Target Region\t500X Base Coverage of Target Region after Duplicates Removed\t1000X Base Coverage of Target Region\t1000X Base Coverage of Target Region after Duplicates Removed\t2000X Base Coverage of Target Region\t2000X Base Coverage of Target Region after Duplicates Removed\t4000X Base Coverage of Target Region\t4000X Base Coverage of Target Region after Duplicates Removed\tAverage Depth of Target Region\tAverage Depth of Target Region after Duplicates Removed\tRemark\n";
while (<LIST>){
	chomp;
	my @cut = split /\t/, $_;
	my ($sample, $bam, $rmdup_bam, $r1_fastqc, $r2_fastqc) = ($cut[0], $cut[1], $cut[2], $cut[3], $cut[4]);
	my $out = &lib_qc($sample, $bam, $rmdup_bam, $r1_fastqc, $r2_fastqc);
	print OUT $out;
}
close LIST;
close OUT;
`mv $outdir/$prefix.library_QC.tmp $outdir/$prefix.library_QC.xls`;

##################################### sub model ######################################
sub panel_count{
	open PANEL, '<', $bed or die $!;
	while (<PANEL>){
		chomp;
		my @cut = split /\t/, $_;
		$panel_size += $cut[2] - $cut[1];
	}
	close PANEL;
}

sub lib_qc{
	my ($id, $b, $rb, $r1, $r2) = ($_[0], $_[1], $_[2], $_[3], $_[4]);

	##### calculate map rate(before rmdup), total reads(before & after rmdup), total duplication reads and duplication rate #####
	my (@flagstat_b, @flagstat_rb, $total_reads, $total_reads_rmdup, $total_dup_reads, $dup_rate, $map_reads_rate);
	@flagstat_b = `samtools flagstat $b`;
	@flagstat_rb = `samtools flagstat $rb`;
	for (@flagstat_b){
		if (/(\d+) \+ \d+ in total/){
			$total_reads = $1;
			next;
		}elsif (/mapped \((.+%) *:.+\)/){
			$map_reads_rate = $1;
			next;
		}
	}
	for (@flagstat_rb){
		if (/(\d+) \+ \d+ in total/){
			$total_reads_rmdup = $1;
			next;
		}
	}
	$total_dup_reads = $total_reads - $total_reads_rmdup;
	$dup_rate = (sprintf "%0.2f", ($total_dup_reads / $total_reads * 100))."\%";

	##### calculate total base mapped to genome(before rmdup), total base mapped to target region(before & after rmdup), fraction of base mapped to target region(before & after rmdup), coverage of target region(before & after rmdup), average depth of target region(before & after rmdup) #####
	my ($total_base_mapgenome, $total_base_targeted, $total_rmdup_base_targeted, $base_targeted_fraction, $rmdup_base_targeted_fraction, $remark);
	my ($coverage_1x, $rmdup_coverage_1x, $coverage_20x, $rmdup_coverage_20x, $coverage_100x, $rmdup_coverage_100x, $coverage_500x, $rmdup_coverage_500x, $coverage_1000x, $rmdup_coverage_1000x, $coverage_2000x, $rmdup_coverage_2000x, $coverage_4000x, $rmdup_coverage_4000x) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	`bedtools bamtobed -i $b >$outdir/$id.bamtobed.tmp`;
	open BAMTOBED, '<', "$outdir/$id.bamtobed.tmp" or die $!;
	while (<BAMTOBED>){
		my @cut = split /\t/, $_;
		$total_base_mapgenome += $cut[2] - $cut[1];
	}
	close BAMTOBED;
	`rm $outdir/$id.bamtobed.tmp`;
	`samtools depth -d 0 -b $bed $b >$outdir/$id.coverage.tmp`;
	open COVERAGE1, '<', "$outdir/$id.coverage.tmp" or die $!;
	while (<COVERAGE1>){
		chomp;
		my @cut = split /\t/, $_;
		$total_base_targeted += $cut[-1];
		$coverage_1x ++ if $cut[-1] >= 1;
		$coverage_20x ++ if $cut[-1] >= 20;
		$coverage_100x ++ if $cut[-1] >= 100;
		$coverage_500x ++ if $cut[-1] >= 500;
		$coverage_1000x ++ if $cut[-1] >= 1000;
		$coverage_2000x ++ if $cut[-1] >= 2000;
		$coverage_4000x ++ if $cut[-1] >= 4000;
	}
	close COVERAGE1;
	`rm $outdir/$id.coverage.tmp`;
	`samtools depth -d 0 -b $bed $rb >$outdir/$id.rmdup_coverage.tmp`;
	open COVERAGE2, '<', "$outdir/$id.rmdup_coverage.tmp" or die $!;
	while (<COVERAGE2>){
		chomp;
		my @cut = split /\t/, $_;
		$total_rmdup_base_targeted += $cut[-1];
		$rmdup_coverage_1x ++ if $cut[-1] >= 1;
		$rmdup_coverage_20x ++ if $cut[-1] >= 20;
		$rmdup_coverage_100x ++ if $cut[-1] >= 100;
		$rmdup_coverage_500x ++ if $cut[-1] >= 500;
		$rmdup_coverage_1000x ++ if $cut[-1] >= 1000;
		$rmdup_coverage_2000x ++ if $cut[-1] >= 2000;
		$rmdup_coverage_4000x ++ if $cut[-1] >= 4000;
	}
	close COVERAGE2;
	`rm $outdir/$id.rmdup_coverage.tmp`;
	$base_targeted_fraction = (sprintf "%0.2f", ($total_base_targeted / $total_base_mapgenome * 100))."\%";
	$rmdup_base_targeted_fraction = (sprintf "%0.2f", ($total_rmdup_base_targeted / $total_base_mapgenome * 100))."\%";
	my $fraction_1x = (sprintf "%0.2f", ($coverage_1x / $panel_size * 100))."\%";
	my $fraction_20x = (sprintf "%0.2f", ($coverage_20x / $panel_size * 100))."\%";
	my $fraction_100x = (sprintf "%0.2f", ($coverage_100x / $panel_size * 100))."\%";
	my $fraction_500x = (sprintf "%0.2f", ($coverage_500x / $panel_size * 100))."\%";
	my $fraction_1000x = (sprintf "%0.2f", ($coverage_1000x / $panel_size * 100))."\%";
	my $fraction_2000x = (sprintf "%0.2f", ($coverage_2000x / $panel_size * 100))."\%";
	my $fraction_4000x = (sprintf "%0.2f", ($coverage_4000x / $panel_size * 100))."\%";
	my $rmdup_fraction_1x = (sprintf "%0.2f", ($rmdup_coverage_1x / $panel_size * 100))."\%";
	my $rmdup_fraction_20x = (sprintf "%0.2f", ($rmdup_coverage_20x / $panel_size * 100))."\%";
	my $rmdup_fraction_100x = (sprintf "%0.2f", ($rmdup_coverage_100x / $panel_size * 100))."\%";
	my $rmdup_fraction_500x = (sprintf "%0.2f", ($rmdup_coverage_500x / $panel_size * 100))."\%";
	my $rmdup_fraction_1000x = (sprintf "%0.2f", ($rmdup_coverage_1000x / $panel_size * 100))."\%";
	my $rmdup_fraction_2000x = (sprintf "%0.2f", ($rmdup_coverage_2000x / $panel_size * 100))."\%";
	my $rmdup_fraction_4000x = (sprintf "%0.2f", ($rmdup_coverage_4000x / $panel_size * 100))."\%";
	my $average_depth = sprintf "%0.2f", ($total_base_targeted / $panel_size);
	my $average_rmdup_depth = sprintf "%0.2f", ($total_rmdup_base_targeted / $panel_size);
	if (($rmdup_coverage_1000x / $panel_size) >= 0.5) {
		$remark = "PASS";
	} elsif (($rmdup_coverage_500x / $panel_size) >= 0.5) {
		$remark = "PASS/REJECT";
	} else {
		$remark = "REJECT";
	}

	##### calculate Q20 and Q30 of rawdata #####
	my ($r1_q20_count, $r1_q30_count, $r1_total_count) = &q20_q30($r1);
	my ($r2_q20_count, $r2_q30_count, $r2_total_count) = &q20_q30($r2);
	my $r1_q20 = sprintf "%.2f%%", (($r1_q20_count / $r1_total_count) * 100);
	my $r1_q30 = sprintf "%.2f%%", (($r1_q30_count / $r1_total_count) * 100);
	my $r2_q20 = sprintf "%.2f%%", (($r2_q20_count / $r2_total_count) * 100);
	my $r2_q30 = sprintf "%.2f%%", (($r2_q30_count / $r2_total_count) * 100);
	my $paired_q20 = sprintf "%.2f%%", (($r1_q20_count + $r2_q20_count) / ($r1_total_count + $r2_total_count) * 100);
	my $paired_q30 = sprintf "%.2f%%", (($r1_q30_count + $r2_q30_count) / ($r1_total_count + $r2_total_count) * 100);

	##### report #####
	my $report;
	$report .= $id."\t";
	$report .= $r1_q20."\t";
	$report .= $r2_q20."\t";
	$report .= $paired_q20."\t";
	$report .= $r1_q30."\t";
	$report .= $r2_q30."\t";
	$report .= $paired_q30."\t";
	$report .= $total_reads."\t";
	$report .= $total_dup_reads."\t";
	$report .= $dup_rate."\t";
	$report .= $map_reads_rate."\t";
	$report .= $total_base_mapgenome."\t";
	$report .= $total_base_targeted."\t";
	$report .= $base_targeted_fraction."\t";
	$report .= $total_rmdup_base_targeted."\t";
	$report .= $rmdup_base_targeted_fraction."\t";
	$report .= $fraction_1x."\t";
	$report .= $rmdup_fraction_1x."\t";
	$report .= $fraction_20x."\t";
	$report .= $rmdup_fraction_20x."\t";
	$report .= $fraction_100x."\t";
	$report .= $rmdup_fraction_100x."\t";
	$report .= $fraction_500x."\t";
	$report .= $rmdup_fraction_500x."\t";
	$report .= $fraction_1000x."\t";
	$report .= $rmdup_fraction_1000x."\t";
	$report .= $fraction_2000x."\t";
	$report .= $rmdup_fraction_2000x."\t";
	$report .= $fraction_4000x."\t";
	$report .= $rmdup_fraction_4000x."\t";
	$report .= $average_depth."\t";
	$report .= $average_rmdup_depth."\t";
	$report .= $remark."\n";
	return $report;
}

sub q20_q30{
	my $fastqc = shift;
	my $flag = 0;
	my ($total_count, $q20_count, $q30_count);
	open FASTQC, $fastqc or die $!;
	while (<FASTQC>){
		chomp;
		next if /^#/;
		if (/^>>Per sequence quality scores/){
			$flag = 1;
			next;
		}
		if (/^>>END_MODULE/){
			$flag = 0;
		}
		if ($flag == 1){
			my ($quality, $count) = split /\t/, $_;
			$total_count += $count;
			if ($quality >= 20){
				$q20_count += $count;
			}
			if ($quality >= 30){
				$q30_count += $count;
			}
		}
	}
	return ($q20_count, $q30_count, $total_count);
}

