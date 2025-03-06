#!/usr/bin/env perl

use List::Util qw(min max sum);
use strict;
use warnings;

my $bamlist = shift @ARGV;
my $maflist = shift @ARGV;
my $reference = shift @ARGV;

if (!$bamlist || !$maflist || !$reference) {
	die "USAGE: $0 bam_list maf_list reference.fasta\n";
}

foreach my $input ($bamlist, $maflist, $reference) {
	if (! -e $input) {
		die "File $input does not exist\n";
	}
}


my %seen_var;

# Iterate through MAFs

my @maflist = `cat $maflist`;
chomp @maflist;
foreach my $maf (@maflist) {
	my %idx = %{ get_headers($maf) };
	open MAF, "<$maf" or die "Can't open MAF $maf\n";
	while (<MAF>) {
		next if /Hugo_Symbol/;
		my $line = $_;
		chomp $line;
		my @col = split(/\t/, $line);
		my $chr = $col[$idx{"Chromosome"}];
		my $start = $col[$idx{"Start_Position"}];
		my $ref = $col[$idx{"Reference_Allele"}];
		my $alt = $col[$idx{"Tumor_Seq_Allele2"}];
		my $sample = $col[$idx{"Tumor_Sample_Barcode"}];
		my $type = $col[$idx{"Variant_Type"}];
		my $end = $col[$idx{"End_Position"}];
		my $maf_start = $start;
		my $maf_end = $end;
		if ($type eq 'INS') {
			$end = $start;
		} elsif ($type eq 'DEL') {
			$start = $start - 1;
			$end = $start;
		}
		next unless $type =~ /DEL|INS|DNP|SNP|TNP/;
		my $key = "$chr:$start:$end:$ref:$alt";
		push @{$seen_var{$key}{"samples"}}, $sample;
		$seen_var{$key}{"type"} = $type;
		$seen_var{$key}{"maf_start"} = $maf_start;
		$seen_var{$key}{"maf_end"} = $maf_end;
		$seen_var{$key}{"Hugo_Symbol"} = $col[$idx{"Hugo_Symbol"}];
		$seen_var{$key}{"Gene"} = $col[$idx{"Gene"}];
		$seen_var{$key}{"HGVSp_Short"} = $col[$idx{"HGVSp_Short"}];
	}
}

my @col_headers = ("Hugo_Symbol", "Gene", "HGVSp_Short", "Chromosome", "Start_Position", 
  "End_Position", "Variant_Type", "Reference_Allele", "Tumour_Seq_Allele2", "Tumor_Sample_Barcode",
  "Alt_Count", "Tot_Count", "Alt_Perc", "Var_in_MAF", "Status");

print join("\t", @col_headers) . "\n";

# Check each variant in each sample

my @bamlist = `cat $bamlist`;
chomp @bamlist;
my @samplelist;

foreach my $bam (@bamlist) {
#	my $sample = $1 if $bam =~ /(^[\/]+).sample.dupmarked.bam/;
	my $sample = `samtools samples $bam | cut -f 1`;
	chomp $sample;

	print STDERR "Sample $sample inferred from BAM $bam\n";
	if (!$sample) {
		print STDERR "Can't file sample name from BAM file name.\n";
	};
	push @samplelist, $sample;
}

# Check each variant with mpileup

foreach my $var (sort keys %seen_var) {
	my ($chr, $start, $end, $ref, $alt) = split(/:/, $var);
	my @result = `samtools mpileup -q 0 -Q 0 -d 1000 -b $bamlist -A -r $chr:$start-$end --ignore-overlaps --no-BAQ --no-output-ins --no-output-del --no-output-ends --reference $reference 2>/dev/null`;
	if (!@result) {
		print STDERR "WARNING: $chr\t$start\t$end\t$ref\t$alt\tNO_PILEUP\n";
		next;
	}
	chomp @result;
	my $vartype = $seen_var{"$chr:$start:$end:$ref:$alt"}{"type"};
	my @mnp;
	if ($vartype =~ /NP/ && $vartype ne 'SNP') {
		@mnp = split('', $alt); 
	}
	# Iterate through each position of each variant
	my $totsamples = @samplelist;
	my $k = 3;
	my $bamnum = -1;
	# Iterate through each sample
	foreach my $sample (@samplelist) {
		$bamnum++;
		my $start_col = $k;
		my $end_col = $start_col + 3;
		my @tot_counts;
		my @tot_cov;
		my @tot_perc;
		foreach (my $i = $start_col; $i < $end_col; $i += 3) {
			my $status = "";
			my $found = 0;
			foreach (my $j = 0; $j <= $#result; $j++) {
				my @res = split (/\t/, $result[$j]);
				# Uncomment below to see pileups
				#print STDERR "$result[$j]\n";
				#print STDERR "$res[$i+1]\n";
				my $tot = $res[$i];
				push @tot_cov, $tot;
				#next if $tot == 0;
				my $bam = $bamlist[$bamnum];
				# Make all letters upper case
				$res[$i+1] = uc($res[$i+1]);
				my $alt_allele = $alt;
				if (@mnp) {
					$alt_allele = $mnp[$j];
#					print STDERR "Looking for $alt_allele from $alt\n";
				}
				my $counts = count_alts($vartype, $ref, $alt_allele, $res[$i+1]);
				push @tot_counts, $counts;
				my $perc = 100*$counts/$tot;
				push @tot_perc, $perc;
				#next if $counts == 0;
				$found++ if $counts > 0;
			}
			if (@mnp && $found > 0 && $found < scalar(@mnp)) {
				$status = "PART_MATCH";
			}
			my $mean_counts = mean(@tot_counts);
			my $mean_cov = mean(@tot_cov);
			my $mean_perc = mean(@tot_perc);
			my $key = "$chr:$start:$end:$ref:$alt";
			my $orig_start = $seen_var{$key}{"maf_start"};
			my $orig_end = $seen_var{$key}{"maf_end"};
			my $samples_in_maf = $seen_var{$key}{"samples"};
			my $type = $seen_var{$key}{"type"};
			my $hugo = $seen_var{$key}{"Hugo_Symbol"};
			my $gene = $seen_var{$key}{"Gene"};
			my $prot = $seen_var{$key}{"HGVSp_Short"};
			my $seen = grep(/^$sample$/, @{$samples_in_maf}) ? "PRESENT_IN_MAF" : "ABSENT_FROM_MAF";
			if (!$status) {
				if ($seen eq "PRESENT_IN_MAF") {
					if ($mean_perc > 0) {
						$status = "TRUE_POSITIVE";
					} else {
						$status = "FALSE_POSITIVE";
					}
				} else {
					if ($mean_perc == 0) {
						$status = "TRUE_NEGATIVE";
					} else {
						$status = "FALSE_NEGATIVE";
					}
				}
			}
			# Change the conditional if you want to exclude TRUE_NEGATIVEs
			# if ($status ne "TRUE_NEGATIVE" && $status ne "TRUE_POSITIVE") {
			# Un/comment the 'if' conditional if you want/don't want to print ALL results
#			if ($status ne "TRUE_NEGATIVE") {
				# "Hugo_Symbol", "Gene", "HGVSp_Short", "Chromosome", "Start_Position"
				# "End_Position", "Variant_Type", "Reference_Allele", "Tumour_Seq_Allele2", "Tumor_Sample_Barcode");
				print "$hugo\t$gene\t$prot\t$chr\t$orig_start\t$orig_end\t$type\t$ref\t$alt\t$sample\t$mean_counts\t$mean_cov\t$mean_perc\t$seen\t$status\n";
#			}
		}
		$k += 3;
	}
}


########## Subroutines ##########

# Get the column indexes from the MAFs
sub get_headers {
	my $maf = shift @_;
	my %idx;
	my $header = `head -n 1 $maf`;
	chomp $header;
	my @header_cols = split(/\t/, $header);
	while (my($idx, $header) = each(@header_cols)){
		$idx{$header} = $idx;
	}
	return \%idx;
}


sub count_alts {
	my ($type, $ref, $alt, $res) = @_;
	my $counts = 0;
	my $size = $type eq 'DEL' ? length($ref) : length($alt);
	if ($type =~ /NP/) {
		$counts = () = ($res =~ /$alt/g);
	} elsif ($type eq 'DEL') {
		#$counts = () = ($res =~ /\*/g);
		$counts = () = ($res =~ /-$size/g);
	} elsif ($type eq 'INS') {
		$counts = () = ($res =~ /\+$size/g);
	}
	return $counts;
}

sub mean {
	my @num = @_;
	return sum(@_)/ @_;
}
