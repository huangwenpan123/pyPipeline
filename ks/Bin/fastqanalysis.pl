#!/usr/bin/perl

###########################################################
## Author: Kxyao
## Email : 70910595@qq.com
## Date  : `date`
## Ver   : 0.99.1
## Usage : This script is use for ...
###########################################################

#use utf8;
use strict;
use warnings;
use Carp qw/carp croak/;
use Data::Dumper;

#use Bio::Seq;
#use Bio::SeqIO;

#use Cwd qw/abs_path/;
use FindBin qw/$Bin $Script/;
use File::Basename qw/dirname basename/;
#use FileHandle;
#use Spreadsheet::XLSX;
#use List::Util qw/max min sum/;
#use List::MoreUtils qw/any all uniq/;

#use threads;
#use threads::shared;
#use Thread::Semaphore;
#use Thread::Queue;

my %opts;
use Getopt::Long;
GetOptions (\%opts,
	"fq1=s",
	"fq2=s",
	"primer=s",
	"olp=s",
	"l=i",
	"lr=i",
	"fa=s",
	"plasmidfa=s",
	"h37rvfa=s",
	"out=s",
	"JHann=s",
	"database=s",
	"samtools=s",
	"fastp=s",
	"bwa=s",
	"name=s",
	"SE=s",
	"model=s",
	"IAC=s",
	"help" => \&help,
);

###########################################################

#die &help unless ($opts{in} && $opts{out});
#die &help if ($opts{help} );
## Your code go here
my $name=$opts{name};
$opts{model}||="full";
system("mkdir -p $opts{out}/$opts{name}/align");
system("mkdir -p $opts{out}/$opts{name}/clean");
if (defined$opts{fq2} && -e $opts{fq2})
{
	system("$opts{fastp} -i $opts{fq1} -I $opts{fq2} -o $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz -O $opts{out}/$opts{name}/clean/${name}_clean_2.fq.gz -l 30 -c -q 20 -u 20 -n 6  --adapter_sequence TGGAATTCTCGGGTGCCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGT -j $opts{out}/$opts{name}/clean/${name}.json -h $opts{out}/$opts{name}/clean/${name}.html");

	if (defined$opts{primer} and -e "$opts{primer}")
	{
		&primer($opts{primer}, $opts{fq1}, $opts{fq2}, "$opts{out}/$opts{name}/clean/$opts{name}", $opts{l},$opts{olp},$opts{lr},"PE");
	}
	else
	{
		system("cp -f $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz $opts{out}/$opts{name}/clean/${name}_1.fq.gz");
		system("cp -f $opts{out}/$opts{name}/clean/${name}_clean_2.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz");
	}
	#######
	system("rm -rf $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz $opts{out}/$opts{name}/clean/${name}_clean_2.fq.gz\n");
	##########
	system("$opts{bwa} mem -t 16 $opts{fa} $opts{out}/$opts{name}/clean/${name}_1.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz > $opts{out}/$opts{name}/align/$name.sam");
	system("$opts{bwa} mem -t 16 $opts{plasmidfa} $opts{out}/$opts{name}/clean/${name}_1.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz > $opts{out}/$opts{name}/align/$name.plasmid.sam");
	system("$opts{samtools} view -F 8 $opts{out}/$opts{name}/align/$name.sam > $opts{out}/$opts{name}/align/$name.align.sam");
	system("$opts{samtools} view -F 8 $opts{out}/$opts{name}/align/$name.plasmid.sam > $opts{out}/$opts{name}/align/$name.align.plasmid.sam");
	&sam2result("$opts{out}/$opts{name}/align/$name.align.sam","$opts{out}/$opts{name}/align/$name.result.xls","$opts{out}/$opts{name}/align/$name.align.plasmid.sam",25,"$opts{fa}.ann","no",$opts{IAC});
	######
	if ($opts{model} eq "silent")
	{
		system("rm -rf $opts{out}/$opts{name}/align/$name.sam $opts{out}/$opts{name}/align/$name.plasmid.sam $opts{out}/$opts{name}/align/$name.align.sam $opts{out}/$opts{name}/align/$name.align.plasmid.sam");
	}
	######
	if ($opts{model} eq "silent")
	{
		system("rm -rf $opts{out}/$opts{name}/align/${name}.JHalign.bam $opts{out}/$opts{name}/align/${name}.JHsort.bam $opts{out}/$opts{name}/align/${name}.JH.filter.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam.bai $opts{out}/$opts{name}/align/${name}.JH.pileup $opts{out}/$opts{name}/align/${name}.JH.pileup.out $opts{out}/$opts{name}/align/${name}.JH.fileter.pileup2ann.txt $opts{out}/$opts{name}/align/${name}.JH.refGene.variant_function $opts{out}/$opts{name}/align/${name}.JH.refGene.exonic_variant_function $opts{out}/$opts{name}/align/${name}.JH.refGene.log $opts{out}/$opts{name}/align/${name}.JH.JH_multianno.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_ann.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_positive.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_all.xls $opts{out}/$opts{name}/clean/${name}_1.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz");
	}
	else
	{
		system("rm -rf $opts{out}/$opts{name}/align/${name}.JHalign.bam $opts{out}/$opts{name}/align/${name}.JHsort.bam $opts{out}/$opts{name}/align/${name}.JH.filter.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam.bai $opts{out}/$opts{name}/align/${name}.JH.pileup $opts{out}/$opts{name}/align/${name}.JH.pileup.out $opts{out}/$opts{name}/align/$name.sam $opts{out}/$opts{name}/align/$name.plasmid.sam $opts{out}/$opts{name}/align/$name.align.sam $opts{out}/$opts{name}/align/$name.align.plasmid.sam $opts{out}/$opts{name}/clean/${name}_1.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz");
	}
}
elsif (!defined$opts{fq2} || $opts{SE} =~/yes/i)
{
	system("$opts{fastp} -i $opts{fq1} -o $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz -l 30 -q 20 -u 20 -n 6  --adapter_sequence TGGAATTCTCGGGTGCCA -j $opts{out}/$opts{name}/clean/${name}.json -h $opts{out}/$opts{name}/clean/${name}.html");

	if (defined$opts{primer} && -e "$opts{primer}")
	{
		&primer($opts{primer}, $opts{fq1}, "", "$opts{out}/$opts{name}/clean/$opts{name}", $opts{l},$opts{olp},$opts{lr},"SE");
	}
	else
	{
		system("cp -f $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz $opts{out}/$opts{name}/clean/${name}_single_1.fq.gz");
	}
	#################
	system("rm -rf $opts{out}/$opts{name}/clean/${name}_clean_1.fq.gz");
	#################
	system("$opts{bwa} mem -t 16 $opts{fa} $opts{out}/$opts{name}/clean/${name}_single_1.fq.gz > $opts{out}/$opts{name}/align/$name.sam");
	system("$opts{bwa} mem -t 16 $opts{plasmidfa} $opts{out}/$opts{name}/clean/${name}_single_1.fq.gz > $opts{out}/$opts{name}/align/$name.plasmid.sam");
	system("$opts{samtools} view -F 4 $opts{out}/$opts{name}/align/$name.sam > $opts{out}/$opts{name}/align/$name.align.sam");
	system("$opts{samtools} view -F 4 $opts{out}/$opts{name}/align/$name.plasmid.sam > $opts{out}/$opts{name}/align/$name.align.plasmid.sam");
	&sam2result("$opts{out}/$opts{name}/align/$name.align.sam","$opts{out}/$opts{name}/align/$name.result.xls","$opts{out}/$opts{name}/align/$name.align.plasmid.sam",25,"$opts{fa}.ann","yes",$opts{IAC});
	########
	if ($opts{model} eq "silent")
	{
		system("rm -rf $opts{out}/$opts{name}/align/$name.sam $opts{out}/$opts{name}/align/$name.plasmid.sam $opts{out}/$opts{name}/align/$name.align.sam $opts{out}/$opts{name}/align/$name.align.plasmid.sam ");
	}
	########
	####################
	if ($opts{model} eq "silent")
	{
		system("rm -rf $opts{out}/$opts{name}/align/${name}.JHalign.bam $opts{out}/$opts{name}/align/${name}.JHsort.bam $opts{out}/$opts{name}/align/${name}.JH.filter.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam.bai $opts{out}/$opts{name}/align/${name}.JH.pileup $opts{out}/$opts{name}/align/${name}.JH.pileup.out $opts{out}/$opts{name}/align/${name}.JH.fileter.pileup2ann.txt $opts{out}/$opts{name}/align/${name}.JH.refGene.variant_function $opts{out}/$opts{name}/align/${name}.JH.refGene.exonic_variant_function $opts{out}/$opts{name}/align/${name}.JH.refGene.log  $opts{out}/$opts{name}/align/${name}.JH.JH_multianno.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_ann.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_positive.xls $opts{out}/$opts{name}/align/${name}.JH.hotspot_all.xls $opts{out}/$opts{name}/clean/${name}_single_1.fq.gz");
	}
	else
	{
		system("rm -rf $opts{out}/$opts{name}/align/${name}.JHalign.bam $opts{out}/$opts{name}/align/${name}.JHsort.bam $opts{out}/$opts{name}/align/${name}.JH.filter.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam $opts{out}/$opts{name}/align/${name}.JH.sort.bam.bai $opts{out}/$opts{name}/align/${name}.JH.pileup $opts{out}/$opts{name}/align/${name}.JH.pileup.out $opts{out}/$opts{name}/align/$name.sam $opts{out}/$opts{name}/align/$name.plasmid.sam $opts{out}/$opts{name}/align/$name.align.sam $opts{out}/$opts{name}/align/$name.align.plasmid.sam $opts{out}/$opts{name}/clean/${name}_1.fq.gz $opts{out}/$opts{name}/clean/${name}_2.fq.gz $opts{out}/$opts{name}/clean/${name}_single_1.fq.gz");
	}
}

######################### Subs ############################
sub help{
	print <<HELP;
Name
	$0
Description
	This script is used for ...
Options
	-in        <str>     ** input file
	-out       <str>     ** output file
	-help                   help information
Author
	kxyao	70910595\@qq.com
Version
	0.99.1
Usage
	perl $0 [options] -in <input> -out <output>
HELP
	exit 0;
}

sub primer{
	my %optsprimer;
	($optsprimer{p},$optsprimer{fq1},$optsprimer{fq2},$optsprimer{o},$optsprimer{l},$optsprimer{olp},$optsprimer{lr},$optsprimer{fqtype})=@_;

	$optsprimer{seed}=7;                 ######引物识别seed长度
	$optsprimer{l} ||= 50;               ######reads dimer 过滤长度
	$optsprimer{m} ||= 2 ;               ######引物识别mismtach
	$optsprimer{lower} ||= 2/3;          ######引物最短识别长度比值。
	$optsprimer{unmap_primer}||="no";    ######输出无法识别的引物
	$optsprimer{nonspe} ||="no";         ######输出引物非特异结合reads 的fastq
	$optsprimer{fqtype} ||= "PE";        ######针对韦翰斯150+40数据
	$optsprimer{lr}||=$optsprimer{l};
	$optsprimer{cut}||="no";

	my %primer_result;
	my %primer_expand;
	my %subprimer;

	my $basename=basename($optsprimer{o});
###########read primer
	my %primer;
	my %primer_re;
	my %fseed;my %rseed;my %fseed1;my %rseed1;my %hash_fseed;my %hash_rseed;

	my %primer_2_ref;
	my %error_primer;
	my %expression;
	my $total_primer;

	my %overlap;
	if (defined$optsprimer{olp})
	{
		open OVERLAP,"$Bin/configd.dist -z $optsprimer{olp}|" or die $!;
		while (<OVERLAP>)
		{
			chomp;
			my @inf=split;
			my @overlapperime=split /;/,$inf[1];
			for (@overlapperime)
			{
				$overlap{$_}=$inf[0];
			}
			$expression{$inf[0]}=0;
		}
		close OVERLAP;
	}
	open PRIMER,"$Bin/configd.dist -z $optsprimer{p}|" or die "please check the $optsprimer{p} file";
	while (<PRIMER>)
	{
		chomp;
		$total_primer++;
		my @inf=split;
		next if $inf[0] =~/^id/i;
		next if /^#/;
		my $len1=length($inf[1]);
		my $len2=length($inf[2]);
		next if $len1 <=13;
		next if $len2 <=13;
		####
		$primer{$inf[1]}{$inf[2]}=$inf[0];
		push@{$primer_re{$inf[0]}{1}},$inf[1];
		push@{$primer_re{$inf[0]}{2}},$inf[2];
		$primer_2_ref{F}{$inf[1]}=$inf[0];
		$primer_2_ref{R}{$inf[2]}=$inf[0];
		$expression{$inf[0]}=0;

		###make 8bp seed for search
		####pattern
		my $seed1=substr($inf[1],$len1-$optsprimer{seed},$optsprimer{seed});
		my $seed2=substr($inf[2],$len2-$optsprimer{seed},$optsprimer{seed});
		my $seed3=substr($inf[1],$len1-$optsprimer{seed}-$optsprimer{seed},$optsprimer{seed});
		my $seed4=substr($inf[2],$len2-$optsprimer{seed}-$optsprimer{seed},$optsprimer{seed});
		push @{$fseed{$seed1}},$inf[1];
		push @{$rseed{$seed2}},$inf[2];
		push @{$fseed1{$seed3}},$inf[1];
		push @{$rseed1{$seed4}},$inf[2];
		#####hash
		my $seed5=substr($inf[1],1,$optsprimer{seed});
		my $seed6=substr($inf[2],1,$optsprimer{seed});
		my $seed7=substr($inf[1],1+$optsprimer{seed},$optsprimer{seed});
		my $seed8=substr($inf[2],1+$optsprimer{seed},$optsprimer{seed});
		push @{$hash_fseed{$seed5}},$inf[1];
		push @{$hash_rseed{$seed6}},$inf[2];
		push @{$hash_fseed{$seed7}},$inf[1];
		push @{$hash_rseed{$seed8}},$inf[2];
	}
	close PRIMER;


	my %reads;
#################fastq read
	if ($optsprimer{fq1}=~/.gz$/)
	{
		open (FQ1,"gzip -dc $optsprimer{fq1}|") or die $!;
	}
	else
	{
		open (FQ1,"<$optsprimer{fq1}") or die $!;
	}
	if (defined$optsprimer{fq2} && -e "$optsprimer{fq2}")
	{
		if ($optsprimer{fq2}=~/.gz$/)
		{
			open (FQ2,"gzip -dc $optsprimer{fq2}|") or die $!;
		}
		else
		{
			open (FQ2,"<$optsprimer{fq2}") or die $!;
		}
	}
	else
	{
		if ($optsprimer{fq1}=~/.gz$/)
		{
			open (FQ2,"gzip -dc $optsprimer{fq1}|") or die $!;
		}
		else
		{
			open (FQ2,"<$optsprimer{fq1}") or die $!;
		}
	}
	open OUT1,"| gzip > $optsprimer{o}_1.fq.gz";

	open OUT2,"| gzip > $optsprimer{o}_2.fq.gz";

	open OUT3,"| gzip > $optsprimer{o}_single_1.fq.gz";

	open OUT4,"| gzip > $optsprimer{o}_single_2.fq.gz";

	if ($optsprimer{nonspe} eq "yes")
	{
		open OUT5,"| gzip > $optsprimer{o}_nonspe_1.fq.gz";

		open OUT6,"| gzip > $optsprimer{o}_nonspe_2.fq.gz";
	}

	my ($unmap,$total,$single,$map,$dimmer,$Non)=(0,0,0,0,0,0);
	my %unmap;
	my %fprimer_check;
	my %rprimer_check;
	my %Non;
##############
	while(<FQ1>)
	{
		if(($.-1)%4==0)
		{
			my ($cut_seq1,$cut_seq2,$cut_qual1,$cut_qual2);
			my $id1=$_;
			my $id2=<FQ2>;
			my $seq1=<FQ1>;
			my $seq2=<FQ2>;
			chomp $seq1;
			chomp $seq2;
			<FQ1>;<FQ2>;
			my $qual1=<FQ1>;
			my $qual2=<FQ2>;
			next unless (defined$seq1  && $seq1 ne "");
#			next unless (defined$seq2  && $seq2 ne "");
			if (!defined$optsprimer{fq2} )
			{
				$seq2="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
				$qual2="FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
			}
			chomp $qual1;chomp $qual2;chomp $id1;chomp $id2;
			$total++;
			my $read1_primer;my $read2_primer;my $real_primer1;my $real_primer2;
			my $match1=0;my $match2=0;my $match1_pos=length($seq1);my $match2_pos=length($seq2);
			##################
			my ($check_read1,$check_read2)=(0,0);
			#####################search fq1
			#####################hash
			my $read1_seed1=substr($seq1,1,$optsprimer{seed});
			my $read1_seed2=substr($seq1,1+$optsprimer{seed},$optsprimer{seed});
			if (defined$hash_fseed{$read1_seed1} || defined$hash_fseed{$read1_seed2})
			{
				my $read1_seed;
				$read1_seed=$read1_seed1 if defined$hash_fseed{$read1_seed1};
				$read1_seed=$read1_seed2 if defined$hash_fseed{$read1_seed2};
				for my $primer(@{$hash_fseed{$read1_seed}})
				{
					my $primer_check=0;
					my $len1=length($primer);
					my $mismatch=0;
					my $seq1_primer=substr($seq1,0,$len1);
					my @align=split//,($seq1_primer ^ $primer);
					for (@align)
					{
						$mismatch++ if ord($_);
					}
					if ($mismatch<=$optsprimer{m})
					{
						my $new_match=length($seq1_primer)-$mismatch;
						next if $new_match < $match1;
						$match1=$new_match;
						$read1_primer=$primer;
						$primer_check=1;
						my $tmp_primer=$seq1_primer;
						$fprimer_check{$primer}{$tmp_primer}+=1;
						$real_primer1=$primer;
						$cut_seq1=$seq1;
						$cut_qual1=$qual1;
					}
					$check_read1+=1 if $primer_check==1;
#				last if $primer_check==1;
				}
			}
			##################pattern
			if ($check_read1==0)
			{
				for my $seed1(keys %fseed)
				{
					my $primer_check=0;
					if ($seq1 =~ /(?:$seed1)/)
					{
						my $len2 = index($seq1,$seed1)+$optsprimer{seed};
						for my $primer(@{$fseed{$seed1}})
						{
							my $len1=length($primer);
							my $mismatch=0;
							if ($len2>=$len1)
							{
								my $seq1_primer=substr($seq1,$len2-$len1,$len1);
								my @align=split//,($seq1_primer ^ $primer);
								for (@align)
								{
									$mismatch++ if ord($_);
								}
								if ($mismatch<=$optsprimer{m})
								{
									my $new_match=length($seq1_primer)-$mismatch;
									my $new_matchpos=$len2-$len1;
									next if $new_match < $match1+1;
									next if $new_matchpos > $match1_pos+15;
									$match1_pos=$len2-$len1;
									$match1=$new_match;
									$read1_primer=$primer;
									$primer_check=1;
									my $tmp_primer=substr($seq1,0,$len2);
									$fprimer_check{$primer}{$tmp_primer}+=1;
									$real_primer1=$tmp_primer;
									$cut_seq1=substr($seq1,$len2-$len1);
									$cut_qual1=substr($qual1,$len2-$len1);
								}
							}
							else
							{
								next if $len2 ==0;
								my $seq1_primer=substr($seq1,0,$len2);
								my $cut_primer=substr($primer,$len1-$len2,$len2);
								my @align=split//,($seq1_primer ^ $cut_primer);
								for (@align)
								{
									$mismatch++ if ord($_);
								}
								if ($mismatch<=$optsprimer{m})
								{
									my $new_match=length($seq1_primer)-$mismatch;
									my $new_matchpos=0;
									next if $new_match < $match1+1;
									next if $new_matchpos > $match1_pos;
									$match1_pos=0;
									$match1=$new_match;
									$read1_primer=$primer;
									$primer_check=1;
									$fprimer_check{$primer}{$seq1_primer}+=1;
									$real_primer1=$seq1_primer;
									$cut_seq1=substr($seq1,0);
									$cut_qual1=substr($qual1,0);
								}
							}
						}
						$check_read1+=1 if $primer_check==1;
#					last if $primer_check==1;
					}
				}
			}
			#################################
			if ($check_read1==0)
			{
				for my $seed1(keys %fseed1)
				{
					my $primer_check=0;
					if ($seq1 =~ /(?:$seed1)/)
					{
						my $len2 = index($seq1,$seed1)+$optsprimer{seed};
						for my $primer(@{$fseed1{$seed1}})
						{
							my $len1=length($primer);
							my $mismatch=0;
							if ($len2+$optsprimer{seed}>=$len1)
							{
								my $seq1_primer=substr($seq1,$len2-$len1+$optsprimer{seed},$len1);
								my @align=split//,($seq1_primer ^ $primer);
								for (@align)
								{
									$mismatch++ if ord($_);
								}
								if ($mismatch<=$optsprimer{m})
								{
									my $new_match=length($seq1_primer)-$mismatch;
									my $new_matchpos=$len2-$len1+$optsprimer{seed};
									next if $new_match < $match1+1;
									next if $new_matchpos > $match1_pos+15;
									$match1_pos=$new_matchpos;
									$match1=$new_match;
									$read1_primer=$primer;
									$primer_check=1;
									my $tmp_primer=substr($seq1,0,$len2+$optsprimer{seed});
									$fprimer_check{$primer}{$tmp_primer}+=1;
									$real_primer1=$tmp_primer;
									$cut_seq1=substr($seq1,$len2-$len1+$optsprimer{seed});
									$cut_qual1=substr($qual1,$len2-$len1+$optsprimer{seed});
								}
							}
							else
							{
								next if $len2 ==0;
								my $seq1_primer=substr($seq1,0,$len2+$optsprimer{seed});
								my $cut_primer=substr($primer,$len1-$len2-$optsprimer{seed},$len2+$optsprimer{seed});
								my @align=split//,($seq1_primer ^ $cut_primer);
								for (@align)
								{
									$mismatch++ if ord($_);
								}
								if ($mismatch<=$optsprimer{m})
								{
									my $new_match=length($seq1_primer)-$mismatch;
									my $new_matchpos=0;
									next if $new_match < $match1+1;
									next if $new_matchpos > $match1_pos;
									$match1_pos=0;
									$match1=$new_match;
									$read1_primer=$primer;
									$primer_check=1;
									$fprimer_check{$primer}{$seq1_primer}+=1;
									$real_primer1=$seq1_primer;
									$cut_seq1=substr($seq1,0);
									$cut_qual1=substr($qual1,0);
								}
							}
						}
						$check_read1+=1 if $primer_check==1;
#					last if $primer_check==1;
					}
				}
			}

			####################################
			my $ummap_read1=substr($seq1,0,16);
#		$unmap{1}{$ummap_read1}=0 if $check_read1==0;
			#####################search fq2
			######################hash
			if (defined$optsprimer{fq2} and -e "$optsprimer{fq2}")
			{
				my $read2_seed1=substr($seq2,1,$optsprimer{seed});
				my $read2_seed2=substr($seq2,1+$optsprimer{seed},$optsprimer{seed});
				if (defined$hash_rseed{$read2_seed1} || defined$hash_rseed{$read2_seed2})
				{
					my $read2_seed;
					$read2_seed=$read2_seed1 if defined$hash_rseed{$read2_seed1};
					$read2_seed=$read2_seed2 if defined$hash_rseed{$read2_seed2};
					for my $primer(@{$hash_rseed{$read2_seed}})
					{
						my $primer_check=0;
						my $len2=length($primer);
						my $mismatch=0;
						my $seq2_primer=substr($seq2,0,$len2);
						my @align=split//,($seq2_primer ^ $primer);
						for (@align)
						{
							$mismatch++ if ord($_);
						}
						if ($mismatch<=$optsprimer{m})
						{
							my $new_match=length($seq2_primer)-$mismatch;
							next if $new_match < $match2;
							$match2=$new_match;
							$read2_primer=$primer;
							$primer_check=1;
							my $tmp_primer=$seq2_primer;
							$fprimer_check{$primer}{$tmp_primer}+=1;
							$real_primer2=$primer;
							$cut_seq2=$seq2;
							$cut_qual2=$qual2;
						}
						$check_read2+=1 if $primer_check==1;
#				last if $primer_check==1;
					}
				}
				#################pattern
				if ($check_read2==0)
				{
					for my $seed2(keys %rseed)
					{
						my $primer_check=0;
						if ($seq2 =~ /(?:$seed2)/)
						{
							my $len2 = index($seq2,$seed2)+$optsprimer{seed};
							for my $primer(@{$rseed{$seed2}})
							{
								my $len1=length($primer);
								my $mismatch=0;
								if ($len2>=$len1)
								{
									my $seq2_primer=substr($seq2,$len2-$len1,$len1);
									my @align=split//,($seq2_primer ^ $primer);
									for (@align)
									{
										$mismatch++ if ord($_);
									}
									if ($mismatch<=$optsprimer{m})
									{
										my $new_match=length($seq2_primer)-$mismatch;
										my $new_matchpos=$len2-$len1;
										next if $new_match < $match2+1;
										next if $new_matchpos > $match2_pos+15;
										$match2_pos=$len2-$len1;
										$match2=$new_match;
										$read2_primer=$primer;
										$primer_check=1;
										my $tmp_primer=substr($seq2,0,$len2);
										$rprimer_check{$primer}{$tmp_primer}+=1;
										$real_primer2=$tmp_primer;
										$cut_seq2=substr($seq2,$len2-$len1);
										$cut_qual2=substr($qual2,$len2-$len1);
									}
								}
								else
								{
									next if $len2 ==0;
									my $seq2_primer=substr($seq2,0,$len2);
									my $cut_primer=substr($primer,$len1-$len2,$len2);
									my @align=split//,($seq2_primer ^ $cut_primer);
									for (@align)
									{
										$mismatch++ if ord($_);
									}
									if ($mismatch<=$optsprimer{m})
									{
										my $new_match=length($seq2_primer)-$mismatch;
										my $new_matchpos=0;
										next if $new_match < $match2;
										next if $new_matchpos > $match2_pos;
										$match2_pos=0;
										$match2=$new_match;
										$read2_primer=$primer;
										$primer_check=1;
										$rprimer_check{$primer}{$seq2_primer}+=1;
										$real_primer2=$seq2_primer;
										$cut_seq2=substr($seq2,0);
										$cut_qual2=substr($qual2,0);
									}
								}
							}
							$check_read2+=1 if $primer_check==1;
#					last if $primer_check==1;
						}
					}
				}
				##############3
				if ($check_read2==0)
				{
					for my $seed2(keys %rseed1)
					{
						my $primer_check=0;
						if ($seq2 =~ /(?:$seed2)/)
						{
							my $len2 = index($seq2,$seed2)+$optsprimer{seed};
							for my $primer(@{$rseed1{$seed2}})
							{
								my $len1=length($primer);
								my $mismatch=0;
								if ($len2+$optsprimer{seed}>=$len1)
								{
									my $seq2_primer=substr($seq2,$len2-$len1+$optsprimer{seed},$len1);
									my @align=split//,($seq2_primer ^ $primer);
									for (@align)
									{
										$mismatch++ if ord($_);
									}
									if ($mismatch<=$optsprimer{m})
									{
										my $new_match=length($seq2_primer)-$mismatch;
										my $new_matchpos=$len2-$len1+$optsprimer{seed};
										next if $new_match < $match2;
										next if $new_matchpos > $match2_pos +15;
										$match2_pos=$new_matchpos;
										$match2=$new_match;
										$read2_primer=$primer;
										$primer_check=1;
										my $tmp_primer=substr($seq2,0,$len2+$optsprimer{seed});
										$fprimer_check{$primer}{$tmp_primer}+=1;
										$real_primer2=$tmp_primer;
										$cut_seq2=substr($seq2,$len2-$len1+$optsprimer{seed});
										$cut_qual2=substr($qual2,$len2-$len1+$optsprimer{seed});
									}
								}
								else
								{
									next if $len2 ==0;
									my $seq2_primer=substr($seq2,0,$len2+$optsprimer{seed});
									my $cut_primer=substr($primer,$len1-$len2-$optsprimer{seed},$len2+$optsprimer{seed});
									my @align=split//,($seq2_primer ^ $cut_primer);
									for (@align)
									{
										$mismatch++ if ord($_);
									}
									if ($mismatch<=$optsprimer{m})
									{
										my $new_match=length($seq2_primer)-$mismatch;
										my $new_matchpos=0;
										next if $new_match < $match2;
										next if $new_matchpos > $match2_pos;
										$match2_pos=0;
										$match2=$new_match;
										$read2_primer=$primer;
										$primer_check=1;
										$fprimer_check{$primer}{$seq2_primer}+=1;
										$real_primer2=$seq2_primer;
										$cut_seq2=substr($seq2,0);
										$cut_qual2=substr($qual2,0);
									}
								}
							}
							$check_read2+=1 if $primer_check==1;
#					last if $primer_check==1;
						}
					}
				}
				####
			}
			my $ummap_read2=substr($seq2,0,16);
#		$unmap{2}{$ummap_read2}=0 if $check_read2==0;
			######################################################
			################modify 20181106 最短比对长度
			if (defined$read1_primer)
			{
				my $primerl1=length$real_primer1;
				undef$read1_primer if $primerl1 / length($read1_primer) <$optsprimer{lower};
			}
			if (defined$read2_primer)
			{
				my $primerl2=length$real_primer2; 
				undef$read2_primer if $primerl2 / length($read2_primer) <$optsprimer{lower};
			}
			################去除引物
			my $tmpseq1len=length($cut_seq1);
			my $tmpseq2len=length($cut_seq2);
			if ($optsprimer{cut} =~/yes/i)
			{
				if (defined$read2_primer && defined$read1_primer)
				{

					my $reads1len;
					my $reads2len;

					if (length($read1_primer) > length($real_primer1))
					{
						$reads1len=length($real_primer1);
					}
					else
					{
						$reads1len=length($read1_primer);
					}
					if (length($read2_primer) > length($real_primer2))
					{
						$reads2len=length($real_primer2);
					}
					else
					{
						$reads2len=length($read2_primer);
					}

					$cut_seq1=substr($cut_seq1,$reads1len);
					$cut_seq2=substr($cut_seq2,$reads2len);
					$cut_qual1=substr($cut_qual1,$reads1len);
					$cut_qual2=substr($cut_qual2,$reads2len);
				}
			}
			##############################################
			$optsprimer{lr}=0 if $optsprimer{fqtype} eq "SE";
			if (defined$read2_primer && defined$read1_primer)
			{
				my $fprimer=$primer_2_ref{F}{$read1_primer};
				my $rprimer=$primer_2_ref{R}{$read2_primer};
				if (defined$primer{$read1_primer}{$read2_primer} && $tmpseq1len >= $optsprimer{l} && $tmpseq2len >= $optsprimer{lr})
				{
					$primer_result{$primer{$read1_primer}{$read2_primer}} ++ ;
					$subprimer{$primer{$read1_primer}{$read2_primer}}{$read1_primer}{$read2_primer}++ ;
					$primer_expand{"$primer{$read1_primer}{$read2_primer}"}{$read1_primer}{$real_primer1}++;
					$primer_expand{"$primer{$read1_primer}{$read2_primer}"}{$read2_primer}{$real_primer2}++;
					$map++;
					print OUT1 "$id1\n$cut_seq1\n+\n$cut_qual1\n";
					print OUT2 "$id2\n$cut_seq2\n+\n$cut_qual2\n";
				}
				elsif ($tmpseq1len < $optsprimer{l} || $tmpseq2len < $optsprimer{lr})
				{
					$error_primer{F}{$fprimer}{$real_primer1}++;
					$error_primer{R}{$rprimer}{$real_primer2}++;
					$dimmer++;
				}
				else
				{
					if (defined$overlap{$fprimer} && defined$overlap{$rprimer})
					{
						if ($overlap{$fprimer} eq $overlap{$rprimer})
						{
							$primer_result{$overlap{$fprimer}} ++ ;
							$subprimer{$overlap{$fprimer}}{$read1_primer}{$read2_primer}++ ;
							$primer_expand{$overlap{$fprimer}}{$read1_primer}{$real_primer1}++;
							$primer_expand{$overlap{$fprimer}}{$read2_primer}{$real_primer2}++;
							$map++;
							print OUT1 "$id1\n$cut_seq1\n+\n$cut_qual1\n";
							print OUT2 "$id2\n$cut_seq2\n+\n$cut_qual2\n";
						}
						else
						{
							$error_primer{F}{$fprimer}{$real_primer1}++;
							$error_primer{R}{$rprimer}{$real_primer2}++;
							$Non++;
							my ($nsfp,$nsrp)=("Unknow","Unknow");
							$nsfp="$primer_2_ref{F}{$read1_primer}_F" if defined$primer_2_ref{F}{$read1_primer};
							$nsrp="$primer_2_ref{R}{$read2_primer}_R" if defined$primer_2_ref{R}{$read2_primer};
							if ($optsprimer{nonspe} eq "yes")
							{
								print OUT5 "$id1\n$cut_seq1\n+\n$cut_qual1\n";
								print OUT6 "$id2\n$cut_seq2\n+\n$cut_qual2\n";
							}
							$Non{"$read1_primer\t$read2_primer\t$nsfp\t$nsrp"}++;
						}
					}
					else
					{
						$error_primer{F}{$fprimer}{$real_primer1}++;
						$error_primer{R}{$rprimer}{$real_primer2}++;
						$Non++;
						my ($nsfp,$nsrp)=("Unknow","Unknow");
						$nsfp="$primer_2_ref{F}{$read1_primer}_F" if defined$primer_2_ref{F}{$read1_primer};
						$nsrp="$primer_2_ref{R}{$read2_primer}_R" if defined$primer_2_ref{R}{$read2_primer};
						if ($optsprimer{nonspe} eq "yes")
						{
							print OUT5 "$id1\n$cut_seq1\n+\n$cut_qual1\n";
							print OUT6 "$id2\n$cut_seq2\n+\n$cut_qual2\n";
						}
						$Non{"$read1_primer\t$read2_primer\t$nsfp\t$nsrp"}++;
					}
				}
			}
			elsif(defined$read1_primer || defined$read2_primer)
			{
				if (defined$read1_primer && $tmpseq1len >=$optsprimer{l})
				{
					my $fprimer=$primer_2_ref{F}{$read1_primer};
					$error_primer{F}{$fprimer}{$real_primer1}++;
					$single++;
					print OUT3 "$id1\n$cut_seq1\n+\n$cut_qual1\n";
					print OUT4 "$id2\n$seq2\n+\n$qual2\n";
				}
				elsif(defined$read2_primer && $tmpseq2len >=$optsprimer{lr})
				{
					my $rprimer=$primer_2_ref{R}{$read2_primer};
					$error_primer{R}{$rprimer}{$real_primer2}++;
					$single++;
					print OUT3 "$id1\n$seq1\n+\n$qual1\n";
					print OUT4 "$id2\n$cut_seq2\n+\n$cut_qual2\n";
				}
				else
				{
					$dimmer++;
				}
			}
			else
			{
				$unmap++;
				$unmap{"$ummap_read1\t$ummap_read2"}++;
			}
		}
	}
	close FQ1;close FQ2;
	close OUT1;close OUT2;
	close OUT3;close OUT4;
	if ($optsprimer{nonspe} eq "yes")
	{
		close OUT5;close OUT6;
	}
	open STA,">$optsprimer{o}.primer_search.xls";
#print STA "primerID\treads_number\n";
	for my $id(keys %primer_result)
	{
		print STA "$id\t$primer_result{$id}\n";
		#####
		my $pri_num=0;
		for my $ref_primer(@{$primer_re{$id}{1}})
		{
			my $ref_primer2=${$primer_re{$id}{2}}[$pri_num];
			print STA "primer${pri_num}_F\t$ref_primer:$subprimer{$id}{$ref_primer}{$ref_primer2}\n";
			my $i=0;
			for (sort {$primer_expand{$id}{$ref_primer}{$b} <=>$primer_expand{$id}{$ref_primer}{$a}} keys %{$primer_expand{$id}{$ref_primer}})
			{
				$i++;
				print STA "$_:$primer_expand{$id}{$ref_primer}{$_}";
				last if $i>3;
				print STA "\t";
			}
			print STA "\n";
			print STA "primer${pri_num}_R\t$ref_primer2:$subprimer{$id}{$ref_primer}{$ref_primer2}\n";
			$i=0;
			for (sort {$primer_expand{$id}{$ref_primer2}{$b} <=>$primer_expand{$id}{$ref_primer2}{$a}} keys %{$primer_expand{$id}{$ref_primer2}})
			{
				$i++;
				print STA "$_:$primer_expand{$id}{$ref_primer2}{$_}";
				last if $i>3;
				print STA "\t";
			}
			print STA "\n";
			$i=0;
			print STA "Error primer F\t$ref_primer:$subprimer{$id}{$ref_primer}{$ref_primer2}\n";
			for (sort {$error_primer{F}{$id}{$b} <=> $error_primer{F}{$id}{$a}} keys %{$error_primer{F}{$id}})
			{
				$i++;
				print STA "$_:$error_primer{F}{$id}{$_}";
				last if $i > 3;
				print STA "\t";
			}
			print STA "\n";
			$i=0;
			print STA "Error primer R\t$ref_primer2:$subprimer{$id}{$ref_primer}{$ref_primer2}\n";
			for (sort {$error_primer{R}{$id}{$b} <=> $error_primer{R}{$id}{$a}} keys %{$error_primer{R}{$id}})
			{
				$i++;
				print STA "$_:$error_primer{R}{$id}{$_}";
				last if $i > 3;
				print STA "\t";
			}
			print STA "\n";
			$pri_num++;
		}
		####
		print STA "\n";
	}
	close STA;

	my $all_match;my $all_num;
	open IN,"<$optsprimer{o}.primer_search.xls";
	while (<IN>)
	{
		chomp;
		my @inf=split;
		next unless defined$inf[1];
		next unless $inf[1]=~/^\d+$/;
		$all_match+=$inf[1];
		$all_num++;
		$expression{$inf[0]}=$inf[1];
	}
	close IN;
	open EXP,">$optsprimer{o}.expression.xls";
	print EXP "ID\t${basename}_Count\t${basename}_Exp\t${basename}_Count/aver_dep\t${basename}_aver_dep\n";
	my $aver_depth;
	my $unpass_primer=0;
	$all_match||=1;
	$all_num||=0;
	for (keys %expression)
	{
		my $exp=sprintf "%.4f",$expression{$_}*$all_num*10000/$all_match;
		my $aver_dep=sprintf "%.4f",$expression{$_}*$all_num/$all_match;
		my $aver=sprintf "%.2f",$all_match/$total_primer;
		$aver_depth=$aver;
		$unpass_primer++ if $aver_dep <= 0.1;
		print EXP "$_\t$expression{$_}\t$exp\t$aver_dep\t$aver\n";
	}
	close EXP;

	my $top10=0;
	open Non,">$optsprimer{o}.Nonspecific.top10.xls" or die $!;
	for (sort {$Non{$b} <=> $Non{$a}} keys %Non)
	{
		$top10++;
		print Non "$_\t$Non{$_}\n";
		last if $top10 >=10;
	}
	close Non;
	if ($optsprimer{unmap_primer} eq "yes")
	{
		open OUT,">$optsprimer{o}.unmap.primer.xls" or die $!;
		for (sort{$unmap{$b}<=>$unmap{$a}} keys %unmap)
		{
			print OUT "$_\t$unmap{$_}\n";
		}
		close OUT;
	}

	open LOG,">$optsprimer{o}.primer_search.log" or die $!;
	print LOG "Total pair reads\t$total\n";
	print LOG "Pair end mapped reads\t$map\n";
	print LOG "Nonspecific amplification\t$Non\n";
	print LOG "Primer dimer\t$dimmer\n";
	print LOG "Single end mapped reads\t$single\n";
	print LOG "Unmap reads\t$unmap\n";
	print LOG "Average Depth\t$aver_depth\n";
	print LOG "Unqualified Primer\t$unpass_primer\n";
	my $map_ration= sprintf "%.2f",(( $map / $total) * 100 );
	print LOG "Mapped ratio\t$map_ration%\n";
	close LOG;
}




sub sam2result{
	my %optssam2result;
	($optssam2result{in},$optssam2result{out},$optssam2result{plasmid},$optssam2result{l},$optssam2result{p},$optssam2result{se},$optssam2result{IAC})=@_;


	$optssam2result{l}||=60;
	$optssam2result{se}||="no";

	my $reads=2;
	if ($optssam2result{se} eq "yes")
	{
		$reads=1;
	}
	
	my $samplename=basename($optssam2result{out});
	$samplename=~s/.result.xls//;
	
	my %DXcode;
	my %DXtarget;
	if (defined$optssam2result{IAC} && -e "$optssam2result{IAC}")
	{
		open IN,"<$optssam2result{IAC}" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split /\t/,$_;
			next unless $inf[0] eq $samplename;
			$inf[2]="DX0000";
			$DXcode{$inf[0]}=$inf[2];
			if (-e "$Bin/../database/dxdir/$inf[2].result")
			{
				open RE,"$Bin/.configx.dist -z $Bin/../database/dxdir/$inf[2].result|" or die $!;
				while (<RE>)
				{
					chomp;
					my @tmp=split /\t/,$_;
					$DXtarget{$inf[0]}{$tmp[0]}=1;
				}
				close RE;
			}
		}
		close IN;
	}
	###############
	my %hash;my %plasmid;my @sort;
	my %sorthash;
	my %sortuniq;
	if (defined$optssam2result{p})
	{
		open IN,"<$optssam2result{p}" or die $!;
		while (<IN>)
		{
			chomp;
			next unless /null/;
			my @inf=split; 
			my @name=split /,/,$inf[1];
			my $name="$name[0],$name[1]";
			next if defined$sortuniq{$name};
			push @sort,$name;
			$sortuniq{$name}=1;
		}
		close IN;
	}

	my %plasmid_seq;
	open IN,"<$optssam2result{plasmid}" or die $!;
	while (<IN>)
	{
		next if /^@/;
		chomp;
		my @inf=split;
		next if $inf[2] eq "*";
		my $len=length($inf[9]);
		next if $len <= $optssam2result{l};
		my ($S,$D,$H)=(0,0,0);
		if ($inf[5]=~/(\d+)S/)
		{
			$S=$1;
			while ($inf[5]=~/(\d+)S/g)
			{
				$S=$1 if $1 > $S;
			}
		}
		if ($inf[5]=~/(\d+)D/)
		{
			$D=$1;
			while ($inf[5]=~/(\d+)D/g)
			{
				$D=$1 if $1 > $D;
			}
		}
		if ($inf[5] =~ /(\d+)H/)
		{
			$H=$1;
			while ($inf[5]=~/(\d+)H/g)
			{
				$H=$1 if $1 > $H;
			}
		}
		next if $S > 7 || $D > 7 || $H > 7;
		$plasmid_seq{$inf[0]}+=1;
	}
	close IN;

	my %log;

	my %aligncheck;
	my %alignmatch;
	open IN,"<$optssam2result{in}" or die $!;
	while(<IN>)
	{
		next if /^@/;
		chomp;
		my @inf=split;
		$log{"total"}++;
		next if $inf[2] eq "*";
		$log{"allalign"}++;
		(my $name=$inf[2])=~s/,\d+$//;
		$sorthash{$name}=1;
		$aligncheck{$inf[0]}{$name}++;
		$alignmatch{$inf[0]}{"len"}+=length($inf[9]);
		while ($inf[5]=~/(\d+)M/g)
		{
			$alignmatch{$inf[0]}{"match"}+=$1;
		}
		while ($_=~ /NM:i:(\d+)/g)
		{
			$alignmatch{$inf[0]}{"mismatch"}+=$1;
		}
		$alignmatch{$inf[0]}{"match"}||=0;
		$alignmatch{$inf[0]}{"mismatch"}||=0;
		if (defined$plasmid_seq{$inf[0]} && $plasmid_seq{$inf[0]} >=$reads)
		{
			$plasmid{$name}{$inf[0]}+=1;
			$log{"align"}+=1;
		}
		else
		{
			if ($aligncheck{$inf[0]}{$name} == $reads)
			{
				next if ($alignmatch{$inf[0]}{"match"}-$alignmatch{$inf[0]}{"mismatch"})/$alignmatch{$inf[0]}{"len"} <0.9;
				$hash{$name}{$inf[0]}=$reads;
				$log{"align"}+=$reads;
			}
		}
	}
	close IN;

	unless (defined$sort[0])
	{
		for (sort keys %sorthash)
		{
			push @sort,$_ if /内参/;
		}
		for (sort keys %sorthash)
		{
			push @sort,$_ unless (/内参/ || /IAC/i); 
		}
	}

	open OUT,">$optssam2result{out}" or die $!;
	for (my $i=0;$i<@sort;$i++)
	{
		if(exists $hash{$sort[$i]})
		{
			my $num=0;
			foreach my $key2(keys %{$hash{$sort[$i]}})
			{
				if($hash{$sort[$i]}{$key2}==$reads)
				{
					$num++;
				}
			}
			my $tmptarget=(split/,/,$sort[$i])[0];
			if (defined$DXcode{$samplename})
			{
				if ($DXcode{$samplename} eq "DX0000")
				{
					$num=$num;
				}
				elsif(defined$DXtarget{$samplename}{$tmptarget} || $sort[$i]=~/内参/ || $sort[$i]=~/errorbarcode/)
				{
					$num=$num;
				}
				else
				{
					$num=0;
				}
			}
			else
			{
				$num=$num;
			}
			print OUT "$sort[$i]\t$num\n";
		}
		else
		{
			print OUT "$sort[$i]\t0\n";
		}
		#print OUT "$_\t$hash{$_}\n";
	}
	close OUT;

	open OUT,">$optssam2result{out}.plasmid" or die $!;
	for (my $i=0;$i<@sort;$i++)
	{
		if(exists $plasmid{$sort[$i]})
		{
			my $num=0;
			foreach my $key2(keys %{$plasmid{$sort[$i]}})
			{
				if($plasmid{$sort[$i]}{$key2}==$reads)
				{
					$num++;
				}
			}
			my $tmptarget=(split/,/,$sort[$i])[0];
			if (defined$DXcode{$samplename})
			{
				if ($DXcode{$samplename} eq "DX0000")
				{
					$num=$num;
				}
				elsif(defined$DXtarget{$samplename}{$tmptarget} || $sort[$i]=~/内参/ || $sort[$i]=~/errorbarcode/)
				{
					$num=$num;
				}
				else
				{
					$num=0;
				}
			}
			else
			{
				$num=$num;
			}
			print OUT "$sort[$i]\t$num\n";
		}
		else
		{
			print OUT "$sort[$i]\t0\n";
		}
		#print OUT "$_\t$plasmid{$_}\n";
	}
	close OUT;
	open OUT,">$optssam2result{out}.align.log" or die $!;
	$log{total}||=0;
	$log{allalign}||=0;
	$log{"align"}||=0;
	print OUT "Total reads num\t$log{total}\nAll align reads\t$log{allalign}\nReal align reads\t$log{align}\n";
	close OUT;
}

sub bamfilter{
	my ($bamfile,$pre,$samtools)=@_;

	my %hashlow;
	my %hashtotal;
	open IN,"$samtools view $bamfile|";
	while(<IN>)
	{
		my $M=0;
		my $S=0;
		my $I=0;
		my $D=0;
		my @ar=split;
		$hashtotal{$ar[0]}++;
		while($ar[5]=~/(\d+)S/g)
		{
			$S+=$1;
		}
		while($ar[5]=~/(\d+)I/g)
		{
			$I+=$1;
		}
		while($ar[5]=~/(\d+)I/g)
		{
			$D+=$1;
		}
		if($ar[12]=~/MD:Z:(.*)/)
		{
			my $seq=$1;
			while($seq=~/(\d+)/g)
			{
				$M+=$1;
			}
			$M=$M-$S-$I-$D;
			if($M<length($ar[9])*0.93) #正确比对上的小于93%被过滤掉
			{
				$hashlow{$ar[0]}++;#基因组层面正确配对小于90%被标记
				#print "soga$seq\t$M\t$ar[0]\n";
			}
			#print "$seq\t$M\n";
		}
	}
	close IN;

	open OUT,"> $pre.filter.sam";
	open IN,"$samtools view -h $bamfile|";
	while(<IN>)
	{
		my $li=$_;
		my @ar=split;
		next if(exists $hashlow{$ar[0]});
		print OUT "$li";

	}
	close IN;

	open OUT,"> $pre.bamfilter.sta";
	my $total=keys %hashtotal;
	my $filter=keys %hashlow;
	print OUT "rawreads\t$total\nlowread\t$filter\b";
	close OUT;

	system("$samtools view -S $pre.filter.sam -b > $pre.filter.bam && rm $pre.filter.sam");
}

