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

use Cwd qw/abs_path/;
use FindBin qw/$Bin $Script/;
use File::Basename qw/dirname basename/;
use lib "$Bin";
use Spreadsheet::XLSX; 
use File::Basename;
use File::Path;
#use FileHandle;

#use List::Util qw/max min sum/;
#use List::MoreUtils qw/any all uniq/;

#use threads;
#use threads::shared;
#use Thread::Semaphore;
#use Thread::Queue;

my %optsnew;
use Getopt::Long;
GetOptions (\%optsnew,
	"BYdirin=s",
	"BYdirIAC=s",
	"BYdirinf=s",
	"sample=s",
	"database=s",
	"out=s",
	"out2=s",
	"model=s",
	"help" => \&help,
);

###########################################################
## Your code go here
$optsnew{model}||="full";
&BY_dir($optsnew{BYdirin}, "$optsnew{out}/pathogen", 1, "", $optsnew{BYdirIAC}, $optsnew{BYdirinf});
if ($optsnew{model} eq "silent")
{
	system("rm -rf $optsnew{out}/pathogen.plasmid $optsnew{out}/pathogen.raw $optsnew{out}/pathogen.uniq $optsnew{out}/pathogen.IAC $optsnew{out}/pathogen.barcode.pollute $optsnew{out}/pathogen.sampleiac $optsnew{out}/pathogen.rawdorclient $optsnew{out}/pathogen.h2o $optsnew{out}/pathogen.errorbarcode $optsnew{out}/pathogen.data.stat.xls $optsnew{out}/pathogen_sample.data.stat.xls.forclient $optsnew{out}/pathogen.raw1 $optsnew{out}/pathogen.barcode.pollute1 $optsnew{out}/pathogen.h2o1 $optsnew{out}/pathogen.uniq1 $optsnew{out}/pathogen.IAC1 $optsnew{out}/pathogen.sampleiac1 $optsnew{out}/pathogen.plasmid1");
}
&sample_positive_stat("$optsnew{out}/pathogen.positive", $optsnew{sample}, "$optsnew{out}/sample.positive.stat.xlsx");

#&positive2reporttxt("$optsnew{out}/pathogen.positive",$optsnew{sample},"$optsnew{database}/dxdir/", "$optsnew{database}/restriction.txt", "$optsnew{database}/explain.txt", "$optsnew{database}/drug_annot.txt","$optsnew{out}/report.stat.txt" ,"$optsnew{database}/pathogen.type.txt");
&positive2reporttxtv11("$optsnew{out}/pathogen.positive",$optsnew{sample},"$optsnew{database}/dxdir/", "$optsnew{database}/restriction_v2.txt", "$optsnew{out2}/report.stat.txt", "$optsnew{out}/pathogen_sample.data.stat.xls");
#&make_mutation_stat($optsnew{BYdirin},"$optsnew{out2}/JH.mutation.stat.xls");

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

sub BY_dir{
	my %opts;
	($opts{in},$opts{out},$opts{n},$opts{x},$opts{IAC},$opts{sample})=@_;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $reporttime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);

	my $xiaolv=$opts{x};
	$opts{n}||=1;

	my %scaler;
	my %average;
	my $samplenum;
	my %tartotal;
	my %errorbarcode;



	my @file=`ls $opts{in}/*/align/*result.xls`;
	my @plasmid=`ls $opts{in}/*/align/*result.xls.plasmid`;
	my @resultsort=`ls $opts{in}/*/align/*.result.xls`;

	my %sampletype;
	my %sampleerror;
	open IN,"<$opts{sample}" or die $!;
	while (<IN>)
	{
		next if $.==1;
		chomp;
		my @inf=split /\t/,$_;
		$inf[1]||="-";
		$sampletype{$inf[0]}=$inf[1];
		$sampletype{$inf[0]}="痰液" if $inf[1]=~/痰液/ || $inf[1] eq "-";
		$sampleerror{$inf[0]}=$inf[4] if defined$inf[4] && $inf[4] ne "-";
	}
	close IN;

	die "There are no file for stat!!\n" unless $file[0];
	my $i=0;my @sort;my %hash;my %plasmid;
	my %hashclient;
	my %sampletar;
	for (@file)
	{
		$i++;
		chomp;
		my $name=basename($_);
		$name=~s/.result.xls$//;
		open IN,"<$_" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split;
			my @primer=(split /,/,$inf[0]);
			my $primer=join ",",($primer[0],$primer[1]);
			if ($i==1)
			{
				push @sort,$primer;
			}
			else
			{
				push @sort,$primer unless defined$tartotal{$inf[0]};
			}
			$hash{$name}{$primer}=$inf[1];
			$hashclient{$name}{$primer[0]}+=$inf[1];
			$sampletar{$name}{$primer}=1;
			$hash{$name}{total}+=$inf[1];
			$hash{$name}{"定量内参"}+=$inf[1] if $inf[0]=~/内参/;
			$hashclient{$name}{"定量内参"}+=$inf[1] if $inf[0]=~/内参/;
			if ($primer=~/errorbarcode/)
			{
				(my $code=$primer)=~s/.errorbarcode//;
				if (defined$errorbarcode{$name}{barcode})
				{
					my $oldprimer=$errorbarcode{$name}{barcode}.",errorbarcode";
					$errorbarcode{$name}{barcode}=$code if $inf[1] > $hash{$name}{$oldprimer};
				}
				else
				{
					$errorbarcode{$name}{barcode}=$code if $inf[1] >0;
				}
			}
			unless (defined$tartotal{$primer})
			{
				$tartotal{$primer}=$inf[1];
			}
			else
			{
				$tartotal{$primer}=$inf[1] if $inf[1]>$tartotal{$primer};
			}
		}
		close IN;
	}
	my @newsort;
	push @newsort,"定量内参";
	for (@sort)
	{
		push @newsort,$_ if $_=~/内参/;
	}
	for (@sort)
	{
		push @newsort,$_ unless $_=~/内参/;
	}
	@sort=@newsort;



	for (@plasmid)
	{
		$i++;
		chomp;
		my $name=basename($_);
		$name=~s/.result.xls.plasmid$//;
		open IN,"<$_" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split;
			my @primer=(split /,/,$inf[0]);
			my $primer=join ",",($primer[0],$primer[1]);
			$plasmid{$name}{$primer}=$inf[1];
			$plasmid{$name}{total}+=$inf[1];
		}
		close IN;
	}

	my %sampletarget;
	my %alltar;
	my @othersort;
	for (@sort)
	{
		my $name=(split /,/,$_)[0];
		$alltar{$name}=1;
	}
	for my $sample(sort keys %hash)
	{
		for my $file(@resultsort)
		{
			chomp $file;
			my $filename=(split /\//,$file)[-3];
			next unless $sample eq $filename;
			open IN,"<$file" or die $!;
			while (<IN>)
			{
				chomp;
				my @inf=split;
				$sampletarget{$sample}{$inf[0]}=1;
				push @othersort,$inf[0] unless defined$alltar{$inf[0]};
				$alltar{$inf[0]}=1;
			}
			close IN;
		}
	}
####
	my %IAC;
	if (defined$opts{IAC})
	{
		open IN,"<$opts{IAC}" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split /\t/,$_;
			$IAC{$inf[0]}=$inf[1];
		}
		close IN;
	}




############zhili
	open OUT ,">$opts{out}.plasmid" or die $!;
	print OUT "样本";
	for (sort{$b cmp $a} keys %plasmid)
	{
		print OUT "\t$_";
	}
	print OUT "\n";
	for my $tar(@sort)
	{
		print OUT "$tar";
		for (sort{$b cmp $a} keys %plasmid)
		{
			$plasmid{$_}{$tar}||=0;
			print OUT "\t$plasmid{$_}{$tar}";
		}
		print OUT "\n";
	}
	close OUT;


#####JHNY
	open OUT,">$opts{out}.JHNY.txt" or die $!;
	for (sort keys %sampletype)
	{
		my $samplename=$_;
		if (-s "$opts{in}/$_/align/$_.JH.hotspot_positive.xls")
		{
			open IN,"<$opts{in}/$_/align/$_.JH.hotspot_positive.xls" or die $!;
			while (<IN>)
			{
				chomp;
				next if $.==1;
				print OUT "$samplename\t$_\n";
			}
			close IN;
		}
	}
	close OUT;

###############原始数据\归一化\水过滤\扩增效率校正\IAC归一化
	open OUT1 ,">$opts{out}.raw" or die $!;
	open OUT2 ,">$opts{out}.barcode.pollute" or die $!;
	open OUT3 ,">$opts{out}.h2o" or die $!;
	open OUT4 ,">$opts{out}.uniq" or die $!;
	open OUT5 ,">$opts{out}.IAC" or die $!;
	open OUT6 ,">$opts{out}.sampleiac" or die $!;
	open OUT7 ,">$opts{out}.rawforclient" or die $!;
	open OUT8 ,">$opts{out}.errorbarcode" or die $!;
	print OUT1 "样本";print OUT2 "样本";print OUT3 "样本";print OUT4 "样本";print OUT5 "样本";print OUT6 "样本";print OUT7 "样本";print OUT8 "样本";
	for (sort{$b cmp $a} keys %hash)
	{

		print OUT1 "\t$_";print OUT2 "\t$_";print OUT3 "\t$_";print OUT4 "\t$_";print OUT5 "\t$_";print OUT6 "\t$_";print OUT7 "\t$_";print OUT8 "\t$_";
	}
	print OUT1 "\n";print OUT2 "\n";print OUT3 "\n";print OUT4 "\n";print OUT5 "\n";print OUT6 "\n";print OUT7 "\n";print OUT8 "\n";

######
	my $water_name="water_undefined";
	my @water_name;
	for (sort{$b cmp $a} keys %hash)
	{
		push @water_name,$_ if $_=~ /water|^shui|^h2o/i;
		$samplenum++ unless $_=~ /water|^shui|^h2o/i;
	}
	push @water_name,$water_name unless defined$water_name[0];
####

#####防错标签
	print OUT8 "实际标签";
	for (sort{$b cmp $a} keys %hash)
	{
		$errorbarcode{$_}{barcode}||="-";
		print OUT8 "\t$errorbarcode{$_}{barcode}";
	}
	print OUT8 "\n";
	print OUT8 "登记标签";
	for (sort{$b cmp $a} keys %hash)
	{
		$sampleerror{$_}||="-";
		print OUT8 "\t$sampleerror{$_}";
	}
	print OUT8 "\n";

	for my $tar(@sort)
	{
		next unless $tar=~/errorbarcode/;
		(my $code=$tar)=~s/.errorbarcode//;
		print OUT8 "$code";
		for (sort{$b cmp $a} keys %hash)
		{
			print OUT8 "\t$hash{$_}{$tar}";
		}
		print OUT8 "\n";
	}

######barcode漂移过滤
	for my $tar(@sort)
	{
		next if $tar=~/errorbarcode/;
		for (sort{$b cmp $a} keys %hash)
		{
			##判断水样本 不做barcode漂移去除
			my $watercheck=0;
			$hash{$_}{$tar}||=0;
			$hash{$_}{total}||=1;
			$hash{$_}{"定量内参"}||=1;
			$tartotal{$tar}||=0;
			for my $watername(@water_name){$watercheck++ if $_ eq $watername};
			$watercheck=0;
			if ($watercheck >0  || $tar =~/内参/i)
			{
				$scaler{barcode}{$_}{$tar}=$hash{$_}{$tar};
			}
			else
			{
				$scaler{barcode}{$_}{$tar}=int($hash{$_}{$tar} - ($tartotal{$tar} * 0.001));#barcode漂移1/1000
				$scaler{barcode}{$_}{$tar}=0 if $scaler{barcode}{$_}{$tar} <0;
			}
			$hash{$_}{newtotal}+=$scaler{barcode}{$_}{$tar};
		}
	}

	my %uniq;
	my %water;
########
	for my $tar(@sort)
	{
		next if $tar=~/errorbarcode/;
		print OUT1 "$tar";print OUT2 "$tar";print OUT3 "$tar";print OUT5 "$tar";print OUT6 "$tar";
		my $wuran=0;
		my $pathogenname=(split /,/,$tar)[0];
		for (sort{$b cmp $a} keys %hash)
		{
			$hash{$_}{$tar}||=0;
			$hash{$_}{total}||=1;
			$hash{$_}{"定量内参"}||=1;
			$hash{$_}{newtotal}||=1;
		}
		for (sort{$b cmp $a} keys %hash)
		{
			####水过滤
			if ($tar =~ /内参/i)
			{
				$scaler{h2o}{$_}{$tar}=$scaler{barcode}{$_}{$tar};
			}
			else
			{
				my $watern=$water_name[0];
				my $tmpwater=$scaler{barcode}{$water_name[0]}{$tar};
				for my $water(@water_name)
				{
					$scaler{barcode}{$water}{$tar}||=0;
					$scaler{barcode}{$water}{"定量内参"}||=1;
					if ($scaler{barcode}{$water}{$tar}>$tmpwater)
					{
						$tmpwater=$scaler{barcode}{$water}{$tar};
						$watern=$water;
					}
					$wuran=$scaler{barcode}{$water}{$tar}/$scaler{barcode}{$water}{"定量内参"} if $scaler{barcode}{$water}{$tar}/$scaler{barcode}{$water}{"定量内参"} > $wuran;
				}
				$scaler{barcode}{$_}{"定量内参"}||=1;
				if (defined$plasmid{$_}{$tar} && $plasmid{$_}{$tar} > 0 && defined$plasmid{$watern}{$tar} && $plasmid{$watern}{$tar} >0)
				{
					$scaler{h2o}{$_}{$tar}=$scaler{barcode}{$_}{$tar}-50*($scaler{barcode}{$watern}{$tar}*$plasmid{$_}{$tar}/$plasmid{$watern}{$tar});#10倍水
					$scaler{h2o}{$_}{$tar}=0 if $scaler{h2o}{$_}{$tar}<0;
				}
				else
				{
					$scaler{h2o}{$_}{$tar}=$scaler{barcode}{$_}{$tar}-($wuran * 50 *  $scaler{barcode}{$_}{"定量内参"}); #10倍水
					$scaler{h2o}{$_}{$tar}=0 if $scaler{h2o}{$_}{$tar}<0;
				}
			}
			##########################
			#######菌株唯一引物
			if (defined$scaler{uniq}{$_}{$pathogenname})
			{
				if ($scaler{h2o}{$_}{$tar}>$scaler{uniq}{$_}{$pathogenname})
				{
					$scaler{uniq}{$_}{$pathogenname}=$scaler{h2o}{$_}{$tar};
					$uniq{$_}{$pathogenname}=$tar;
				}
			}
			else
			{
				$scaler{uniq}{$_}{$pathogenname}=$scaler{h2o}{$_}{$tar};
				$uniq{$_}{$pathogenname}=$tar;
			}
			#####################################################
			####PCR管内定量
			if (defined$plasmid{$_}{$tar} && $plasmid{$_}{$tar} > 0)
			{
				$scaler{IAC}{$_}{$tar}=$scaler{h2o}{$_}{$tar}*50/$plasmid{$_}{$tar};
			}
			else
			{
				$scaler{h2o}{$_}{"定量内参"}||=1;
				$scaler{IAC}{$_}{$tar}=$scaler{h2o}{$_}{$tar}*50/$scaler{h2o}{$_}{"定量内参"};
			}
##################
			#样本定量
			$sampletype{$_}||="痰液";
			if ($sampletype{$_}=~/肺泡灌洗液|脑脊液|血|胸水|腹水|胸腹水|尿液/i)
			{
				$IAC{$_}||=1;
				$scaler{sampleiac}{$_}{$tar}=$scaler{IAC}{$_}{$tar}*50/(5*$IAC{$_});
			}
			elsif($sampletype{$_}=~/房水|玻璃体液/)
			{
				$scaler{sampleiac}{$_}{$tar}=$scaler{IAC}{$_}{$tar}*20/(5*0.1);
			}
			else
			{
				$scaler{sampleiac}{$_}{$tar}=$scaler{IAC}{$_}{$tar}*80*2/(5*0.2);
			}

			##################
			$scaler{h2o}{$_}{$tar}=int($scaler{h2o}{$_}{$tar});
			$scaler{h2o}{$_}{total}+=$scaler{h2o}{$_}{$tar};
			$scaler{uniq}{$_}{$pathogenname}=int($scaler{uniq}{$_}{$pathogenname});
			$scaler{IAC}{$_}{$tar}=int($scaler{IAC}{$_}{$tar});
			$scaler{sampleiac}{$_}{$tar}=int($scaler{sampleiac}{$_}{$tar});
			$average{$tar}+=$scaler{h2o}{$_}{$tar} unless $_=~ /water|^shui|^h2o/i;
			print OUT1 "\t$hash{$_}{$tar}";
			print OUT2 "\t$scaler{barcode}{$_}{$tar}";
			print OUT3 "\t$scaler{h2o}{$_}{$tar}";
			#print OUT4 "\t$scaler{uniq}{$_}{$pathogenname}";
			print OUT5 "\t$scaler{IAC}{$_}{$tar}";
			print OUT6 "\t$scaler{sampleiac}{$_}{$tar}";
		}
		print OUT1 "\n";print OUT2 "\n";print OUT3 "\n";print OUT5 "\n";print OUT6 "\n";
	}


	my %uniqout;
	for my $tar(@sort)
	{
		my $pathogenname=(split /,/,$tar)[0];
		for (sort{$b cmp $a} keys %hash)
		{
			if (defined $uniq{$_}{$pathogenname} &&  $uniq{$_}{$pathogenname} eq $tar)
			{
				$uniqout{$pathogenname}{$_}=$scaler{h2o}{$_}{$tar};
			}
#		else
#		{
#			$uniqout{$pathogenname}{$_}="-";
#		}
		}
	}

	my %q;
	for my $tar(@sort)
	{
		my $newtar=(split /,/,$tar)[0];
		next if defined$q{$newtar};
		$q{$newtar}=1;
		print OUT4 "$newtar";
		for (sort{$b cmp $a} keys %{$uniqout{$newtar}})
		{
			$uniqout{$newtar}{$_}||=0;
			print OUT4 "\t$uniqout{$newtar}{$_}";
		}
		print OUT4 "\n";
	}
	close OUT1;close OUT2;close OUT4;close OUT5;

	my %p;
	for my $tar(@sort)
	{
		my $newtar=(split /,/,$tar)[0];
		next if defined$p{$newtar};
		next if $tar=~/扩增内参/;
		$p{$newtar}=1;
		print OUT7 "$newtar";
		for (sort{$b cmp $a} keys %hash)
		{
			print OUT7 "\t$hashclient{$_}{$newtar}";
		}
		print OUT7 "\n";
	}

	open OUT ,">$opts{out}.positive" or die $!;
	print OUT "样本\t检测目标\tIAC\t引物扩增效率校正结果\t置信度\t相对含量\n";
	for my $tar(@sort)
	{
		next if $tar=~/扩增内参/;
		next if $tar=~/乙肝病毒/;
		next if $tar=~/errorbarcode/;
		for my $sample(sort{$b cmp $a} keys %hash)
		{
			next if $sample=~ /water|^shui|^h2o/i;
			my $pathogenname=(split /,/,$tar)[0];
#		next unless defined$sampletarget{$sample}{$pathogenname};
			$scaler{h2o}{$sample}{total}||=1;
			my $percentage=sprintf "%.3f",$scaler{h2o}{$sample}{$tar}*100/$scaler{h2o}{$sample}{total};
			next unless $uniq{$sample}{$pathogenname} eq $tar;
			if($scaler{h2o}{$sample}{total} > 100000)
			{
				if ($sampletype{$sample} =~ /痰液|肺泡灌洗液|咽拭子|鼻咽拭子/)
				{
					if ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 200)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t高\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 60)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 0)
					{
						if ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > $average{$tar}*3/$samplenum && $scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 5) #平均深度3倍，且大于5条
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
						}
						else
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t低\t$percentage\n";
						}
					}
					else
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t无\t$percentage\n";
					}
				}
				else
				{
					if ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 200)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t高\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 50)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 0)
					{

						if ($scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > $average{$tar}*3/$samplenum && $scaler{h2o}{$sample}{$tar}*100000/$scaler{h2o}{$sample}{total} > 5) #平均深度3倍，且大于5条
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
						}
						else
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t低\t$percentage\n";
						}
					}
					else
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t无\t$percentage\n";
					}
				}
			}
			else
			{
				if ($sampletype{$sample} =~/痰液|肺泡灌洗液/)
				{
					if ($scaler{h2o}{$sample}{$tar} > 200)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t高\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar} > 60)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
					}
					elsif  ($scaler{h2o}{$sample}{$tar} > 0)
					{
						if ($scaler{h2o}{$sample}{$tar} > $average{$tar}*3/$samplenum && $scaler{h2o}{$sample}{$tar} > 5) #平均深度3倍，且大于5条
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
						}
						else
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t低\t$percentage\n";
						}
					}
					else
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t无\t$percentage\n";
					}
				}
				else
				{
					if ($scaler{h2o}{$sample}{$tar} > 200)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t高\t$percentage\n";
					}
					elsif ($scaler{h2o}{$sample}{$tar} > 50)
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
					}
					elsif  ($scaler{h2o}{$sample}{$tar} > 0)
					{
						if ($scaler{h2o}{$sample}{$tar} > $average{$tar}*3/$samplenum && $scaler{h2o}{$sample}{$tar} > 5) #平均深度3倍，且大于5条
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t中\t$percentage\n";
						}
						else
						{
							print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t低\t$percentage\n";
						}
					}
					else
					{
						print OUT "$sample\t$pathogenname\t$scaler{sampleiac}{$sample}{$tar}\t$scaler{h2o}{$sample}{$tar}\t无\t$percentage\n";
					}
				}
			}
		}
	}
	for my $tar(@othersort)
	{
		next if $tar=~/内参/;
		next if $tar=~/乙肝病毒/;
		for my $sample(sort{$b cmp $a} keys %hash)
		{
			next if $sample=~ /water|^shui|^h2o/i;
			next unless defined$sampletarget{$sample}{$tar};
			my $pathogenname=(split /,/,$tar)[0];
			my $percentage=sprintf "%.3f",0*100/$scaler{h2o}{$sample}{total};
#		print OUT "$sample\t$pathogenname\t0\t0\t无\t$percentage\n";
		}
	}
	close OUT;

################################
	open OUT,">$opts{out}.data.stat.xls" or die $!;
	print OUT "Sample\n实际标签\n登记标签\nTotal pair reads\nhigh quality pair reads(Q30)\nAdapter filter\nDimmer\nNonspecific amplification\nSingle end mapped reads\nUnknow Primer\nClean pair reads\nTotal Align\nReal Align\n";
	close OUT;
	for (sort{$b cmp $a} keys %hash)
	{
		next if /water_undefined/;
		open OUT,">$opts{out}.tmp.head" or die $!;
		print OUT "Sample\t$_\n实际标签\t$errorbarcode{$_}{barcode}\n登记标签\t$sampleerror{$_}\n";
		close OUT;
		&quality_stat("$opts{in}/$_/clean/$_.json","$opts{in}/$_/clean/${_}.primer_search.log","$opts{in}/$_/align/${_}.result.xls.align.log","$opts{out}.stat.xls.tmp");
		system("cat $opts{out}.tmp.head $opts{out}.stat.xls.tmp > $opts{out}.stat.xls.tmp2");
		&add_desc("$opts{out}.stat.xls.tmp2","1","$opts{out}.data.stat.xls","0"," $opts{out}.stat.xls.tmp1");
		system("mv -f $opts{out}.stat.xls.tmp1 $opts{out}.data.stat.xls");
		system("rm -rf $opts{out}.stat.xls.tmp2");
		&merge_sampleinf_datainf("$opts{sample}","$opts{out}.data.stat.xls","$opts{out}_sample.data.stat.xls");
		system("cp $opts{out}.data.stat.xls $opts{out}_sample.data.stat.xls") unless defined$opts{sample};
		system("rm -rf $opts{out}.stat.xls.tmp $opts{out}.tmp.head");
	}
	system("awk 'NR==1||NR==2||NR==3||NR==4||NR==5||NR==13||NR==15||NR==18||NR==19' $opts{out}_sample.data.stat.xls > $opts{out}_sample.data.stat.xls.forclient");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.raw","$opts{out}.raw1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.barcode.pollute","$opts{out}.barcode.pollute1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.h2o","$opts{out}.h2o1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.uniq","$opts{out}.uniq1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.IAC","$opts{out}.IAC1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.sampleiac","$opts{out}.sampleiac1");
	&pathogenstatmerge("$opts{out}_sample.data.stat.xls","$opts{out}.plasmid","$opts{out}.plasmid1");
	system("$Bin/tabtk_xlsx $opts{out}.$reporttime.xlsx 原始reads数据:$opts{out}.raw1 barcode漂移结果:$opts{out}.barcode.pollute1 水过滤:$opts{out}.h2o1 引物唯一:$opts{out}.uniq1 PCR定量:$opts{out}.IAC1 样本定量:$opts{out}.sampleiac1 errorbarcode:$opts{out}.errorbarcode 质粒结果:$opts{out}.plasmid1");
	system("$Bin/tabtk_xlsx $opts{out}.$reporttime.stat.xlsx 原始reads数据:$opts{out}.rawforclient 引物原始reads数据:$opts{out}.raw 引物唯一:$opts{out}.uniq errorbarcode:$opts{out}.errorbarcode 统计:$opts{out}_sample.data.stat.xls.forclient");
}


sub confidence_level{
	my ($in)=@_;
	my $out="unknow";
	if ($in==1)
	{
		$out="低（小于70%）";
	}
	elsif($in==2)
	{
		$out="中（高于70%）";
	}
	elsif($in==3)
	{
		$out="中高（大于90%）";
	}
	elsif($in==4)
	{
		$out="高（大于99%）";
	}
	return $out;
}
sub checkNumber{    
	return shift =~ /^[+\-]?([1-9]\d*|0)(\.\d+)?([eE][+\-]?([1-9]\d*|0)(\.\d+)?)?$/;
}

sub positive2reporttxt{
	my %opts;
	($opts{positive},$opts{sample},$opts{dxdir},$opts{restriction},$opts{explain},$opts{drug_annot},$opts{out},$opts{pathogentype})=@_;
	

	my %positive;
	my %dpositive;
	my %sample;
	my %backgroup;
	my %lowcopy;
	my %restriction;

############sample inf
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $reporttime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
	if (defined$opts{sample})
	{
		my $excel = Spreadsheet::XLSX -> new ($opts{sample});
		my $sheet = ${$excel -> {Worksheet}}[0];
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		$sheet -> {MaxCol} ||= $sheet -> {MinCol};
		my $mytime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
		my ($sampleID,$collectdata,$company,$room,$samplename,$type,$dataID,$nongdu,$tiji,$librarynongdu,$check,$reportcode,$charge);
		my ($bed,$bingli,$linchuang,$sex,$age,$number,$ID,$wbc,$lbxb,$zxl,$crp,$pct,$pyjg,$jdjg,$jjjg,$lczd,$daili,$zdgz,$sjys,$cyrq);
		$reportcode=10000;
		$charge||=10000;
		$tiji||=10000;
		$librarynongdu||=10000;
		$cyrq||=10000;
		foreach my $col ($sheet -> {MinCol} .. $sheet -> {MaxCol})
		{
			my $b=$sheet -> {Cells} [0] [$col];
			next unless ($b);
			my $a=$b -> {Val};
			$sampleID     = $col if $a =~/样本编号/;
			$collectdata  = $col if $a =~/收样日期/;
			$company      = $col if $a =~/送检单位/;
			$room         = $col if $a =~/送检科室/;
			$samplename   = $col if $a =~/患者姓名/;
			$type         = $col if $a =~/样本类型/;
			$dataID       = $col if $a =~/实验编号/;
			$nongdu       = $col if $a =~/DNA浓度/;
			$tiji         = $col if $a =~/样本体积/;
			$check        = $col if $a=~/检测内容/;
			$librarynongdu= $col if $a =~/文库浓度/;
			$reportcode   = $col if $a =~/检测项目/;
			$charge       = $col if $a =~/付费/;
			$bed          = $col if $a=~/床号/;
			$bingli       = $col if $a=~/病历号/;
			$linchuang    = $col if $a=~/医院检测/;
			$sex          = $col if $a=~/性别/;
			$age          = $col if $a=~/年龄/;
			$number       = $col if $a=~/电话/;
			$ID           = $col if $a=~/身份证/;
			$wbc          = $col if $a=~/WBC/;
			$lbxb         = $col if $a=~/淋巴细胞/;
			$zxl          = $col if $a=~/中性粒/;
			$crp          = $col if $a=~/CRP/;
			$pct          = $col if $a=~/PCT/;
			$daili        = $col if $a=~/代理商/;
			$zdgz         = $col if $a=~/重点关注/;
			$sjys         = $col if $a=~/送检医生/;
			$cyrq         = $col if $a=~/采样日期/;
		}
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow})
		{
			my $sampleID_inf      ||= $sheet -> {Cells} [$row] [$sampleID];
			my $collectdata_inf   ||= $sheet -> {Cells} [$row] [$collectdata]; 
			my $company_inf       ||= $sheet -> {Cells} [$row] [$company];
			my $room_inf          ||= $sheet -> {Cells} [$row] [$room];
			my $samplename_inf    ||= $sheet -> {Cells} [$row] [$samplename];
			my $type_inf          ||= $sheet -> {Cells} [$row] [$type];
			my $dataID_inf        ||= $sheet -> {Cells} [$row] [$dataID];
			my $nongdu_inf        ||= $sheet -> {Cells} [$row] [$nongdu];
			my $tiji_inf          ||= $sheet -> {Cells} [$row] [$tiji];
			my $librarynongdu_inf ||= $sheet -> {Cells} [$row] [$librarynongdu];
			my $check_inf         ||= $sheet -> {Cells} [$row] [$check];
			my $reportcode_inf    ||= $sheet -> {Cells} [$row] [$reportcode];
			my $charge_inf        ||= $sheet -> {Cells} [$row] [$charge];
			my $bed_inf           ||= $sheet -> {Cells} [$row] [$bed];
			my $bingli_inf        ||= $sheet -> {Cells} [$row] [$bingli];
			my $linchuang_inf     ||= $sheet -> {Cells} [$row] [$linchuang];
			my $sex_inf           ||= $sheet -> {Cells} [$row] [$sex];
			my $age_inf           ||= $sheet -> {Cells} [$row] [$age];
			my $number_inf        ||= $sheet -> {Cells} [$row] [$number];
			my $ID_inf            ||= $sheet -> {Cells} [$row] [$ID];
			my $wbc_inf           ||= $sheet -> {Cells} [$row] [$wbc];
			my $lbxb_inf          ||= $sheet -> {Cells} [$row] [$lbxb];
			my $zxl_inf           ||= $sheet -> {Cells} [$row] [$zxl];
			my $crp_inf           ||= $sheet -> {Cells} [$row] [$crp];
			my $pct_inf           ||= $sheet -> {Cells} [$row] [$pct];
			my $daili_inf         ||= $sheet -> {Cells} [$row] [$daili];
			my $zdgz_inf          ||= $sheet -> {Cells} [$row] [$zdgz];
			my $sjys_inf          ||= $sheet -> {Cells} [$row] [$sjys];
			my $cyrq_inf          ||= $sheet -> {Cells} [$row] [$cyrq];
			################
			next unless ($sampleID_inf);
			my $name=$sampleID_inf -> {Val};
			######
			$sample{$name}{name} ||= $samplename_inf -> {Val};
			$sample{$name}{sex} ||= $sex_inf -> {Val};
			$sample{$name}{age} ||= $age_inf -> {Val};
			$sample{$name}{collect} ||= $collectdata_inf ->{Val};
			$sample{$name}{company} ||= $company_inf -> {Val};
			$sample{$name}{room} ||= $room_inf -> {Val};
			$sample{$name}{type} ||= $type_inf -> {Val};
			$sample{$name}{ID} ||= $sampleID_inf -> {Val};
			$sample{$name}{doctor}  ||= $bed_inf -> {Val};
			$sample{$name}{rep}   = $mytime;
			$sample{$name}{reportcode} ||= $reportcode_inf -> {Val};
			$sample{$name}{charge} ||= $charge_inf -> {Val};
			$sample{$name}{bingli} ||= $bingli_inf -> {Val};
			$sample{$name}{linchuang} ||= $linchuang_inf -> {Val};
			$sample{$name}{number} ||= $number_inf -> {Val};
			$sample{$name}{sfzID} ||= $ID_inf -> {Val};
			$sample{$name}{wbc} ||= $wbc_inf -> {Val};
			$sample{$name}{lbxb} ||= $lbxb_inf -> {Val};
			$sample{$name}{zxl} ||= $zxl_inf -> {Val};
			$sample{$name}{crp} ||= $crp_inf -> {Val};
			$sample{$name}{pct} ||= $pct_inf -> {Val};
			$sample{$name}{tiji} ||= $tiji_inf -> {Val};
			$sample{$name}{daili} ||= $daili_inf -> {Val};
			$sample{$name}{zdgz} ||= $zdgz_inf -> {Val};
			$sample{$name}{sjys} ||= $sjys_inf -> {Val};
			$sample{$name}{cyrq} ||= $cyrq_inf -> {Val};
#		$sample{$name}{name} =~s/\s//g;
#		$sample{$name}{name} =~s/[\(\)（）]//g;
		}
	}

	for my $samplename(sort keys %sample)
	{
		for (sort keys %{$sample{$samplename}})
		{
			$sample{$samplename}{$_}||="";
			$sample{$samplename}{$_}=~s/\s+//g;
		}
	}

###########################################################

#	open IN,"<$opts{restriction}" or die $!;
	open IN,"$Bin/configd.dist -z $opts{restriction}|" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
		$restriction{$inf[0]}{naiyao}=$inf[1];
		$restriction{$inf[0]}{chuanran}=$inf[2];
		$restriction{$inf[0]}{lowcopy}=$inf[3];
		$restriction{$inf[0]}{backgroup}=$inf[4];
		$restriction{$inf[0]}{latin}=$inf[5];
	}
	close IN;
	my %DXcode;
	my %alldx;
	for my $sample(sort keys %sample)
	{
		next if $sample{$sample}{reportcode}=~/检测项目/;
		open IN,"$Bin/.configx.dist -z $opts{dxdir}/$sample{$sample}{reportcode}.result|";
		while (<IN>)
		{
			chomp;
			my @inf=split;
			$DXcode{$sample{$sample}{reportcode}}{$inf[0]}=$inf[1];
			my $type;
			if ($inf[1] eq "D")
			{
				$type="DNA病毒";
			}
			elsif($inf[1] eq "drug")
			{
				$type="耐药基因";
			}
			elsif($inf[1] eq "F")
			{
				$type="真菌";
			}
			elsif($inf[1] eq "G")
			{
				$type="菌属";
			}
			elsif($inf[1] eq "H")
			{
				$type="螺旋体";
				$type="其他";
			}
			elsif($inf[1] eq "J")
			{
				$type="寄生虫";
				$type="其他";
			}
			elsif($inf[1] eq "N")
			{
				$type="革兰氏阴性菌";
			}
			elsif($inf[1] eq "O")
			{
				$type="其他";
			}
			elsif($inf[1] eq "P")
			{
				$type="革兰氏阳性菌";
			}
			elsif($inf[1] eq "R")
			{
				$type="RNA病毒";
			}
			elsif($inf[1] eq "ALL")
			{
				$type="内参";
			}
			$alldx{$inf[0]}=$type;
		}
		close IN;
	}

	open IN,"$Bin/configd.dist -z $opts{pathogentype}|" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split;
		my $type;
		if ($inf[1] eq "GramP")
		{
			$type="革兰氏阳性菌";
		}
		elsif($inf[1] eq "GramN")
		{
			$type="革兰氏阴性菌";
		}
		elsif($inf[1] eq "fungi")
		{
			$type="真菌";
		}
		elsif($inf[1] eq "Genus")
		{
			$type="菌属";
		}
		elsif($inf[1] eq "other")
		{
			$type="其他";
		}
		elsif($inf[2] eq "tmpDNAvirus")
		{
			$type="DNA病毒";
		}
		elsif($inf[2] eq "tmpRNAvirus")
		{
			$type="RNA病毒";
		}
		elsif($inf[2] eq "tmpdrug")
		{
			$type="耐药基因";
		}
		$alldx{$inf[0]}=$type unless defined$alldx{$inf[0]};
	}
	close IN;


	my %explain;
	open IN,"$Bin/configd.dist -z $opts{explain}|" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
#	my $name=$inf[0];
#	shift @inf;
#	my $out=join " ",@inf;
		$explain{$inf[0]}=$inf[1];
	}
	close IN;

	my %drugannot;
	open IN,"$Bin/configd.dist -z $opts{drug_annot}|" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
		$drugannot{$inf[1]}{annot}=$inf[2];
		$drugannot{$inf[1]}{target}=$inf[0];
	}
	close IN;

	my %druggene;
############主报告阳性菌确认耐药基因
	open IN,"<$opts{positive}" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split;next if /^样本/;
		$restriction{$inf[1]}{backgroup}||="-";
		$sample{$inf[0]}{type}||="肺泡灌洗液";
		if ($inf[4]=~/高/ || $inf[4]=~/中/)
		{
			if ($restriction{$inf[1]}{backgroup} ne "-" &&  $restriction{$inf[1]}{backgroup}=~/$sample{$inf[0]}{type}/)##背景菌
			{
				if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})
				{
					for my $gene(sort keys %drugannot)
					{
						if ($drugannot{$gene}{target}=~/$inf[1]/ || $drugannot{$gene}{target} eq "ALL")
						{
							$druggene{$inf[0]}{$gene}=1;
						}
					}
				}
			}
			else
			{
				for my $gene(sort keys %drugannot)
				{
					if ($drugannot{$gene}{target}=~/$inf[1]/ || $drugannot{$gene}{target} eq "ALL")
					{
						$druggene{$inf[0]}{$gene}=1;
					}
				}
			}
		}
	}
	close IN;
	#
	#
##################################################

	open IN,"<$opts{positive}" or die $!;
	open OUT,">$opts{out}" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split;
		next if /^样本/;
		next if /内参/;
		$explain{$inf[1]}||="-";
		$restriction{$inf[1]}{latin}||="-";
		if ($inf[4]=~/高/ || $inf[4]=~/中/)
		{
			if ($restriction{$inf[1]}{backgroup} ne "-" &&  $restriction{$inf[1]}{backgroup}=~/$sample{$inf[0]}{type}/)##背景菌
			{
				if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})##是否在DX内
				{
					print OUT "$inf[0]\t背景菌\t$alldx{$inf[1]}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$explain{$inf[1]}\t$sample{$inf[0]}{name}\n";
				}
				else
				{
					print OUT "$inf[0]\t不展示(背景非DX)\t$alldx{$inf[1]}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$explain{$inf[1]}\t$sample{$inf[0]}{name}\n";

				}
			}
			else
			{
				my $explain;
				if (defined$druggene{$inf[0]}{$inf[1]})
				{
					$explain=$drugannot{$inf[1]}{annot};
				}
				elsif(defined$explain{$inf[1]})
				{
					$explain=$explain{$inf[1]};
				}
				else
				{
					$explain="-";
				}

				if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})
				{
					print OUT "$inf[0]\t主报告\t$alldx{$inf[1]}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$explain\t$sample{$inf[0]}{name}\n";
#					print "$inf[1]\t$sample{$inf[0]}{reportcode}\n" unless defined$alldx{$inf[1]};
				}
				else
				{
					next if $alldx{$inf[1]} eq "耐药基因";
					print OUT "$inf[0]\t主报告*\t$alldx{$inf[1]}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$explain\t$sample{$inf[0]}{name}\n";
				}
			}
		}
		elsif($inf[4]=~/低/)
		{
			if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})
			{
				if (defined$restriction{$inf[1]}{lowcopy} && $restriction{$inf[1]}{lowcopy} ne "no")
				{
					next if $alldx{$inf[1]} eq  "耐药基因";
					next if($restriction{$inf[1]}{backgroup}=~/$sample{$inf[0]}{type}/);
					print OUT "$inf[0]\t灰区\t$alldx{$inf[1]}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$explain{$inf[1]}\t$sample{$inf[0]}{name}\n";
				}
			}
		}
	}
	close IN;
	close OUT;
}


sub sample_positive_stat{
	my %opts;
	($opts{positive},$opts{sample},$opts{out})=@_;


	my %sample;
#	$opts{out}=abs_path($opts{out});
############sample inf
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $reporttime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
	if (defined$opts{sample})
	{
		my $excel = Spreadsheet::XLSX -> new ($opts{sample});
		my $sheet = ${$excel -> {Worksheet}}[0];
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		$sheet -> {MaxCol} ||= $sheet -> {MinCol};
		my $mytime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
		my ($sampleID,$collectdata,$company,$room,$samplename,$type,$dataID,$nongdu,$tiji,$librarynongdu,$check,$reportcode,$charge);
		my ($bed,$bingli,$linchuang,$sex,$age,$number,$ID,$wbc,$lbxb,$zxl,$crp,$pct,$pyjg,$jdjg,$jjjg,$lczd,$daili);
		$reportcode=10000;
		$charge||=10000;
		foreach my $col ($sheet -> {MinCol} .. $sheet -> {MaxCol})
		{
			my $b=$sheet -> {Cells} [0] [$col];
			next unless ($b);
			my $a=$b -> {Val};
			$sampleID     = $col if $a =~/样本编号/;
			$collectdata  = $col if $a =~/收样日期/;
			$company      = $col if $a =~/送检单位/;
			$room         = $col if $a =~/送检科室/;
			$samplename   = $col if $a =~/患者姓名/;
			$type         = $col if $a =~/样本类型/;
			$dataID       = $col if $a =~/实验编号/;
			$nongdu       = $col if $a =~/DNA浓度/;
			$tiji         = $col if $a =~/DNA体积/;
			$check        = $col if $a=~/检测内容/;
			$librarynongdu= $col if $a =~/文库浓度/;
			$reportcode   = $col if $a =~/检测项目/;
			$charge       = $col if $a =~/付费/;
			$bed          = $col if $a=~/床号/;
			$bingli       = $col if $a=~/病历号/;
			$linchuang    = $col if $a=~/医院检测/;
			$sex          = $col if $a=~/性别/;
			$age          = $col if $a=~/年龄/;
			$number       = $col if $a=~/电话/;
			$ID           = $col if $a=~/身份证/;
			$wbc          = $col if $a=~/WBC/;
			$lbxb         = $col if $a=~/淋巴细胞/;
			$zxl          = $col if $a=~/中性粒/;
			$crp          = $col if $a=~/CRP/;
			$pct          = $col if $a=~/PCT/;
			$daili        = $col if $a=~/代理商/;
		}
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow})
		{
			my $sampleID_inf      ||= $sheet -> {Cells} [$row] [$sampleID];
			my $collectdata_inf   ||= $sheet -> {Cells} [$row] [$collectdata]; 
			my $company_inf       ||= $sheet -> {Cells} [$row] [$company];
			my $room_inf          ||= $sheet -> {Cells} [$row] [$room];
			my $samplename_inf    ||= $sheet -> {Cells} [$row] [$samplename];
			my $type_inf          ||= $sheet -> {Cells} [$row] [$type];
			my $dataID_inf        ||= $sheet -> {Cells} [$row] [$dataID];
			my $nongdu_inf        ||= $sheet -> {Cells} [$row] [$nongdu];
			my $tiji_inf          ||= $sheet -> {Cells} [$row] [$tiji];
			my $librarynongdu_inf ||= $sheet -> {Cells} [$row] [$librarynongdu];
			my $check_inf         ||= $sheet -> {Cells} [$row] [$check];
			my $reportcode_inf    ||= $sheet -> {Cells} [$row] [$reportcode];
			my $charge_inf        ||= $sheet -> {Cells} [$row] [$charge];
			my $bed_inf           ||= $sheet -> {Cells} [$row] [$bed];
			my $bingli_inf        ||= $sheet -> {Cells} [$row] [$bingli];
			my $linchuang_inf     ||= $sheet -> {Cells} [$row] [$linchuang];
			my $sex_inf           ||= $sheet -> {Cells} [$row] [$sex];
			my $age_inf           ||= $sheet -> {Cells} [$row] [$age];
			my $number_inf        ||= $sheet -> {Cells} [$row] [$number];
			my $ID_inf            ||= $sheet -> {Cells} [$row] [$ID];
			my $wbc_inf           ||= $sheet -> {Cells} [$row] [$wbc];
			my $lbxb_inf          ||= $sheet -> {Cells} [$row] [$lbxb];
			my $zxl_inf           ||= $sheet -> {Cells} [$row] [$zxl];
			my $crp_inf           ||= $sheet -> {Cells} [$row] [$crp];
			my $pct_inf           ||= $sheet -> {Cells} [$row] [$pct];
			my $daili_inf           ||= $sheet -> {Cells} [$row] [$daili];
			################
			next unless ($sampleID_inf);
			my $name=$sampleID_inf -> {Val};
			######
			$sample{$name}{name} ||= $samplename_inf -> {Val};
			$sample{$name}{sex} ||= $sex_inf -> {Val};
			$sample{$name}{age} ||= $age_inf -> {Val};
			$sample{$name}{collect} ||= $collectdata_inf ->{Val};
			$sample{$name}{company} ||= $company_inf -> {Val};
			$sample{$name}{room} ||= $room_inf -> {Val};
			$sample{$name}{type} ||= $type_inf -> {Val};
			$sample{$name}{ID} ||= $sampleID_inf -> {Val};
			$sample{$name}{doctor}  ||= $bed_inf -> {Val};
			$sample{$name}{rep}   = $mytime;
			$sample{$name}{reportcode} ||= $reportcode_inf -> {Val};
			$sample{$name}{charge} ||= $charge_inf -> {Val};
			$sample{$name}{bingli} ||= $bingli_inf -> {Val};
			$sample{$name}{linchuang} ||= $linchuang_inf -> {Val};
			$sample{$name}{number} ||= $number_inf -> {Val};
			$sample{$name}{sfzID} ||= $ID_inf -> {Val};
			$sample{$name}{wbc} ||= $wbc_inf -> {Val};
			$sample{$name}{lbxb} ||= $lbxb_inf -> {Val};
			$sample{$name}{zxl} ||= $zxl_inf -> {Val};
			$sample{$name}{crp} ||= $crp_inf -> {Val};
			$sample{$name}{pct} ||= $pct_inf -> {Val};
			$sample{$name}{tiji} ||= $tiji_inf -> {Val};
			$sample{$name}{daili} ||= $daili_inf -> {Val};
#		$sample{$name}{name} =~s/\s//g;
#		$sample{$name}{name} =~s/[\(\)（）]//g;
		}
	}

	my %report;
############target type  读取属性。细菌、真菌、DNA病毒、耐药,RNA病毒,其他病原体
	if (defined$opts{positive} && -s "$opts{positive}")
	{
		open IN,"<$opts{positive}" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split;
			next if $inf[4]=~/低|无/;
			my $zhixin=(split /（/,$inf[4])[0];
			push @{$report{$inf[0]}},"$inf[1]\t$zhixin($inf[2],$inf[3],$inf[5])";
		}
		close IN;
	}
######################


	my $count=0;
	for (sort keys %report)
	{
		my $tmp=@{$report{$_}};
		$count=$tmp if $tmp> $count;
	}

######################
	open OUT,">$opts{out}.tmp.file" or die $!;
	print OUT "分析日期\t样本编号\t样本名称\t样本类型\t检测项目\t代理商";
	print OUT "\t检出菌\t置信度(IAC,Reads,相对含量)" x $count;
	print OUT "\n";
	for (sort keys %sample)
	{
		my $pos="";
		if (defined$report{$_})
		{
			$pos=join "\t",@{$report{$_}};
		}
		next if $_=~/样本编号/;
		$sample{$_}{name}||="";
		$sample{$_}{type}||="";
		$sample{$_}{reportcode}||="";
		$sample{$_}{daili}||="";
		print OUT "$reporttime\t$_\t$sample{$_}{name}\t$sample{$_}{type}\t$sample{$_}{reportcode}\t$sample{$_}{daili}\t$pos\n";
	}
	close OUT;
	system("$Bin/tabtk_xlsx $opts{out} sheet1:$opts{out}.tmp.file;");
}

sub pathogenstatmerge{
	my %opts;
	($opts{inf},$opts{in},$opts{out})=@_;
	
	
	my @sort;
	my %hash;
	open IN,"<$opts{inf}" or die $!;
	my $head=<IN>;
	chomp $head;
	my @head=split /\t/,$head;
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
		push @sort,$inf[0];
		for (1..$#inf)
		{
			$hash{$inf[0]}{$head[$_]}=$inf[$_];
		}
	}
	close IN;

	open IN,"<$opts{in}" or die $!;
	open OUT,">$opts{out}" or die $!;
	$head=<IN>;
	chomp $head;
	@head=split /\t/,$head;
	print OUT "$head\n";
	for my $sort(@sort)
	{
		print OUT "$sort";
		for (1..$#head)
		{
			$hash{$sort}{$head[$_]}||="-";
			print OUT "\t$hash{$sort}{$head[$_]}";
		}
		print OUT "\n";
	}
	while (<IN>)
	{
		chomp;
		print OUT "$_\n";
	}
	close OUT;
	close IN;
}


sub quality_stat{
	my ($fastp,$primerlog,$alignlog,$out)=@_;
	my ($raw,$clean,$adapter,$dimmer,$pass,$single);
	open IN,"<$fastp" or die $!;
	while (<IN>)
	{
		chomp;
		if ($.==4 && /total_reads":(\d+),/)
		{
			$raw=int($1/2);
		}
		elsif ($.>12 && (/"q30_rate":(0\.\d+),/ || /"q30_rate":(0),/))
		{
			$clean=int($raw * $1);
		}
		elsif ( /adapter_trimmed_reads": (\d+),/)
		{
			$adapter=int($1/2);
		}
		elsif (/"too_short_reads": (\d+)/)
		{
			$dimmer=int($1/2);
		}
	}
	close IN;

	my $nspe=0;
	my $umap=0;
	if (defined$primerlog && -e "$primerlog")
	{
		open IN,"<$primerlog" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split /\t/,$_;
			next unless defined$inf[0];
			if ($inf[0] =~/Pair end mapped reads/)
			{
				$pass=$inf[1];
			}
			elsif($inf[0]=~/Primer dimer/)
			{
#				$dimmer+=$inf[1];
			}
			elsif($inf[0]=~/Nonspecific amplification/)
			{
				$nspe+=$inf[1];
			}
			elsif($inf[0]=~/Unmap reads/)
			{
				$umap+=$inf[1];
			}
			elsif($inf[0]=~/Single end mapped reads/)
			{
				$single+=$inf[1];
			}
		}
		close IN;
	}
	else
	{
		$pass=$clean;
		$nspe=0;
		$umap=0;
		$single=0;
	}

	my ($totalalign,$realalign)=(0,0);
	if (defined$alignlog)
	{
		$totalalign=`awk 'NR==2' $alignlog|cut -f2`;
		chomp $totalalign;
		$totalalign=$totalalign/2;
		$realalign=`awk 'NR==3' $alignlog|cut -f2`;
		chomp $realalign;
		$realalign=$realalign/2;
	}

	open OUT,">$out" or die $!;
	$raw||=1;
	$adapter||=0;

	my ($raw_p,$clean_p,$dimmer_p,$pass_p,$adapter_p,$nspe_p,$unmap_p,$single_p);
	$clean_p=sprintf "%.2f",$clean * 100 /$raw;
	$adapter_p=sprintf "%.2f",$adapter * 100 /$raw;
	$dimmer_p=sprintf "%.2f",$dimmer * 100 /$raw;
	$pass_p=sprintf "%.2f",$pass * 100 /$raw;
	$nspe_p=sprintf "%.2f",$nspe * 100 /$raw;
	$unmap_p=sprintf "%.2f",$umap * 100 /$raw;
	my $totalalign_p=sprintf "%.2f",$totalalign * 100 /$raw;
	my $realalign_p=sprintf "%.2f",$realalign * 100 /$raw;
	$single_p=sprintf "%.2f",$single * 100 /$raw;


	print OUT "Total pair reads\t$raw(100%)\n";
	print OUT "high quality pair reads(Q30)\t$clean($clean_p%)\n";
	print OUT "Adapter filter\t$adapter($adapter_p%)\n";
	print OUT "Dimmer\t$dimmer($dimmer_p%)\n";
	print OUT "Nonspecific amplification\t$nspe($nspe_p%)\n";
	print OUT "Single end mapped reads\t$single($single_p%)\n";
	print OUT "Unknow Primer\t$umap($unmap_p%)\n";
	print OUT "Clean pair reads\t$pass($pass_p%)\n";
	print OUT "Total Align\t$totalalign($totalalign_p%)\n" if defined$alignlog;
	print OUT "Real Align\t$realalign($realalign_p%)\n" if defined$alignlog;
	close OUT;
}


sub add_desc{
	my ($desc,$line,$input,$NA,$output)=@_;
	$NA ||= "0";
	open IN,"<$desc" or die $!;
	my %hash;
	my $n=0;
	$line=$line-1;
	my $sep||='\t';
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
		my $b=$inf[0];
		shift @inf;
		my $a=join "\t",@inf;
		$n=@inf;
		$hash{$b}=$a;
	}
	close IN;
	open IN,"<$input";
	open OUT,">$output";
	while (<IN>)
	{
		chomp;
		my @inf=split/$sep/,$_;
		print OUT "$_\t$hash{$inf[$line]}\n" if defined$hash{$inf[$line]};
		print OUT "$_"."\t$NA" x $n ."\n" unless defined$hash{$inf[$line]};
	}
	close OUT;
	close IN;
}

sub merge_sampleinf_datainf{
	my ($sample,$data,$out)=@_;
	my %inf;
	open IN,"<$sample" or die $!;
	my $head=<IN>;
	chomp $head;
	my @head=split /\t/,$head;

	while (<IN>)
	{
		chomp;
		my @inf=split/\t/,$_;
		for my $col(1..$#inf)
		{
			$inf{$head[$col]}{$inf[0]}=$inf[$col];
		}
	}
	close IN;


	my @sort;
	open OUT,">$out" or die $!;
	open IN,"<$data" or die $!;
	while (<IN>)
	{
		chomp;
		if($.==1)
		{
			@sort=split /\t/,$_;
		}
		print OUT "$_\n";
	}
	close IN;

	for my $information(sort keys %inf)
	{
		print OUT "$information";
		for my $col(1..$#sort)
		{
			print OUT "\t$inf{$information}{$sort[$col]}" if defined$inf{$information}{$sort[$col]};
			print OUT "\t-" unless defined$inf{$information}{$sort[$col]};
		}
		print OUT "\n";
	}
	close OUT;
}


sub positive2reporttxtv11{
	my %opts;
	($opts{positive},$opts{sample},$opts{dxdir},$opts{restriction},$opts{out},$opts{stat})=@_;

	my %positive;
	my %dpositive;
	my %sample;
	my %backgroup;
	my %lowcopy;
	my %restriction;

############sample inf
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $reporttime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
	if (defined$opts{sample})
	{
		my $excel = Spreadsheet::XLSX -> new ($opts{sample});
		my $sheet = ${$excel -> {Worksheet}}[0];
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		$sheet -> {MaxCol} ||= $sheet -> {MinCol};
		my $mytime = sprintf("%04d%02d%02d", $year+1900,$mon+1,$mday);
		my ($sampleID,$collectdata,$company,$room,$samplename,$type,$dataID,$nongdu,$tiji,$librarynongdu,$check,$reportcode,$charge);
		my ($bed,$bingli,$linchuang,$sex,$age,$number,$ID,$wbc,$lbxb,$zxl,$crp,$pct,$pyjg,$jdjg,$jjjg,$lczd,$daili,$zdgz,$sjys,$cyrq);
		$reportcode=10000;
		$charge||=10000;
		$tiji||=10000;
		$librarynongdu||=10000;
		$cyrq||=10000;
		foreach my $col ($sheet -> {MinCol} .. $sheet -> {MaxCol})
		{
			my $b=$sheet -> {Cells} [0] [$col];
			next unless ($b);
			my $a=$b -> {Val};
			$sampleID     = $col if $a =~/样本编号/;
			$collectdata  = $col if $a =~/收样日期/;
			$company      = $col if $a =~/送检单位/;
			$room         = $col if $a =~/送检科室/;
			$samplename   = $col if $a =~/患者姓名/;
			$type         = $col if $a =~/样本类型/;
			$dataID       = $col if $a =~/实验编号/;
			$nongdu       = $col if $a =~/DNA浓度/;
			$tiji         = $col if $a =~/样本体积/;
			$check        = $col if $a=~/检测内容/;
			$librarynongdu= $col if $a =~/文库浓度/;
			$reportcode   = $col if $a =~/检测项目/;
			$charge       = $col if $a =~/付费/;
			$bed          = $col if $a=~/床号/;
			$bingli       = $col if $a=~/病历号/;
			$linchuang    = $col if $a=~/医院检测/;
			$sex          = $col if $a=~/性别/;
			$age          = $col if $a=~/年龄/;
			$number       = $col if $a=~/电话/;
			$ID           = $col if $a=~/身份证/;
			$wbc          = $col if $a=~/WBC/;
			$lbxb         = $col if $a=~/淋巴细胞/;
			$zxl          = $col if $a=~/中性粒/;
			$crp          = $col if $a=~/CRP/;
			$pct          = $col if $a=~/PCT/;
			$daili        = $col if $a=~/代理商/;
			$zdgz         = $col if $a=~/重点关注/;
			$sjys         = $col if $a=~/送检医生/;
			$cyrq         = $col if $a=~/采样日期/;
		}
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow})
		{
			my $sampleID_inf      ||= $sheet -> {Cells} [$row] [$sampleID];
			my $collectdata_inf   ||= $sheet -> {Cells} [$row] [$collectdata]; 
			my $company_inf       ||= $sheet -> {Cells} [$row] [$company];
			my $room_inf          ||= $sheet -> {Cells} [$row] [$room];
			my $samplename_inf    ||= $sheet -> {Cells} [$row] [$samplename];
			my $type_inf          ||= $sheet -> {Cells} [$row] [$type];
			my $dataID_inf        ||= $sheet -> {Cells} [$row] [$dataID];
			my $nongdu_inf        ||= $sheet -> {Cells} [$row] [$nongdu];
			my $tiji_inf          ||= $sheet -> {Cells} [$row] [$tiji];
			my $librarynongdu_inf ||= $sheet -> {Cells} [$row] [$librarynongdu];
			my $check_inf         ||= $sheet -> {Cells} [$row] [$check];
			my $reportcode_inf    ||= $sheet -> {Cells} [$row] [$reportcode];
			my $charge_inf        ||= $sheet -> {Cells} [$row] [$charge];
			my $bed_inf           ||= $sheet -> {Cells} [$row] [$bed];
			my $bingli_inf        ||= $sheet -> {Cells} [$row] [$bingli];
			my $linchuang_inf     ||= $sheet -> {Cells} [$row] [$linchuang];
			my $sex_inf           ||= $sheet -> {Cells} [$row] [$sex];
			my $age_inf           ||= $sheet -> {Cells} [$row] [$age];
			my $number_inf        ||= $sheet -> {Cells} [$row] [$number];
			my $ID_inf            ||= $sheet -> {Cells} [$row] [$ID];
			my $wbc_inf           ||= $sheet -> {Cells} [$row] [$wbc];
			my $lbxb_inf          ||= $sheet -> {Cells} [$row] [$lbxb];
			my $zxl_inf           ||= $sheet -> {Cells} [$row] [$zxl];
			my $crp_inf           ||= $sheet -> {Cells} [$row] [$crp];
			my $pct_inf           ||= $sheet -> {Cells} [$row] [$pct];
			my $daili_inf         ||= $sheet -> {Cells} [$row] [$daili];
			my $zdgz_inf          ||= $sheet -> {Cells} [$row] [$zdgz];
			my $sjys_inf          ||= $sheet -> {Cells} [$row] [$sjys];
			my $cyrq_inf          ||= $sheet -> {Cells} [$row] [$cyrq];
			################
			next unless ($sampleID_inf);
			my $name=$sampleID_inf -> {Val};
			######
			if($samplename_inf){
				$sample{$name}{name} = $samplename_inf -> {Val};
			}else{
				$sample{$name}{name}="-";
			}
			if($sex_inf){
				$sample{$name}{sex} = $sex_inf -> {Val};
			}else{
				$sample{$name}{sex}="-";
			}
			if($age_inf){
				$sample{$name}{age} = $age_inf -> {Val};
			}else{
				$sample{$name}{age} ="-";
			}
			if($collectdata_inf){
				$sample{$name}{collect} = $collectdata_inf ->{Val};
			}else{
				$sample{$name}{collect} ="";
			}
			if($company_inf){
				$sample{$name}{company} ||= $company_inf -> {Val};
			}else{
				$sample{$name}{company} ="-";
			}
			if($room_inf){
				$sample{$name}{room} = $room_inf -> {Val};
			}else{
				$sample{$name}{room} ="-";
			}
			if($type_inf){
				$sample{$name}{type} = $type_inf -> {Val};
			}
			else{
				$sample{$name}{type} ="-";
			}
			$sample{$name}{ID} ||= $sampleID_inf -> {Val};
			$sample{$name}{doctor}  ||= $bed_inf -> {Val};
			$sample{$name}{rep}   = $mytime;
			$sample{$name}{reportcode} ||= $reportcode_inf -> {Val};
			$sample{$name}{charge} ||= $charge_inf -> {Val};
			$sample{$name}{bingli} ||= $bingli_inf -> {Val};
			if($linchuang_inf){
				$sample{$name}{linchuang} = $linchuang_inf -> {Val};
			}
			else{
				$sample{$name}{linchuang} ="-";
			}
			$sample{$name}{number} ||= $number_inf -> {Val};
			$sample{$name}{sfzID} ||= $ID_inf -> {Val};
			$sample{$name}{wbc} ||= $wbc_inf -> {Val};
			$sample{$name}{lbxb} ||= $lbxb_inf -> {Val};
			$sample{$name}{zxl} ||= $zxl_inf -> {Val};
			$sample{$name}{crp} ||= $crp_inf -> {Val};
			$sample{$name}{pct} ||= $pct_inf -> {Val};
			$sample{$name}{tiji} ||= $tiji_inf -> {Val};
			$sample{$name}{daili} ||= $daili_inf -> {Val};
			$sample{$name}{zdgz} ||= $zdgz_inf -> {Val};
			$sample{$name}{sjys} ||= $sjys_inf -> {Val};
			$sample{$name}{cyrq} ||= $cyrq_inf -> {Val};
			if($nongdu_inf){
				$sample{$name}{nongdu} ||= $nongdu_inf -> {Val};
			}
			else{
				$sample{$name}{nongdu} ="-";
			}
			if($librarynongdu_inf){
				$sample{$name}{librarynongdu} ||= $librarynongdu_inf -> {Val};
			}
			else{
				$sample{$name}{librarynongdu} ="-";
			}
#		$sample{$name}{name} =~s/\s//g;
#		$sample{$name}{name} =~s/[\(\)（）]//g;
		}
	}


	my $total_sample=keys %sample;
	for my $samplename(sort keys %sample)
	{
		for (sort keys %{$sample{$samplename}})
		{
			$sample{$samplename}{$_}||="-";
			$sample{$samplename}{$_}=~s/\s+//g;
		}
	}



###########################################################

	open IN,"$Bin/configd.dist -z $opts{restriction}|" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split /\t/,$_;
		$restriction{$inf[0]}{naiyao}=$inf[1];
		$restriction{$inf[0]}{chuanran}=$inf[2];
		$restriction{$inf[0]}{lowcopy}=$inf[3];
		$restriction{$inf[0]}{backgroup}=$inf[4];
		$restriction{$inf[0]}{latin}=$inf[5];
		$restriction{$inf[0]}{fenlei1}=$inf[6];
		$restriction{$inf[0]}{fenlei2}=$inf[7];
		$restriction{$inf[0]}{explain}=$inf[8];


	}
	close IN;

	my (@head,@shiji,@dengji,@real,%hashzhikong,%hashalinread,%hashqua,@quality);
	open IN,"<$opts{stat}";
	while(<IN>)
	{
		chomp;
		my $li=$_;
		my @ar=split(/\t/);
		if($ar[0] eq "Sample")
		{
			@head=split(/\t/,$li);
		}
		if($ar[0] eq "实际标签")
		{
			@shiji=split(/\t/,$li);
		}
		if($ar[0] eq "登记标签")
		{
			@dengji=split(/\t/,$li);
		}
		if($ar[0] eq "Real Align")
		{
			@real=split(/\t/,$li);
		}
		if($ar[0]=~/high/)
		{
			@quality=split(/\t/,$li);
		}

	}
	close IN;

	for(my $i=1;$i<@head;$i++)#增加DMB
	{
		$hashalinread{$head[$i]}=$real[$i];
		$hashqua{$head[$i]}=$quality[$i];
		$real[$i]=~/(\d+)\((\d+)/;
		if($1>10000 && $2>0.6)
		{
			$hashzhikong{$head[$i]}.="Y,";
		}
		else
		{
			$hashzhikong{$head[$i]}.="D,";
		}
		if($shiji[$i] eq "-")
		{
			$hashzhikong{$head[$i]}.="M,";
		}
		elsif($shiji[$i] ne $dengji[$i])
		{
			$hashzhikong{$head[$i]}.="B,";
		}
		else
		{
			$hashzhikong{$head[$i]}.="Y,";
		}
	}

	my %DXcode;
	for my $sample(sort keys %sample)
	{
		next if $sample{$sample}{reportcode}=~/检测项目/;
		unless(-e "$opts{dxdir}/$sample{$sample}{reportcode}.result")
		{
			next;
		}
		my $pass=0;
		open IN,"$Bin/.configx.dist -z $opts{dxdir}/$sample{$sample}{reportcode}.result|";
		while (<IN>)
		{
			chomp;
			my @inf=split;
			next unless($inf[0]);#有些有空行
			$DXcode{$sample{$sample}{reportcode}}{$inf[0]}++;
			if($pass==0)#判断该DX是否有耐药，用来判断报告是否需要加耐药框框
			{
				if($restriction{$inf[0]}{fenlei1} && $restriction{$inf[0]}{fenlei1}=~/毒力|耐药/)
				{
					$hashzhikong{$sample}.="耐药,";
					$pass=1;
				}
			}
		}
		close IN;
	}

	my %hashpathsta;
	my %druggene;
	my %hashrawreads;#用于生成原始数据表格
	my %hashchuanran_judge;#判断是否输出传染病表格
	my %hashchuanran;
############主报告阳性菌确认耐药基因
	open IN,"<$opts{positive}" or die $!;
	while (<IN>)
	{
		chomp;
		my @inf=split;next if /^样本/;
		$hashpathsta{$inf[1]}=0 unless(exists $hashpathsta{$inf[1]});
		#生成原始数据表格
		if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})#且在DX范围内
		{
			$hashrawreads{$inf[0]}.="$inf[0]\t$inf[1]\t$inf[3]\n";
			#	print "$inf[0]\t$inf[1]\t$inf[3]\n";

		}
		#
		next if($inf[1] eq "定量内参");
#	print "soga\t$inf[1]\t$restriction{$inf[1]}{chuanran}\n";
		unless(exists $restriction{$inf[1]}{chuanran})
		{
#		print STDERR "$inf[1]病原体注释缺失\n";#展示restrict表里面没有的病原体
			next;

		}
		if($restriction{$inf[1]}{chuanran} eq "yes")
		{
			if($inf[3]>0)
			{
				$hashchuanran_judge{$inf[0]}=1;
			}
			$hashchuanran{$inf[0]}.="$inf[0]\t$inf[1]\t$inf[3]\n";
		}
#	next if($inf[1] eq "定量内参");
		if ($inf[4]=~/高/ || $inf[4]=~/中/)
		{
			$hashpathsta{$inf[1]}++;
			if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})#且在DX范围内
			{
				unless($restriction{$inf[1]}{chuanran})
				{
					print "soga$inf[1] 病原体未发现\n";
				}
				unless($restriction{$inf[1]}{chuanran} eq "yes")
				{
					unless($restriction{$inf[1]}{fenlei1})
					{
						print "soga\t$inf[1]\n";
					}
					if($restriction{$inf[1]}{fenlei1}=~/革兰氏|属|其他/)#DX范围内，展示在报告中的，记录编号，细菌名称，用于判断后面是否报耐药或者毒力
					{
						$druggene{"$inf[0],$inf[1]"}++;
						$druggene{"$inf[0],ALL"}++;
					}
				}
			}
		}

	}
	close IN;

#输出原始数据表
	my $dir=dirname("$opts{out} ");
	system("mkdir -p $dir");
	unless(-e "$dir/singlesample"){`mkdir $dir/singlesample`;}
	foreach my $key(sort keys %hashrawreads)
	{
		open OUT,"> $dir/singlesample/$key.txt";
		print OUT "$hashrawreads{$key}";
		close OUT;
		#print "$hashrawreads{$key}";
	}
#输出传染病表格
	foreach my $key(sort keys %hashchuanran_judge)
	{
		unless(-e "$dir/infectious"){`mkdir $dir/infectious`;}
		open OUT,"> $dir/infectious/$key.txt";
		print OUT "$hashchuanran{$key}\n";
		close OUT;
	}
	#
	#
##################################################


	my %hashout;
	open IN,"<$opts{positive}" or die $!;
	open OUT,">$opts{out}" or die $!;

	print OUT "样本编号\t患者姓名\t患者年龄\t样本类型\t核酸浓度\t临床信息\t报告标签\t报告区域\t病原体分类\t病原体中文名\t病原体英文文名\treads数\t病原体拷贝数\t信号强度\t病原体注\t代理商名称\t送检单位名称\t性别\t科室\t文库浓度\t正确配对reads\tQ30过滤reads\t该病原同批检出数量\n";
	while (<IN>)
	{
		chomp;
		my @inf=split;
		if(length($inf[2])>6)
		{
			$inf[2]=substr($inf[2],0,6);#不超过10^6
		}
		next if /^样本/;
		$restriction{$inf[1]}{explain}||="-";
		$restriction{$inf[1]}{latin}||="-";
		#定量内参都展示一行，避免有些是0不展示
		if($inf[1]=~/定量内参/)
		{
			$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"不展示"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t不展示\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
		}
		else
		{
			if ($inf[4]=~/高/ || $inf[4]=~/中/)
			{
				next unless($restriction{$inf[1]}{fenlei1});#如果没有这个菌名的话，跳过
				if($restriction{$inf[1]}{fenlei1}=~/耐药|毒力/)#如果匹配到耐药毒力，单独处理，然后下一行
				{
					my @gene=split(/,/,$restriction{$inf[1]}{fenlei2});
					for(my $i=0;$i<@gene;$i++)
					{
						if(exists $druggene{"$inf[0],$gene[$i]"})
						{
							$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"主报告"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t主报告\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
							last;
						}
					}
					next;
				}
				if ($restriction{$inf[1]}{backgroup} ne "-" &&  $restriction{$inf[1]}{backgroup}=~/$sample{$inf[0]}{type}/)##背景菌
				{
					if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})##是否在DX内
					{
						$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"背景菌"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t背景菌\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
					}
					else
					{
						$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"不展示(背景非DX)"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t不展示(背景非DX)\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";

					}
				}
				else
				{
					if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})#未在背景菌，且在DX范围内
					{
						if($restriction{$inf[1]}{chuanran} eq "yes")
						{
							$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"不展示(DX内传染病)"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t不展示(DX内传染病)\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
						}
						else
						{
							$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"主报告"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t主报告\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
						}
					}
					else	#未在背景菌，且不在DX范围内
					{
						next if ($restriction{$inf[1]}{fenlei1}=~/耐药|毒力/);
						$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"不展示(非DX)"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t不展示(非DX)\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
					}
				}
			}
			elsif($inf[4]=~/低/)
			{
				if (defined$DXcode{$sample{$inf[0]}{reportcode}}{$inf[1]})
				{
					next if($restriction{$inf[1]}{backgroup}=~/$sample{$inf[0]}{type}/);#背景菌也不展示在灰区中
					if ( $restriction{$inf[1]}{lowcopy} ne "no")
					{
						next if ($restriction{$inf[1]}{fenlei1}=~/耐药|毒力/);#耐药和毒力都不展示
						$hashout{$sample{$inf[0]}{daili}}{$sample{$inf[0]}{name}}{$inf[0]}{"灰区"}{$inf[3]}.="$inf[0]\t$sample{$inf[0]}{name}\t$sample{$inf[0]}{age}\t$sample{$inf[0]}{type}\t$sample{$inf[0]}{nongdu}\t$sample{$inf[0]}{linchuang}\t$hashzhikong{$inf[0]}\t灰区\t$restriction{$inf[1]}{fenlei1}\t$inf[1]\t$restriction{$inf[1]}{latin}\t$inf[3]\t$inf[2]\t$inf[4]\t$restriction{$inf[1]}{explain}\t$sample{$inf[0]}{daili}\t$sample{$inf[0]}{company}\t$sample{$inf[0]}{sex}\t$sample{$inf[0]}{room}\t$sample{$inf[0]}{librarynongdu}\t$hashalinread{$inf[0]}\t$hashqua{$inf[0]}\t$hashpathsta{$inf[1]}\/$total_sample\n";
					}
				}
			}
		}	
	}
	close IN;


	foreach my $key1(sort keys %hashout)
	{
		foreach my $key2(sort keys %{$hashout{$key1}})
		{
			foreach my $key3(sort keys %{$hashout{$key1}{$key2}})
			{
				foreach my $key4(sort keys  %{$hashout{$key1}{$key2}{$key3}})
				{
					foreach my $key5(sort {$b<=>$a} keys %{$hashout{$key1}{$key2}{$key3}{$key4}})
					{
						print OUT $hashout{$key1}{$key2}{$key3}{$key4}{$key5};
					}
				}
			}
		}
	}

	close OUT;

#添加结核突变文件
	system("$Bin/tabtk_xlsx  $opts{out}.xlsx sheet1:$opts{out}");
}


sub make_mutation_stat{
	my %optsmutstat;
	($optsmutstat{indir},$optsmutstat{out})=@_;
	


	my %sample;
	my %alldrug;

	my @sampledir=`ls -d $optsmutstat{indir}/*`;
	for (@sampledir)
	{
		chomp;
		my $dirname=basename($_);
		next if $dirname=~/water|h2o|shui/i;
		$sample{$dirname}=0;
	}

	my @file=`ls $optsmutstat{indir}/*/align/*JH.drug_class.xls`;
	my %drug_stat;
	for my $file(@file)
	{
		chomp $file;
		open IN,"<$file" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split /\t/,$_;
			push @{$drug_stat{$inf[0]}{$inf[2]}{$inf[1]}},$inf[3];
			$alldrug{$inf[2]}{$inf[3]}=1;
		}
		close IN;
	}

	my $len=15;
#######
	open OUT,">$optsmutstat{out}"or die $!;
	@file=`ls $optsmutstat{indir}/*/align/*.JH.mutation.class`;
	for my $file(@file)
	{
		chomp $file;
		open IN,"<$file" or die $!;
		while (<IN>)
		{
			chomp;
			my @inf=split /\t/,$_;
			$len=$#inf;
			$sample{$inf[0]}=1;
			if ($inf[14] eq "可能" )
			{
				$inf[14]="-";
			}
			if ($inf[15] eq "-" && $inf[1] eq "hot_pos")
			{
				$inf[15]="该突变未有相关文献报道，但该位置发生的其他基因型突变有相关文献报道过耐药现象";
			}
			my ($onelinep,$onelinen,$twolinep,$twolinen,$ntmp,$ntmn)=("-","-","-","-","-","-");
			if (defined$drug_stat{$inf[0]}{'TB一线'}{'检出'})
			{
				$onelinep=join ",",@{$drug_stat{$inf[0]}{'TB一线'}{'检出'}};
			}
			if (defined$drug_stat{$inf[0]}{'TB一线'}{'未检出'})
			{
				$onelinen=join ",",@{$drug_stat{$inf[0]}{'TB一线'}{'未检出'}};
			}
			if (defined$drug_stat{$inf[0]}{'TB二线'}{'检出'})
			{
				$twolinep=join ",",@{$drug_stat{$inf[0]}{'TB二线'}{'检出'}};
			}
			if (defined$drug_stat{$inf[0]}{'TB二线'}{'未检出'})
			{
				$twolinen=join ",",@{$drug_stat{$inf[0]}{'TB二线'}{'未检出'}};
			}
			if (defined$drug_stat{$inf[0]}{'NTM药物'}{'检出'})
			{
				$ntmp=join ",",@{$drug_stat{$inf[0]}{'NTM药物'}{'检出'}};
			}
			if (defined$drug_stat{$inf[0]}{'NTM药物'}{'未检出'})
			{
				$ntmn=join ",",@{$drug_stat{$inf[0]}{'NTM药物'}{'未检出'}};
			}
			my $out=join "\t",@inf;
			print OUT "$out\t$onelinep\t$onelinen\t$twolinep\t$twolinen\t$ntmp\t$ntmn\n";
		}
	}


	my @oneline=sort keys %{$alldrug{'TB一线'}};
	my @twoline=sort keys %{$alldrug{'TB二线'}};
	my @ntm=sort keys %{$alldrug{'NTM药物'}};

	my $oneline=join ",",@oneline;
	my $twoline=join ",",@twoline;
	my $ntm=join ",",@ntm;
	for (sort keys %sample)
	{
		if ($sample{$_}==0)
		{
			print OUT "$_"."\t-" x $len ."\t-\t$oneline\t-\t$twoline\t-\t$ntm\n";
		}
	}

	close OUT;
}
