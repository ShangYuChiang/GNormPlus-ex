#!/usr/bin/env perl

#===================================================
# Author: Chih-Hsuan Wei
# Software: MuExt
# Description: extract article via E-utilities
#===================================================

use LWP::Simple;
use HTML::Entities;
sub OutputFormat
{
	my $Article= $_[0];
	my $output= $_[1];
	my %ent_hash=();
	
	my $timestamp;
	my ($sec,$min,$hour,$day,$mon,$year)=localtime(time);
	$year=$year+1900;
	$timestamp=$year;
	$mon=$mon+1;
	if($mon<10){$timestamp=$timestamp."0".$mon;}
	else{$timestamp=$timestamp.$mon;}
	if($day<10){$timestamp=$timestamp."0".$day;}
	else{$timestamp=$timestamp.$day;}
	if($hour<10){$timestamp=$timestamp."0".$hour;}
	else{$timestamp=$timestamp.$hour;}
	if($min<10){$timestamp=$timestamp."0".$min;}
	else{$timestamp=$timestamp.$min;}
	if($sec<10){$timestamp=$timestamp."0".$sec;}
	else{$timestamp=$timestamp.$sec;}
	
	my $identifier=$timestamp;
	
	$ent_hash{"&#x0005b;"}="[";
	$ent_hash{"&#x0005d;"}="]";
	$ent_hash{"&gt;"}=">";
	
	
	foreach $each_ent (keys %ent_hash)
	{
		$Article =~s/$each_ent/$ent_hash{$each_ent}/g;
	}
	$Article =~s/\&\#x[a-z0-9][a-z0-9]+;/ /g;
	
	#=====
	#PMID
	if($Article =~/<article\-id pub\-id\-type=\"pmid\">(.+?)<\/article\-id>/ || $Article =~/<PMID.+?>(.+?)<\/PMID>/)
	{
		$identifier=$1;
	}
	
	open output,">>".$output;

	#=====
	#TITLE
	if($Article =~/<title\-group>[ ]*<article\-title>(.+?)<\/article\-title>/  || $Article =~/<[Aa]rticleTitle.*?>(.+?)<\/[Aa]rticleTitle>/ )
	{
		my $sentence=$1;
		$sentence =~s/<[^\<\>]+?>/ /g;
		$sentence =~s/[ ]+$//g;$sentence =~s/^[ ]+//g;
		print output $identifier."\|t\|".$sentence."\n";
	}
	else
	{
		print output $identifier."\|t\|-No title-\n";
	}
	
	#=====
	#ABSTRACT
	if($Article =~/<[Aa]bstract.*?>(.*?)<\/[Aa]bstract>/)
	{
		my $sentence=$1;
		$sentence =~s/<[^\<\>]+?>/ /g;
		$sentence =~s/[ ]+$//g;$sentence =~s/^[ ]+//g;
		print output $identifier."\|a\|".$sentence."\n";
	}
	elsif($Article =~/<[Aa]bstractText.*?>(.+?)<\/[Aa]bstractText>/ )
	{
		my $sentence=$1;
		$sentence =~s/<[^\<\>]+?>/ /g;
		$sentence =~s/[ ]+$//g;$sentence =~s/^[ ]+//g;
		print output $identifier."\|a\|".$sentence."\n";
	}
	else
	{
		print output $identifier."\|a\|-No abstract-\n";
	}
	
	print output "\n";
	close output;
}
sub GetResultViaEutilities 
{
	my $type= $_[0];
	my $input= $_[1];
	my $output= $_[2];
	my $Article="";

	if($type eq "PMID")
	{
		my $efetch="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=".$input."&retmode=xml";
		$Article = get($efetch);
		$Article =~s/[\n\r\t]//g;
		decode_entities($Article);
		#$Article =~s/→/ > /g;
		#$Article =~s/↑/\(\+\)/g;
		#$Article =~s/↓/\(\-\)/g;
		#$Article =~s/←/ < /g;
		#$Article =~s/Δ/\/\\/g;
		#$Article =~s/α/a/g;
		#$Article =~s/β/�/g;
		#$Article =~s/×/�/g;
		#$Article =~s/κ/k/g;
		#$Article =~s/γ/y/g;
		#$Article =~s/�/b/g;
		#$Article =~s/�/\-/g;
		#$Article =~s/≥/>= /g;
		#$Article =~s/��/<=/g;
		#$Article =~s/ //g;
		#$Article =~s/�/ /g;
		#$Article =~s/ /   /g;
		#$Article =~s/��/  /g;
		#$Article =~s/�/ /g;
		#$Article =~s/�/ /g;
		#$Article =~s/μ/�/g;
		#$Article =~s/\&lt\;/</g;
		#$Article =~s/\&gt\;/>/g;	
		#$Article =~s/\&quot\;/"/g;
		#$Article =~s/\&[a-z]{1,4}\;/\&/g;
		#$Article =~s/\&[a-z]{1,4}\;/\&/g;
		#$Article =~s/\&/\_/g;
		##$Article =~s/\+/ /g;
		#$Article=~s/[^A-Za-z0-9~!@#^%&*()-_+={}:<>\[\]\;\"',.\/]/ /g;
		#$Article =~s/[ ]+/ /g;
		
		open output,">".$output;
		close output;
		OutputFormat($Article,$output);
	}
	elsif($type eq "PMCID")
	{
		my $efetch="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=".$input."";
		$Article = get($efetch);
		$Article =~s/[\n\r\t]//g;
		decode_entities($Article);
		#$Article =~s/→/ > /g;
		#$Article =~s/↑/\(\+\)/g;
		#$Article =~s/↓/\(\-\)/g;
		#$Article =~s/←/ < /g;
		#$Article =~s/Δ/\/\\/g;
		#$Article =~s/α/a/g;
		#$Article =~s/β/�/g;
		#$Article =~s/×/�/g;
		#$Article =~s/κ/k/g;
		#$Article =~s/γ/y/g;
		#$Article =~s/�/b/g;
		#$Article =~s/�/\-/g;
		#$Article =~s/≥/>= /g;
		#$Article =~s/��/<=/g;
		#$Article =~s/ //g;
		#$Article =~s/�/ /g;
		#$Article =~s/ /   /g;
		#$Article =~s/��/  /g;
		#$Article =~s/�/ /g;
		#$Article =~s/�/ /g;
		#$Article =~s/μ/�/g;
		#$Article =~s/\&lt\;/</g;
		#$Article =~s/\&gt\;/>/g;	
		#$Article =~s/\&quot\;/"/g;
		#Article =~s/\&[a-z]{1,4}\;/\&/g;
		#$Article =~s/\&[a-z]{1,4}\;/\&/g;
		#$Article =~s/\&/\_/g;
		##$Article =~s/\+/ /g;
		#$Article=~s/[^A-Za-z0-9~!@#^%&*()-_+={}:<>\[\]\;\"',.\/]/ /g;
		#$Article =~s/[ ]+/ /g;
		
		open output,">".$output;
		print output $input."|fulltext|".$Article."\n\n";
		close output;
	}
	elsif($type eq "PMIDlist")
	{
		
		open output,">".$output;
		close output;
		#open xmls,">".$output.".xmls";
		#close xmls;
		my $count=0;
		open input,"<".$input;
		while(<input>)
		{
			my $pmid=$_;
			$pmid=~s/[\n\r]//g;
			my $efetch="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=".$pmid."&retmode=xml";
			$Article = get($efetch);
			$Article =~s/[\n\r\t]//g;
			decode_entities($Article);
			#$Article =~s/→/ > /g;
			#$Article =~s/↑/\(\+\)/g;
			#$Article =~s/↓/\(\-\)/g;
			#$Article =~s/←/ < /g;
			#$Article =~s/Δ/\/\\/g;
			#$Article =~s/α/a/g;
			#$Article =~s/β/�/g;
			#$Article =~s/×/�/g;
			#$Article =~s/κ/k/g;
			#$Article =~s/γ/y/g;
			#$Article =~s/�/b/g;
			#$Article =~s/�/\-/g;
			#$Article =~s/≥/>= /g;
			#$Article =~s/��/<=/g;
			#$Article =~s/ //g;
			#$Article =~s/�/ /g;
			#$Article =~s/ /   /g;
			#$Article =~s/��/  /g;
			#$Article =~s/�/ /g;
			#$Article =~s/�/ /g;
			#$Article =~s/μ/�/g;
			#$Article =~s/\&lt\;/</g;
			#$Article =~s/\&gt\;/>/g;	
			#$Article =~s/\&quot\;/"/g;
			#$Article =~s/\&[a-z]{1,4}\;/\&/g;
			#$Article =~s/\&[a-z]{1,4}\;/\&/g;
			#$Article =~s/\&/\_/g;
			##$Article =~s/\+/ /g;
			#$Article=~s/[^A-Za-z0-9~!@#^%&*()-_+={}:<>\[\]\;\"',.\/]/ /g;
			#$Article =~s/[ ]+/ /g;
			
			#open xmls,">>".$output.".xmls";
			#print xmls $pmid."\|xml\|".$Article."\n\n";
			#close xmls;
			
			OutputFormat($Article,$output);
			$count++;
			if($count%100==0){print $count."\n";}
			
		}
		close input;
	}
	elsif($type eq "PMCIDlist")
	{
		open output,">".$output;
		close output;
		
		my $count=0;
		open input,"<".$input;
		while(<input>)
		{
			my $pmcid=$_;
			$pmcid=~s/[\n\r]//g;
			my $efetch="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=".$pmcid."";
			$Article = get($efetch);
			$Article =~s/[\n\r\t]//g;
			decode_entities($Article);
			#$Article =~s/→/ > /g;
			#$Article =~s/↑/\(\+\)/g;
			#$Article =~s/↓/\(\-\)/g;
			#$Article =~s/←/ < /g;
			#$Article =~s/Δ/\/\\/g;
			#$Article =~s/α/a/g;
			#$Article =~s/β/�/g;
			#$Article =~s/×/�/g;
			#$Article =~s/κ/k/g;
			#$Article =~s/γ/y/g;
			#$Article =~s/�/b/g;
			#$Article =~s/�/\-/g;
			#$Article =~s/≥/>= /g;
			#$Article =~s/��/<=/g;
			#$Article =~s/ //g;
			#$Article =~s/�/ /g;
			#$Article =~s/ /   /g;
			#$Article =~s/��/  /g;
			#$Article =~s/�/ /g;
			#$Article =~s/�/ /g;
			#$Article =~s/μ/�/g;
			#$Article =~s/\&lt\;/</g;
			#$Article =~s/\&gt\;/>/g;	
			#$Article =~s/\&quot\;/"/g;
			#$Article =~s/\&[a-z]{1,4}\;/\&/g;
			#$Article =~s/\&[a-z]{1,4}\;/\&/g;
			#$Article =~s/\&/\_/g;
			##$Article =~s/\+/ /g;
			#$Article=~s/[^A-Za-z0-9~!@#^%&*()-_+={}:<>\[\]\;\"',.\/]/ /g;
			#$Article =~s/[ ]+/ /g;
			open output,">>".$output;
			print output $pmcid."|fulltext|".$Article."\n\n";
			close output;
			
			$count++;
			if($count%100==0){print $count."\n";}
		}
		close input;
	}
}
sub main
{
	my $type;
	my $input;
	my $output;
	for(my $i=0;$i<@ARGV;$i++)
	{
		if($ARGV[$i] eq "-t")
		{
			$i++;
			$type=$ARGV[$i];
		}
		elsif($ARGV[$i] eq "-i")
		{
			$i++;
			$input=$ARGV[$i];
		}
		elsif($ARGV[$i] eq "-o")
		{
			$i++;
			$output=$ARGV[$i];
		}
		elsif($ARGV[$i]=~/^-s(.+)$/)
		{
			$setup=$1;
		}
		elsif($ARGV[$i]=~/^-i(.+)$/)
		{
			$input=$1;
		}
		elsif($ARGV[$i]=~/^-o(.+)$/)
		{
			$output=$1;
		}
	}
	my %setup_hash=();
	
	if($input eq "")
	{
		print "Input is empty. Please reinsert again.\n";
		print "Instruction Format:\n\n\tperl PreProcessing.pl -t [type:PMID|PMCID|PMIDlist|PMCIDlist] -i [input] -o [output]\n";
		print "\te.g. perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt\n";
		print "\te.g. perl PreProcessing.pl -t PMCID -i 3190010 -o 3190010.txt\n";
		print "PS: XML format inlcudes PubMed abstract or PubMed Centrol full text.\n";
	}
	elsif($output eq "")
	{
		print "Output is empty. Please reinsert again.\n";
		print "Instruction Format:\n\n\tperl PreProcessing.pl -t [type:PMID|PMCID|PMIDlist|PMCIDlist] -i [input] -o [output]\n";
		print "\te.g. perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt\n";
		print "\te.g. perl PreProcessing.pl -t PMCID -i 3190010 -o 3190010.txt\n";
		print "PS: XML format inlcudes PubMed abstract or PubMed Centrol full text.\n";
	}
	elsif($type!~/^(PMID|PMCID|PMIDlist|PMCIDlist)$/)
	{
		print "Type is wrong. Please reinsert again.\n";
		print "Instruction Format:\n\n\tperl PreProcessing.pl -t [type:PMID|PMCID|PMIDlist|PMCIDlist] -i [input] -o [output]\n";
		print "\te.g. perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt\n";
		print "\te.g. perl PreProcessing.pl -t PMCID -i 3190010 -o 3190010.txt\n";
		print "PS: XML format inlcudes PubMed abstract or PubMed Centrol full text.\n";
	}
	elsif($input eq $output)
	{
		print "Input file name cannot the same to Output file name. Please reinsert again.\n";
		print "Instruction Format:\n\n\tperl PreProcessing.pl -t [type:PMID|PMCID|PMIDlist|PMCIDlist] -i [input] -o [output]\n";
		print "\te.g. perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt\n";
		print "\te.g. perl PreProcessing.pl -t PMCID -i 3190010 -o 3190010.txt\n";
		print "PS: XML format inlcudes PubMed abstract or PubMed Centrol full text.\n";
	}
	else
	{
		&GetResultViaEutilities($type,$input,$output);
	}
}

main();
