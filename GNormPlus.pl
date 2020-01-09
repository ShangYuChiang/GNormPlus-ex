#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

sub GNormPlus
{
	my $setup_hash=@_[0];
	my $input=@_[1];
	my $output=@_[2];
	my $filename=@_[3];
	my $dictionary=@_[4];
	my $format=@_[5];

	#### Sentence Segmentation ####
	if($$setup_hash{"Sentence_Segmentation"} eq "True")
	{
		print "Sentence_Segmentation\n";
		require './Library/Sentence_Segmentation.pm';
		if($format eq "BioC")
		{
			require './Library/Sentence_Segmentation_BioC.pm';
			GNormPlus::Sentence_Segmentation_BioC($input,$filename,$dictionary);
		}
		elsif($format eq "PubTator")
		{
			GNormPlus::Sentence_Segmentation_PubTator($input,$filename,$dictionary);
		}
		else #XML/FreeText
		{
			GNormPlus::Sentence_Segmentation($input,$filename,$dictionary);
		}
		GNormPlus::Ab3PExt($input,$filename,$dictionary); #Abbreviation resolution
	}

	#### Species Name Recognition: SR4GN  ####
	if($$setup_hash{"Species_Name_Recognition"} eq "True")
	{
		print "Species_Name_Recognition\n";
		require './Library/Species_Name_Recognition.pm';
		GNormPlus::Species_Name_Recognition($input,$filename,$dictionary);
		GNormPlus::CellLine_Recognition($input,$filename,$dictionary);
		GNormPlus::Disambiguation_occurrence($input,$filename,$dictionary);
		GNormPlus::Filtering($input,$filename,$dictionary);
	}

	#### Gene Name Recognition ####
	#1.Extracting gene mentions by CRF++ model
	#2.Extracting possible gene identifiers. Using gene/protein identifiers(Such as SwisSprot ID) to identify potential mentions. 
	if($$setup_hash{"Gene_Name_Recognition"} eq "True")
	{	
		print "Gene_Name_Recognition\n";
		require './Library/Gene_Name_Recognition.pm';
		GNormPlus::Gene_Name_Recognition($input,$filename,$dictionary);
		GNormPlus::ID_Extraction($input,$filename,$dictionary);		
	}

	#### Species Assignment: SR4GN ####
	#Assigning Taxonomy ID to gene mentions.(excluding Candidate_ID)
	if($$setup_hash{"Species_Assignment"} eq "True")
	{
		print "Species_Assignment\n";
		require './Library/Species_Assignment.pm';
		GNormPlus::Species_Assignment($input,$filename,$dictionary,$$setup_hash{"Focus_Species"});
	}

	####  Simplifying Composite Genes Mentions: SimConcept ####
	if($$setup_hash{"SimConcept"} eq "True")
	{
		print "SimConcept\n";
		require './Library/SimConcept.pm';
		GNormPlus::SimConcept($input,$filename,$dictionary);
	}

	#### Gene Normalization ####
	if($$setup_hash{"Gene_Normalization"} eq "True")
	{
		print "Gene_Normalization\n";
		require './Library/Gene_Normalization.pm';
		GNormPlus::Gene_Normalization($input,$filename,$dictionary);
		GNormPlus::Gene_Mapping($output,$filename,$dictionary);
	}

	#### Output Result ####
	if($format eq "PubTator") #PubTator
	{
		print "Output\n";
		require './Library/Output_Result.pm';
		GNormPlus::Output_Result_PubTator($input,$output,$filename,$dictionary);
	}
	elsif($format eq "BioC") #BioC
	{
		print "Output\n";
		require './Library/Output_Result_BioC.pm';
		GNormPlus::Output_Result_BioC($input,$output,$filename,$dictionary);
	}
	else{}#XML/FreeText

	open clear,"|rm tmp/".$filename.".*";	close clear;
	#open clear,"|rm input/".$filename;	close clear;
}

sub main
{
	my $setup;
	my $input;
	my $output;
	my $dictionary;
	my $FocusSpecies="";
	
	for(my $i=0;$i<@ARGV;$i++)
	{
		if($ARGV[$i] eq "-s")
		{
			$i++;
			$setup=$ARGV[$i];
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
		elsif($ARGV[$i] eq "-d")
		{
			$i++;
			$dictionary=$ARGV[$i];
		}
		elsif($ARGV[$i] eq "-sp")
		{
			$i++;
			$FocusSpecies=$ARGV[$i];
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
		elsif($ARGV[$i]=~/^-d(.+)$/)
		{
			$dictionary=$1;
		}
		elsif($ARGV[$i]=~/^-sp(.+)$/)
		{
			$FocusSpecies=$1;
		}
	}
	my %setup_hash=();
	
	if($input eq "")
	{
		print "Instruction Format:\n\n\tperl GNormPlus.pl -i [input dir] -o [output dir]\n";
		print "\te.g. perl GNormPlus.pl -i input -o output\n";
	}
	elsif($output eq "")
	{
		print "Instruction Format:\n\n\tperl GNormPlus.pl -i [input dir] -o [output dir]\n";
		print "\te.g. perl GNormPlus.pl -i input -o output\n";
	}
	else
	{
		if($dictionary eq ""){$dictionary="dictionary";}
		if($setup eq "") # default: Using all modules
		{
			$setup_hash{"Sentence_Segmentation"}="True";
			$setup_hash{"Species_Name_Recognition"}="True";
			$setup_hash{"Gene_Name_Recognition"}="True";
			$setup_hash{"Species_Assignment"}="True";
			$setup_hash{"SimConcept"}="True";
			$setup_hash{"Gene_Normalization"}="True";
		}
		else
		{
			open setup,"<$setup";
			while(<setup>)
			{
				my $setup=$_;
				if ($setup=~/(\w+)[\t ]+=[\t ]+(\w+)/)
				{
					$setup_hash{$1}=$2;
				}
			}
			close setup;
		}
		
		if($setup_hash{"Focus_Species"}!~/^[0-9]+$/){$setup_hash{"Focus_Species"}="All";}
		if($FocusSpecies ne "" && $FocusSpecies=~/^[0-9]+$/){$setup_hash{"Focus_Species"}=$FocusSpecies;}
		
		#### [Running module List] ####
		print "[Running module List]:\n\n";
		if($setup_hash{"Sentence_Segmentation"} eq "True") {print "[Module]: Sentence segmentation\n";}
		if($setup_hash{"Species_Name_Recognition"} eq "True") {print "[Module]: Species name recognition (SR) module\n";}
		if($setup_hash{"Gene_Name_Recognition"} eq "True") {print "[Module]: Gene name recognition (GNR) module\n";}
		if($setup_hash{"Species_Assignment"} eq "True") {print "[Module]: Species assignment (SA) module\n";}
		if($setup_hash{"SimConcept"} eq "True") {print "[Module]: SimConcept module\n";}
		if($setup_hash{"Gene_Normalization"} eq "True") {print "[Module]: Gene normalization(GN) module\n";}
		print "Focus Species taxonomy ID: ".$setup_hash{"Focus_Species"}."\n";
		print "\nStart....\n";
		
		$counter=1;
		opendir(DIR, $input);
		@class = grep(/[a-z0-9]/,readdir(DIR));
		closedir(DIR);
		foreach my $filename(@class)
		{
			my ($sec1,$min1,$hour1,$day1,$mon,$year)=localtime(time);
			
			#detect format
			my $format="XML_FreeText";
			open input,"<".$input."/".$filename;
			while(<input>)
			{
				my $tmp=$_;
				if($tmp=~/^.+\|t\|/)
				{
					$format="PubTator";
					last;
				}
				elsif($tmp=~/<collection>/)
				{
					$format="BioC";
					last;
				}
			}
			close input;
			
			my $format_failed="N";
			if($format eq "BioC")
			{
				my $STR="";
				open input,"<".$input."/".$filename;
				while(<input>)
				{
					my $tmp=$_;
					$STR=$STR.$tmp." ";
				}
				close input;
				if($STR =~/<text><\/text>/)
				{
					print "(".$counter.") File:".$filename."\tFormat is failed. At least one passage text is empty. This file wouldn't be processed.\n";
					$format_failed="Y";
				}
			}
			
			if($format_failed eq "N")
			{
				GNormPlus(\%setup_hash,$input,$output,$filename,$dictionary,$format);
				my ($sec2,$min2,$hour2,$day2,$mon,$year)=localtime(time);
				my $hour=($day2-$day1)*24+($hour2-$hour1);
				my $min=$hour*60+($min2-$min1);
				my $timecost=$min*60+($sec2-$sec1);
				print "(".$counter.") File:".$filename."\tFinished in $timecost sec.\n";
			}
			$counter++;
		}
		close pmidlist;
	}
}

main();