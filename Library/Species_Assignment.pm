#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
# Description: Assigns species identifier(Taxonomy ID) to each gene mention
#===================================================

package GNormPlus;

require './Library/Text/English.pm';
my %species_prefix_hash=();
	
sub Species_Assignment
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	my ($Focus_Species)=@_[3];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $gmention_extraction="tmp/".$filename.".gmention.xml";
	my $tax_extraction="tmp/".$filename.".tax.xml";
	my $sa_extraction="tmp/".$filename.".sa.xml";
	
	open output,">$sa_extraction";
	close output;
	
	my %pmidhash=();
	my $pmidcount=0;
	open input,"<$sentence_extraction";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/<TEXT pmid='([^\']+)'/)
		{
			$pmidhash{$1}=$pmidcount++;
		}
	}
	close input;
	my @pmidArr = sort {$pmidhash{$a} <=> $pmidhash{$b}} keys %pmidhash;
	
	open coefficient,"<$dictionary/coefficient.txt";
	while(<coefficient>)
	{
		$coefficient=$_;
		if( $coefficient=~/^(.*)	(.*)	(.*)	(.*)	(.*)$/ )
		{
			$coefficient_hash{$1}=$2."|".$3."|".$4."|".$5;
		}
	}
	
	#Read Ab3P
	my %Abb_hash=();
	my $pmid="";
	open Abb,"<tmp/".$filename.".Ab3P";
	while(<Abb>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^  (.+)\|(.+)\|/)
		{
			my $Abb=lc($1);
			$Abb=~s/[\W\-\_]//g;
			$Abb_hash{$Abb}=1;
		}
	}
	close Abb;
	
	my %wholefreq_hash=();
	open wholefreq,"<$dictionary/taxonomy_freq.txt";
	while(<wholefreq>)
	{
		$taxonomy_freq=$_;
		if( $taxonomy_freq=~/^(.*)	(.*)$/ )
		{
			$wholefreq_hash{$1}=$2/17461751;
		}
	}
	
	foreach my $pmid (@pmidArr)
	{
		%species_prefix_hash=();
		my %voting_tax_hash=();
		my %coefficient_hash=();
		my %tax_location_hash=();
		my %species_location_hash=();
		
		open tax,"<$tax_extraction";
		while(<tax>)
		{
			$tax=$_;
			if( $tax=~/<Tax pmid='$pmid' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.+)'>(.+)<\/Tax>/)
			{
				my $sid=$1;
				my $start=$2;
				my $last=$3;
				my $tax_id=$4;
				my $species_mention=$5;
				my $species_mention_tmp=lc($5);
				$species_mention_tmp=~s/[\W\-\_]//g;
				$species_location_hash{$sid."|".$start."|".$last}=1;
				if(not exists $Abb_hash{$species_mention_tmp}) #Not use abbreviation on voting
				{
					if($sid=~/(\_1$|title)/i){$voting_tax_hash{$tax_id}=$voting_tax_hash{$tax_id}+2;} #title has double weight
					else{$voting_tax_hash{$tax_id}=$voting_tax_hash{$tax_id}+1;}
				}
				if($species_mention=~/^(\S+) (\S+)$/)
				{
					my $pre1=substr($1,0,1);
					if($pre1=~/\W\-\_/)
					{
						$pre1=substr($1,1,1);
					}
					my $pre2=substr($2,0,1);
					$species_prefix_hash{$pre1.$pre2}=$tax_id;
					$species_prefix_hash{uc($pre1).lc($pre2)}=$tax_id;
				}
				if(not exists $Abb_hash{$species_mention_tmp}) #Not use abbreviation on occurring
				{
					if (exists $tax_location_hash{$sid})
					{
						$tax_location_hash{$sid}=$tax_location_hash{$sid}."@".$start."|".$tax_id;
					}
					else
					{
						$tax_location_hash{$sid}=$start."|".$tax_id;
					}
				}
			}
		}
		foreach $tax(keys %voting_tax_hash)
		{
			$voting_tax_hash{$tax}=$voting_tax_hash{$tax}+$wholefreq_hash{$tax};
		}
		
		my @rank = reverse sort {$voting_tax_hash{$a} <=> $voting_tax_hash{$b}} keys %voting_tax_hash;

		$majority_voting_species = $rank[0];
		
		#SR4GN - SRI (Species Representation Indicator)
=pod	
		if (!($majority_voting_species =~/^[0-9]+$/))
		{
			my $abstract="";
			open sentence,"<$sentence_extraction";
			while(<sentence>)
			{
				my $tmp=$_;
				if( $tmp=~/^<TEXT pmid='.+' sid='.+'>(.+)<\/TEXT>$/)
				{
					my $sentence=$1;
					$sentence=lc($sentence);
					$sentence=~s/[\W\(\)\[\]\_\-\+]/ /g;
					$sentence=~s/[ ]+/ /g;
					$abstract=$abstract.$sentence;
				}
			}
			
			my %word_hash=();
			my @split_tokenization_entity = split(" ",$abstract);
			my @tmp_split;
			foreach my $tmp (@split_tokenization_entity)
			{
				if (not exists $stopword_hash{$tmp})
				{
					push @tmp_split,$tmp;
				}
			}
			@split_tokenization_entity=@tmp_split;
			@split_tokenization_entity=Text::English::stem(@split_tokenization_entity);
			foreach $tmp (@split_tokenization_entity)
			{
				$word_hash{$tmp}=$word_hash{$tmp}+1;
			}
			
			my $sum_9606=$sum_7227=$sum_4932=$sum_10090=0;
			foreach $word(keys %word_hash)
			{
				if($coefficient_hash{$word} =~/^(.*)\|(.*)\|(.*)\|(.*)$/)
				{
					$z_9606=$1;
					$z_7227=$2;
					$z_4932=$3;
					$z_10090=$4;
					$sum_9606=$sum_9606+$z_9606*$word_hash{$word};
					$sum_7227=$sum_7227+$z_7227*$word_hash{$word};
					$sum_4932=$sum_4932+$z_4932*$word_hash{$word};
					$sum_10090=$sum_10090+$z_10090*$word_hash{$word};
				}
			}
			if($sum_9606>0){$majority_voting_species="SRI9606";}
			elsif($sum_7227>0 && $sum_7227>$sum_9606 && $sum_7227>$sum_4932 && $sum_7227>$sum_10090){$majority_voting_species="SRI7227";}
			elsif($sum_4932>0 && $sum_4932>$sum_9606 && $sum_4932>$sum_7227 && $sum_4932>$sum_10090){$majority_voting_species="SRI4932";}
			elsif($sum_10090>0 && $sum_10090>$sum_7227 && $sum_10090>$sum_4932 && $sum_10090>$sum_9606){$majority_voting_species="SRI10090";}
			close sentence;
		}
=cut

		open gmention,"<$gmention_extraction";
		while(<gmention>)
		{
			$tmp=$_;
			if($tmp=~/<(.+) pmid='$pmid' sid='(.+)' start='(.+)' end='(.+)' score='(.+)'>(.+)<\/.+>/)
			{
				my $gmention_type = $1;
				my $sid = $2;
				my $start = $3;
				my $last = $4;
				my $score = $5;
				my $entity = $6; 
				if( not exists $species_location_hash{$sid."|".$start."|".$last} )
				{
					#Simple rule-based corrdination 
=pod
					if($entity =~ /^(.+)(\-*[0-9]+)[ ]*,[ ]*(\-*[0-9]+)[ ]*(,|,[ ]*and|and|,[ ]*or|or)[ ]*(\-*[0-9]+)$/) #MMP-2, 7, 9
					{
						$tmp1=$1;
						$num1=$2;
						$num2=$3;
						$num3=$5;
						$tmp_entity=$tmp1."".$num1;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
						$tmp_entity=$tmp1."".$num2;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
						$tmp_entity=$tmp1."".$num3;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
					}
					elsif($entity =~ /^(.+?[^0-9])(\-*[0-9]+)[ ]*([,\/]|and|or)[ ]*(\-*[0-9]+)$/ || $entity =~ /^(.+?)(\-*[A-Z])[ ]*([,\/]|and|or)[ ]*(\-*[A-Z])$/) #XRCC2/3 || CDKN2A/B
					{
						$tmp1=$1;
						$num1=$2;
						$num2=$4;
						$tmp_entity=$tmp1."".$num1;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
						$tmp_entity=$tmp1."".$num2;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
					}
					elsif($entity =~ /^(.+[^0-9])(\-*[A-Z][0-9])[ ]*([,\/]|and|or)[ ]*(\-*[A-Z][0-9])$/) #XRCC2/3 || CDKN2A/B
					{
						$tmp1=$1;
						$num1=$2;
						$num2=$4;
						$tmp_entity=$tmp1."".$num1;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
						$tmp_entity=$tmp1."".$num2;
						insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
					}
					elsif($entity =~ /^(.+[^0-9])([0-9]+)\-([0-9]+)$/)
					{
						$tmp1=$1;
						$num1=$2;
						$num2=$3;
						if(($num2-$num1)<=10)
						{
							for($i=$num1;$i<=$num2;$i++)
							{
								$tmp_entity=$tmp1."".$i;
								insertion($pmid,$gmention_type,$sid,$tmp_entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,$entity,$Focus_Species);
							}
						}
					}
=cut		
					insertion($pmid,$gmention_type,$sid,$entity,$start,$last,$majority_voting_species,$tax_extraction,$tax_location_hash{$sid},$sa_extraction,"",$Focus_Species);
				}
			}
		}	
	}
	return 1;
}
sub insertion
{
	my ($pmid) = @_[0];
	my ($gmention_type) = @_[1];
	my ($sid) = @_[2];
	my ($org_entity) = @_[3];
	my ($start_site) = @_[4];
	my ($end_site) = @_[5];
	my ($majority_voting_species) = @_[6];
	my ($tax_extraction) = @_[7];
	my ($tax_location) = @_[8];
	my ($sa_extraction) = @_[9];
	my ($Org_numEnt) = @_[10];
	my ($Focus_Species) = @_[11];
	
	my $ORG=$org_entity;
	if($Org_numEnt ne ""){$ORG=$Org_numEnt};
	
	$entity = $org_entity;
	$entity =~ s/\-/ /g;
	$entity =~ s/[\[\]\'!@\#\$\&\%:;.,()\\=*><\"]/ /g;
	$entity =~ s/([A-Za-z]+)([0-9]+)/$1 $2/g;
	$entity =~ s/([0-9]+)([A-Za-z]+)/$1 $2/g;
	$entity =~ s/[ ]+/ /g;
	my %species_hash=();
	
	open output,">>$sa_extraction";
	
	open tax,"<$tax_extraction";
	while(<tax>)
	{
		my $tmp=$_;
		if( $tmp=~/<Tax pmid='$pmid' .* tax_id='([0-9]+)'>(.+)<\/Tax>/)
		{
			my $tax_id=$1;
			my $species_mention=$2;
			$species_mention =~ s/\-/ /g;
			$species_mention =~ s/[\[\]\'!@\#\$\&\%:;.,()\\=*><\"]/ /g;
			$species_mention =~ s/([A-Za-z]+)([0-9]+)/$1 $2/g;
			$species_mention =~ s/([0-9]+)([A-Za-z]+)/$1 $2/g;
			$species_mention =~ s/[ ]+/ /g;
			$species_hash{$species_mention} = $tax_id;
		}
	}
	close tax;

	my $mode=0;

	#prefix
	if($org_entity =~/^[h]([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='9606' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}	
	elsif($org_entity =~/^[r]([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='10116' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	elsif($org_entity =~/^[m]([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='10090' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	elsif($org_entity =~/^[y]([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='4932' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	elsif($org_entity =~/^[d]([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='7227' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	elsif($org_entity =~/^(z|[Zz]f)([A-Z].*)$/)
	{
		$clear_entity=$2;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='7955' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	elsif($org_entity =~/^[Aa]t([A-Z].*)$/)
	{
		$clear_entity=$1;
		if($Focus_Species=~/^[0-9]+$/)
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$Focus_Species' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		else
		{
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='3702' org='$ORG'>$clear_entity</$gmention_type>\n";
		}
		$mode=1;
	}
	foreach my $prefix (keys %species_prefix_hash)
	{
		if($prefix=~/^[A-Za-z]+$/ && $org_entity =~/^$prefix([A-Z].*)$/)
		{
			$clear_entity=$1;
			my $target="";
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			else
			{
				$target=$species_prefix_hash{$prefix};
			}
			print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$clear_entity</$gmention_type>\n";
			$mode=1;
		}
	}

	#previous 
	if($mode == 0)
	{
		foreach $each_species (keys %species_hash)
		{
			if($org_entity =~/^$each_species (.*)$/i || $org_entity =~/[\W\_\-]$each_species (.*)$/i)
			{
				$clear_entity=$1;
				if($clear_entity !~/^[0-9]+$/ && $clear_entity !~/^[IV]+$/ )
				{
					my $target="";
					if($Focus_Species=~/^[0-9]+$/)
					{
						$target=$Focus_Species;
					}
					else
					{
						$target=$species_hash{$each_species};
					}
					print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$clear_entity</$gmention_type>\n";
					$mode=1;
				}
			}
		}
	}
	
	#co_occurring 
	if($mode == 0)
	{
		my $target="";
		my $target_species_start=0;
		my @array_speices_loca=split("@",$tax_location);
		foreach my $split_species_loca(@array_speices_loca)
		{
			if($split_species_loca=~/(.+)\|(.+)/)
			{
				my $species_start=$1;
				my $tax_id=$2;
				if($species_start<$start_site && $species_start>=$target_species_start)
				{
					$target=$tax_id;
					$target_species_start=$species_start;
				}
			}
		}
		if(!($target eq ""))
		{
			$mode=1;
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			if($Org_numEnt ne "")
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$org_entity</$gmention_type>\n";
			}
			else
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target'>$org_entity</$gmention_type>\n";
			}
		}
	}
	if($mode == 0)
	{
		my $target="";
		my $target_species_start=1000;
		my @array_speices_loca=split("@",$tax_location);
		foreach my $split_species_loca(@array_speices_loca)
		{
			if($split_species_loca=~/(.+)\|(.+)/)
			{
				my $species_start=$1;
				my $tax_id=$2;
				if($species_start>$start_site && $species_start<$target_species_start)
				{
					$target=$tax_id;
					$target_species_start=$species_start;
				}
			}
		}
		if(!($target eq ""))
		{
			$mode=1;
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			if($Org_numEnt ne "")
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$org_entity</$gmention_type>\n";
			}
			else
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target'>$org_entity</$gmention_type>\n";
			}
		}
	}
	
	#focus_species 
	if($mode == 0)
	{
		if ($majority_voting_species =~/^SRI([0-9]+)$/)
		{
			my $target="";
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			else
			{
				$target="9606";
			}
			if($Org_numEnt ne "")
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$org_entity</$gmention_type>\n";
			}
			else
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target'>$org_entity</$gmention_type>\n";
			}
		}
		elsif ($majority_voting_species =~/^[0-9]+$/)
		{
			my $target="";
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			else
			{
				$target=$majority_voting_species;
			}
			if($Org_numEnt ne "")
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$org_entity</$gmention_type>\n";
			}
			else
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target'>$org_entity</$gmention_type>\n";
			}
		}
		else
		{			
			my $target="";
			if($Focus_Species=~/^[0-9]+$/)
			{
				$target=$Focus_Species;
			}
			else
			{
				$target="9606"; #default setting.
			}
			if($Org_numEnt ne "")
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target' org='$ORG'>$org_entity</$gmention_type>\n";
			}
			else
			{
				print output "<$gmention_type pmid='$pmid' sid='$sid' start='$start_site' end='$end_site' tax_id='$target'>$org_entity</$gmention_type>\n";
			}
		}
	}
	
	close output;
	return 1;
}
return 1;