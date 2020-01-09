#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

require './Library/Text/English.pm';
require './Library/Regex/PreSuf.pm';
require './Library/Lingua/EN/Tagger.pm';

my $p = new Lingua::EN::Tagger;

sub round
{
  my($value,$rank) = @_;
  if($value>0){
    return int($value * 10**$rank + 0.5) / 10**$rank;
  }else{
    return int($value * 10**$rank - 0.4 )/ 10**$rank;
  }
}

sub FeatureExtraction
{
	my $filename = $_[0];
	my $inputfile = $_[1];
	my $TrainorTest = $_[2];
	
	my @stemmed_tokens;
	my @POS_tokens;
	
	my %Gene_CTD_locationmap_hash=();
	my %Gene_CTD_lastmap_hash=();
	my %CTD_result_hash=();
	
	my %abb_locationmap_hash=();
	my %abb_lastmap_hash=();
	my %abb_Mention2Type_hash=();
	my %abb_result_hash=();
	my %full_locationmap_hash=();
	my %full_lastmap_hash=();
	my %full_result_hash=();
	my %Result_hash=();
	my %locationmap_hash=();
	my %lastmap_hash=();
	my %Mention2Type_hash=();
	my %type_hash=();
	my $mention_Regex_region_CTD=0;
	my $mention_Regex_region_abb=0;
	my $mention_Regex_region_full=0;
	my $sid="";
	
	my %Begin_hash=();	
	my %Inside_hash=();
	my %End_hash=();
	$Begin_hash{"Gene"}="G";	$Inside_hash{"Gene"}="E";	$End_hash{"Gene"}="N";
	$Begin_hash{"FamilyName"}="F";	$Inside_hash{"FamilyName"}="A";	$End_hash{"FamilyName"}="M";
	$Begin_hash{"DomainMotif"}="D";	$Inside_hash{"DomainMotif"}="I";	$End_hash{"DomainMotif"}="T";
	$Begin_hash{"Cell"}="C";	$Inside_hash{"Cell"}="H";	$End_hash{"Cell"}="L";
	$Begin_hash{"ChromosomeLocation"}="X";	$Inside_hash{"ChromosomeLocation"}="Y";	$End_hash{"ChromosomeLocation"}="Z";

	my %PostProessStopWord_hash=();
	open stopword,"<./Library/Regex/CTDStopWord.txt";	
	while(<stopword>)
	{
		my $tmp=lc($_);
		$tmp=~s/[\n\r\W\-\_]//g;
		$PostProessStopWord_hash{$tmp}=1;
	}
	close stopword;
	
	my %abb_hash=();
	my %full_hash=();
	my %pmid_abb_hash=();
	my %pmid_full_hash=();
	my $pmid="";
	open Abb,"<tmp/".$filename.".Ab3P";
	while(<Abb>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\|]+)$/)
		{
			$pmid=$1;
			$tmp=<Abb>; #sentences
			%abb_hash=();
			%full_hash=();			
		}
		elsif($tmp=~/^  (.+)\|(.+)\|/)
		{
			my $Abb=$1;
			my $FN=$2;
			$abb_hash{$Abb}=length($Abb);
			$full_hash{$FN}=length($FN);
		}
		elsif(length($tmp)==0)
		{
			my @abb_array = reverse sort {$abb_hash{$a} <=> $abb_hash{$b}} keys %abb_hash;
			my $abb_RegEx = Regex::PreSuf::presuf(@abb_array);
			my @full_array = reverse sort {$full_hash{$a} <=> $full_hash{$b}} keys %full_hash;
			my $full_RegEx = Regex::PreSuf::presuf(@full_array);
			if(@abb_array>0)
			{
				$pmid_abb_hash{$pmid}=$abb_RegEx;
			}
			if(@full_array>0)
			{
				$pmid_full_hash{$pmid}=$full_RegEx;
			}
		}
	}
	close Abb;
	
	open CTD_gene,"<./Library/Regex/CTD_gene.RegEx";
	my $CTD_gene_RegEx=<CTD_gene>;
	close CTD_gene;
	
	my $sentence="";
	my $article="";
	$pmid="";
	my $sid="";
	
	open input,"<".$inputfile;
	open output,">".$inputfile.".data";
	close output;
	open locationA,">".$inputfile.".location";
	close locationA;
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^<TEXT pmid='([^<>]+)' sid='([^<>]+)'>(.+)<\/TEXT>$/)
		{
			my $pmid=$1;
			my $sid=$2;
			my $sentence=$3;
			while($sentence =~/^(.*)(<.+?>)(.*)$/)
			{
				my $pre=$1;
				my $tag=$2;
				my $post=$3;
				$tag=~s/./ /g;
				$sentence=$pre.$tag.$post;
			}
			$sentence=~s/[\n\r]//g;
			
			#CTD_gene_RegEx
			my $text=" ".$sentence." ";
			while (	$text=~/^(.*?W)(($CTD_gene_RegEx)(s|)[\W\-\_]*(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega))(\W.*)$/i ||
					$text=~/^(.*\W)(($CTD_gene_RegEx)(s|)[\W\-\_]*([0-9ABCIV\-\/, ]| and| or|){0,})(\W.*)$/i )
			{
				my $pre=$1;
				$pre=substr($pre,1,);
				my $mention=$2;
				$mention=~s/[\., ]+$//g; 
				$mention=~s/[\W\-\_]+(and|or)$//g;
				my $post=substr($text,length($pre)+length($mention)+1,);
				my $tmp=$mention;
				$tmp=~s/./@/g;
				$text=" ".$pre.$tmp.$post;
				if(lc($mention)!~/^[0-9\W\-\_]+$/i && (not exists $stopword_hash{lc($mention)}))
				{
					$Gene_CTD_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$Gene_CTD_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
				}
			}
			$text=" ".$sentence." ";
			while (	$text=~/^(.*\W(human|rat|rats|mouse|mice|fish|cat|dog|fly|yeast|coli|murine|elegans|rosophila) )([^@ ]+)(( gene| protein|[\,\.])\W.*)$/i ||
					$text=~/^(.*\Wgene[ ]+)([^@ ]+)([ ]*[\,\.]\W.*)$/i
					)
			{
				my $pre=$1;
				$pre=substr($pre,1,);
				my $mention=$2;
				my $tmp=$2;
				my $post=substr($text,length($pre)+length($mention)+1,);
				$tmp=~s/./@/g;
				$text=" ".$pre.$tmp.$post;
				if(lc($mention)!~/^[0-9\W\-\_]+$/i && not exists $PostProessStopWord_hash{lc($mention)})
				{
					$Gene_CTD_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$Gene_CTD_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
				}
			}
			
			#abb_RegEx
			$text=" ".$sentence." ";
			my $abb_RegEx=$pmid_abb_hash{$pmid};
			if($abb_RegEx)
			{
				$abb_RegEx=~s/([\W\-\_])/\\$1/g;
				while ($text=~/^(.*\W)($abb_RegEx)(\W.*)$/i)
				{
					my $pre=$1;
					$pre=substr($pre,1,);
					my $mention=$2;
					my $tmp=$2;
					my $post=substr($text,length($pre)+length($mention)+1,length($text)-length($pre)-length($mention));
					$tmp=~s/./@/g;
					$text=" ".$pre.$tmp.$post;
					$abb_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$abb_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
				}
			}
			
			#full_RegEx
			$text=" ".$sentence." ";
			my $full_RegEx=$pmid_full_hash{$pmid};
			if($full_RegEx)
			{
				$full_RegEx=~s/([\W\-\_])/\\$1/g;
				while ($text=~/^(.*\W)($full_RegEx)(\W.*)$/i)
				{
					my $pre=$1;
					$pre=substr($pre,1,);
					my $mention=$2;
					my $tmp=$2;
					my $post=substr($text,length($pre)+length($mention)+1,length($text)-length($pre)-length($mention));
					$tmp=~s/./@/g;
					$text=" ".$pre.$tmp.$post;
					$full_locationmap_hash{length($pre)+1}=(length($pre)+1)." ".(length($pre)+length($mention));
					$full_lastmap_hash{length($pre)+1}=(length($pre)+length($mention));
				}
			}
			
			#recording token's location
			my %LastTokenPosition_hash=();
			my %CurrentStart_hash=();
			my %CurrentLast_hash=();
			my %NextTokenPosition_hash=();
			open locationA,">>".$inputfile.".location";
			my $tmp = $sentence." ";
			my $count_token=0;
			my $start=0;
			my $mention_region=0;
			while($tmp=~/^([A-Z]+)(.*)$/ || $tmp=~/^([a-z]+)(.*)$/ || $tmp=~/^([0-9]+)(.*)$/ || $tmp=~/^([\W\-\_])(.*)$/ )
			{
				my $pre=$1;
				my $post=$2;
				my $end=length($pre)+$start;
				
				if ($pre ne " ")
				{
					#for location
					if(exists $locationmap_hash{$start+1})
					{
						$mention_region=$start+1;
						my $state="I";
						if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
						$Result_hash{$count_token}=$Begin_hash{$type_hash{$locationmap_hash{$mention_region}}};
						print locationA $pmid."	".$sid."	".$pre."	".($start+1)."	".$end."	".$Result_hash{$count_token}."\n";
						$LastTokenPosition_hash{$count_token+1}=$end;
						$NextTokenPosition_hash{$count_token-1}=$start+1;
						$CurrentStart_hash{$count_token}=$start+1;
						$CurrentLast_hash{$count_token}=$end;
					}
					elsif($mention_region!=0)
					{
						if($lastmap_hash{$mention_region}>$end)
						{
							my $state="I";
							if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
							$Result_hash{$count_token}=$Inside_hash{$type_hash{$locationmap_hash{$mention_region}}};
							print locationA $pmid."	".$sid."	".$pre."	".($start+1)."	".$end."	".$Result_hash{$count_token}."\n";
							$LastTokenPosition_hash{$count_token+1}=$end;
							$NextTokenPosition_hash{$count_token-1}=$start+1;
							$CurrentStart_hash{$count_token}=$start+1;
							$CurrentLast_hash{$count_token}=$end;
						}
						elsif($lastmap_hash{$mention_region}==$end)
						{
							my $state="I";
							if(exists $character_hash{$start+1}){$state=$character_hash{$start+1};}
							$Result_hash{$count_token}=$End_hash{$type_hash{$locationmap_hash{$mention_region}}};
							print locationA $pmid."	".$sid."	".$pre."	".($start+1)."	".$end."	".$Result_hash{$count_token}."\n";
							$LastTokenPosition_hash{$count_token+1}=$end;
							$NextTokenPosition_hash{$count_token-1}=$start+1;
							$CurrentStart_hash{$count_token}=$start+1;
							$CurrentLast_hash{$count_token}=$end;
						}
						else
						{
							if(not exists $Result_hash{$count_token-2})
							{
								$Result_hash{$count_token-1}=$Begin_hash{$type_hash{$locationmap_hash{$mention_region}}};
							}
							else
							{
								$Result_hash{$count_token-1}=$End_hash{$type_hash{$locationmap_hash{$mention_region}}};
							}
							print locationA $pmid."	".$sid."	".$pre."	".($start+1)."	".$end."\n";
							$LastTokenPosition_hash{$count_token+1}=$end;
							$NextTokenPosition_hash{$count_token-1}=$start+1;
							$CurrentStart_hash{$count_token}=$start+1;
							$CurrentLast_hash{$count_token}=$end;
							$mention_region=0;
						}
					}
					else
					{
						print locationA $pmid."	".$sid."	".$pre."	".($start+1)."	".$end."\n";
						$LastTokenPosition_hash{$count_token+1}=$end;
						$NextTokenPosition_hash{$count_token-1}=$start+1;
						$CurrentStart_hash{$count_token}=$start+1;
						$CurrentLast_hash{$count_token}=$end;
					}
					
					#Gene CTD RegEX
					if(exists $Gene_CTD_locationmap_hash{$start+1})
					{
						$mention_Regex_region_CTD=$start+1;
						$CTD_result_hash{$count_token}="CTD_gene";#Begin
					}
					elsif($mention_Regex_region_CTD!=0)
					{
						if($Gene_CTD_lastmap_hash{$mention_Regex_region_CTD}>=$end)
						{
							$CTD_result_hash{$count_token}="CTD_gene";#Inside
						}
						else
						{
							$CTD_result_hash{$count_token-1}="CTD_gene";#End
							$mention_Regex_region_CTD=0;
						}
					}
					
					#abb RegEX
					if(exists $abb_locationmap_hash{$start+1})
					{
						$mention_Regex_region_abb=$start+1;
						$abb_result_hash{$count_token}="Abbreviation";#Begin
					}
					elsif($mention_Regex_region_abb!=0)
					{
						if($abb_lastmap_hash{$mention_Regex_region_abb}>=$end)
						{
							$abb_result_hash{$count_token}="Abbreviation";#Inside
						}
						else
						{
							$abb_result_hash{$count_token-1}="Abbreviation";#End
							$mention_Regex_region_abb=0;
						}
					}
					
					#full RegEX
					if(exists $full_locationmap_hash{$start+1})
					{
						$mention_Regex_region_full=$start+1;
						$full_result_hash{$count_token}="Fullname";#Begin
					}
					elsif($mention_Regex_region_full!=0)
					{
						if($full_lastmap_hash{$mention_Regex_region_full}>=$end)
						{
							$full_result_hash{$count_token}="Fullname";#Inside
						}
						else
						{
							$full_result_hash{$count_token-1}="Fullname";#End
							$mention_Regex_region_full=0;
						}
					}
					
					$count_token++;
				}
				
				$tmp=$post;
				$start=$end;
			}
			print locationA "\n";
			close locationA;

			#Featuring
			$text=" ".$sentence." ";
			$text=~s/([0-9])([A-Za-z])/$1 $2/g;
			$text=~s/([A-Z])([a-z])/$1 $2/g;
			$text=~s/([a-z])([A-Z])/$1 $2/g;
			$text=~s/([A-Za-z])([0-9])/$1 $2/g;
			$text=~s/([\W\-\_])/ $1 /g;
			$text =~ s/[ ]+/ /g;
			
			#tokens
			my @tokens = split(" ",$text);
			
			#Stemming
			my @stemmed_tokens = Text::English::stem( @tokens ); 
			
			#POS	
			my @POS_tokens=();
			my $tagged_text = $p->add_tags( $text ); 
			my $count_pos=0;
			while($tagged_text =~/^<([a-z]+?)>(.+?)<\/(.+?)>(.+)$/)
			{
				$postag=$1;
				$POS_tokens[$count_pos]=$postag;
				$count_pos++;
				$tagged_text=$4;
				$tagged_text =~ s/^ //g;
			}
			
			open outputA,">>".$inputfile.".data";
			$count_token=0;
			my $last_token="";
			my $next_token="";
			foreach $token(@tokens)
			{
				if($count_token<@tokens)
				{
					$next_token=$tokens[$count_token+1];
				}
				$insert_query="";
				if($stemmed_tokens[$count_token] eq "")
				{
					$stemmed_tokens[$count_token]=$token;
				}
				if($POS_tokens[$count_token] eq "")
				{
					$POS_tokens[$count_token]=$token;
				}
				$temmed=$stemmed_tokens[$count_token];
				$POStag=$POS_tokens[$count_token];
				
				#white space
				my $WhiteSpaceFront="WSF:Y";
				my $WhiteSpaceBack="WSB:Y";
				if($LastTokenPosition_hash{$count_token} == $CurrentStart_hash{$count_token}-1)
				{
					$WhiteSpaceFront="WSF:N";
				}
				if($NextTokenPosition_hash{$count_token} == $CurrentLast_hash{$count_token}+1)
				{
					$WhiteSpaceBack="WSB:N";
				}
				
				#Number of Numbers [0-9] 
				my $tmp=$token;
				$tmp=~s/[^0-9]//g;
				my $Num_num="";
				if(length($tmp)>3){$Num_num="N:4+";}else{$Num_num="N:".length($tmp);}
				
				#Number of Uppercase [A-Z]
				$tmp=$token;
				$tmp=~s/[^A-Z]//g;
				my $Num_Uc="";
				if(length($tmp)>3){$Num_Uc="U:4+";}else{$Num_Uc="U:".length($tmp);}
				
				#Number of Lowercase [a-z]
				$tmp=$token;
				$tmp=~s/[^a-z]//g;
				my $Num_Lc="";
				if(length($tmp)>3){$Num_Lc="L:4+";}else{$Num_Lc="L:".length($tmp);}
				
				#Number of ALL char
				$tmp=$token;
				$tmp=~s/[^a-z]//g;
				my $Num_All="";
				if(length($tmp)>3){$Num_All="A:4+";}else{$Num_All="A:".length($tmp);}
				
				#specific character (;:,.->+_)
				$tmp=$token;
				my $SpecificC="";
				if($tmp=~/[\;\:\,\.\-\>\+\_]/){$SpecificC="-SpecificC1-";}
				elsif($tmp=~/[\(\)]/){$SpecificC="-SpecificC2-";}
				elsif($tmp=~/[\{\}]/){$SpecificC="-SpecificC3-";}
				elsif($tmp=~/[\[\]]/){$SpecificC="-SpecificC4-";}
				elsif($tmp=~/[\\\/]/){$SpecificC="-SpecificC5-";}
				else{$SpecificC="__nil__";}
				
				#chemical CHEM
				$tmp=$token;
				my $CHEM="";
				if($tmp=~/(yl|ylidyne|oyl|sulfonyl)$/){$CHEM="-CHEMinlineSuffix-";}
				elsif($tmp=~/(meth|eth|prop|tetracos)/){$CHEM="-CHEMalkaneStem-";}
				elsif($tmp=~/(di|tri|tetra)/){$CHEM="-CHEMsimpleMultiplier-";}
				elsif($tmp=~/(benzen|pyridin|toluen)/){$CHEM="-CHEMtrivialRing-";}
				elsif($tmp=~/(one|ol|carboxylic|amide|ate|acid|ium|ylium|ide|uide|iran|olan|inan|pyrid|acrid|amid|keten|formazan|fydrazin)$/){$CHEM="-CHEMsuffix-";}
				else{$CHEM="__nil__";}
				
				#MentionType
				$tmp=lc($token);
				my $MentionType="";
				
				if($tmp eq "to" && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$CTD_result_hash{$count_token}="CTD_gene";}
				if($tmp=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
				elsif($tmp=~/^(or|and|,)$/ && $last_token=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-2} eq "CTD_gene" && $CTD_result_hash{$count_token+1} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
				elsif($tmp=~/^(or|and|,)$/ && $next_token=~/^(or|and|,)$/ && $CTD_result_hash{$count_token-1} eq "CTD_gene" && $CTD_result_hash{$count_token+2} eq "CTD_gene"){$MentionType="-Type_GeneConjunction-";}
				elsif($tmp=~/^(ytochrome|cytochrome)$/){$MentionType="-Type_cytochrome-";}
				elsif($tmp=~/^target$/){$MentionType="-Type_target-";}
				elsif($tmp=~/^(irradiation|hybrid|fusion|experiment|gst|est|gap|antigen)$/){$MentionType="-Type_ExperimentNoun-";}
				elsif($tmp=~/^(disease|disorder|dystrophy|deficiency|syndrome|dysgenesis|cancer|injury|neoplasm|diabetes|diabete)$/){$MentionType="-Type_Disease-";}
				elsif($tmp=~/(motif|domain|omain|binding|site|region|sequence|frameshift|DNA|RNA)/){$MentionType="-Type_DomainMotif-";}
				elsif($tmp eq "-" && $next_token=~/(motif|domain|omain|binding|site|region|sequence|frameshift|DNA|RNA)/){$MentionType="-Type_DomainMotif-";}
				elsif($tmp=~/^[rmc]$/ && ($next_token eq "DNA" || $next_token eq "RNA") ){$MentionType="-Type_DomainMotif-";}
				elsif($tmp=~/(famil|complex|cluster|proteins|genes|factors|transporter|proteinase|membrane|ligand|enzyme|channels|tors$|ase$|ases$)/){$MentionType="-Type_Family-";}
				elsif($tmp=~/^marker/ && ($last_token eq "M" && $tmp=~/arker/) && ($tmp eq "M" && $next_token=~/arker/)){$MentionType="-Type_Marker-";}
				elsif($tmp=~/(cell)/ || ($next_token=~/(cell)/ && $tmp=~/^(T|B|monocytic|cancer|tumor|myeloma|epithelial|crypt)$/) ){$MentionType="-Type_Cell-";}
				#elsif($tmp eq "chromosome" ){$MentionType="-Type_Chromosome-";}
				#elsif($tmp=~/^[pq]$/ && ($next_token=~/^[0-9]+$/ || $last_token=~/^[0-9]+$/)){$MentionType="-Type_ChromosomeStrain-";}
				elsif($tmp=~/(related|regulated|associated|correlated|reactive)/){$MentionType="-Type_relation-";}
				elsif($lc_tmp=~/(polymorphism|mutation|deletion|insertion|duplication|genotype|genotypes)/){$MentionType="-Type_VariationTerms-";}
				elsif($tmp=~/(oxidase|transferase|transferases|kinase|kinese|subunit|unit|receptor|adrenoceptor|transporter|regulator|transcription|antigen|protein|gene|factor|member|molecule|channel|deaminase|spectrin)/){$MentionType="-Type_suffix-";}
				elsif($tmp=~/[\(\-\_]/ && lc($next_token)=~/^(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])$/){$MentionType="-Type_strain-";}
				elsif($tmp=~/^(alpha|beta|gamma|delta|theta|kappa|zeta|sigma|omega|i|ii|iii|iv|v|vi|[abcdefgyr])$/ && $POS_tokens[$count_token] ne "det"){$MentionType="-Type_strain-";}
				else{$MentionType="__nil__";}
				
				#Protein symbols
				my $uc_tmp=$token;
				my $lc_tmp=lc($token);
				my $ProteinSym="";
				if($lc_tmp=~/(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift)/){$ProteinSym="-ProteinSymFull-";}
				elsif($lc_tmp=~/^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)$/){$ProteinSym="-ProteinSymTri-";}
				elsif($lc_tmp=~/^(ys|le|er|ln|et|sn|ro|ys|sp|hr|he|la|ly|is|eu|rg|rp|al|lu|yr)$/ && $last_token=~/^[CISQMNPKDTFAGHLRWVEYX]$/){$ProteinSym="-ProteinSymTriSub-";}
				elsif($uc_tmp=~/^[CISQMNPKDTFAGHLRWVEYX]$/ && $POS_tokens[$count_token] ne "det"){$ProteinSym="-ProteinSymChar-";}
				else{$ProteinSym="__nil__";}
				
				#Patterns
				my $Pattern1=$token;
				if($Pattern1=~/[\W\-\_]/){$Pattern1="__nil__";}
				else{
					$Pattern1=~s/[A-Z]/A/g;
					$Pattern1=~s/[a-z]/a/g;
					$Pattern1=~s/[0-9]/0/g;
					$Pattern1="P1:".$Pattern1;
				}
				my $Pattern2=$token;
				if($Pattern2=~/[\W\-\_]/){$Pattern2="__nil__";}
				else{
					$Pattern2=~s/[A-Za-z]/a/g;
					$Pattern2=~s/[0-9]/0/g;
					$Pattern2="P2:".$Pattern2;
				}
				my $Pattern3=$token;
				if($Pattern3=~/[\W\-\_]/){$Pattern3="__nil__";}
				else{
					$Pattern3=~s/[A-Z]+/A/g;
					$Pattern3=~s/[a-z]+/a/g;
					$Pattern3=~s/[0-9]+/0/g;
					$Pattern3="P3:".$Pattern3;
				}
				my $Pattern4=$token;
				if($Pattern4=~/[\W\-\_]/){$Pattern4="__nil__";}
				else{
					$Pattern4=~s/[A-Za-z]+/a/g;
					$Pattern4=~s/[0-9]+/0/g;
					$Pattern4="P4:".$Pattern4;
				}
				
				#prefix
				my $prefix="";
				$temp=$token;
				if(length($temp)>=1){ $prefix=$prefix.substr($temp,0,1);}else { $prefix=$prefix."__nil__";}
				if(length($temp)>=2){ $prefix=$prefix." ".substr($temp,0,2);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=3){ $prefix=$prefix." ".substr($temp,0,3);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=4){ $prefix=$prefix." ".substr($temp,0,4);}else { $prefix=$prefix." __nil__";}
				if(length($temp)>=5){ $prefix=$prefix." ".substr($temp,0,5);}else { $prefix=$prefix." __nil__";}
				
				#suffix
				my $suffix="";
				$temp=$token;
				if(length($temp)>=1){ $suffix=$suffix.substr($temp,-1,1);}else { $suffix=$suffix."__nil__";}
				if(length($temp)>=2){ $suffix=$suffix." ".substr($temp,-2,2);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=3){ $suffix=$suffix." ".substr($temp,-3,3);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=4){ $suffix=$suffix." ".substr($temp,-4,4);}else { $suffix=$suffix." __nil__";}
				if(length($temp)>=5){ $suffix=$suffix." ".substr($temp,-5,5);}else { $suffix=$suffix." __nil__";}		
				
				if((not exists $Result_hash{$count_token}) || $Result_hash{$count_token} eq ""){$Result_hash{$count_token}="O";}
				print outputA $token." ".$temmed." ".$POStag." ".$WhiteSpaceFront." ".$WhiteSpaceBack." ".$Num_num." ".$Num_Uc." ".$Num_Lc." ".$Num_All." ".$SpecificC." ".$CHEM." ".$MentionType." ".$ProteinSym." ".$Pattern1." ".$Pattern2." ".$Pattern3." ".$Pattern4." ".$prefix." ".$suffix;
				
				if(exists $PostProessStopWord_hash{lc($token)} && $MentionType eq "__nil__" && $ProteinSym eq "__nil__"&& not exists $CTD_result_hash{$count_token}){$CTD_result_hash{$count_token}="StopWord";}
				elsif(not exists $CTD_result_hash{$count_token}){$CTD_result_hash{$count_token}="__nil__";}
				print outputA " ".$CTD_result_hash{$count_token}; #CTD result

				my $ABBtype="__nil__";
				if(exists $full_result_hash{$count_token}){$ABBtype=$full_result_hash{$count_token};}
				elsif(exists $abb_result_hash{$count_token}){$ABBtype=$abb_result_hash{$count_token};}
				print outputA " ".$ABBtype; #Abb result
				
				if ($TrainorTest eq "Test"){print outputA "\n";}
				else {print outputA " ".$Result_hash{$count_token}."\n";}
				
				$count_token++;
				$last_token=$token;
			}
			if($text ne "")
			{
				print outputA "                                 \n";
			}
			close outputA;
			
			%Result_hash=();
			%locationmap_hash=();
			%lastmap_hash=();
			%Mention2Type_hash=();
			%type_hash=();
			
			%Gene_CTD_locationmap_hash=();
			%Gene_CTD_lastmap_hash=();
			%CTD_result_hash=();
			
			%abb_locationmap_hash=();
			%abb_lastmap_hash=();
			%abb_result_hash=();
			%full_locationmap_hash=();
			%full_lastmap_hash=();
			%full_result_hash=();
		}
	}
	close input;
	
	return 1;
}

sub TranslateOutput
{
	my $sentence_extraction=$_[0];
	my $CRF_output=$_[1];
	my $location=$_[2];
	my $gmention_extraction=$_[3];
	my $threshold=$_[4];
	my $threshold_Gene=$_[5];

	my %article_hash=();
	my %pmidsid_hash=();my $pmid_count=0;
	open input,"<".$sentence_extraction;
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^<TEXT pmid='([^<>]+)' sid='([^<>]+)'>(.+)<\/TEXT>$/)
		{
			$pmid=$1;
			$sid=$2;
			$sentence=$3;
			$sentence=~s/[\n\r]//g;
			$article_hash{$pmid."\t".$sid}=$sentence;
			$pmidsid_hash{$pmid."\t".$sid}=$pmid_count++;
		}
	}
	close input;
	my @pmidsidARR = sort {$pmidsid_hash{$a} <=> $pmidsid_hash{$b}} keys %pmidsid_hash;

	my %output0_hash=();
	my %output1_hash=();
	my %output2_hash=();
	my %score0_hash=();
	my %score1_hash=();
	my %score2_hash=();
	my $count0=0;
	my $count1=0;
	my $count2=0;
	my $rank="";
	my $score=0;
	open output,"<".$CRF_output;
	while(<output>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/# ([0-2]) ([0-9\.]+)/)
		{
			$rank=$1;
			$score=$2;
		}
		elsif($rank eq "0")
		{
			$output0_hash{$count0}=$tmp;
			$score0_hash{$count0}=$score;
			$count0++;
		}
		elsif($rank eq "1")
		{
			$output1_hash{$count1}=$tmp;
			$score1_hash{$count1}=$score;
			$count1++;
		}
		elsif($rank eq "2")
		{
			$output2_hash{$count2}=$tmp;
			$score2_hash{$count2}=$score;
			$count2++;
		}
	}
	
	my %location_hash=();
	$count=0;
	open location,"<".$location;
	while(<location>)
	{
		my $location=$_;
		$location=~s/[\n\r]//g;
		$location_hash{$count}=$location;
		$count++;
	}
	
	my %printSTR0_hash=();
	for(my $i=0;$i<$count;$i++)
	{
		my $output=$output0_hash{$i};
		my $location=$location_hash{$i};
		
		my $start=100000;
		my $last=0;
		my $pmidsid="";
		#Gene G/E/N
		if($output=~/	[G]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output0_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[EN]$/)
			{
				while($output=~/	[EN]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output0_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR0_hash{$pmidsid}=$printSTR0_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Gene	score:".$score0_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#FamilyName F/A/M
		if($output=~/	[F]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output0_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[AM]$/)
			{
				while($output=~/	[AM]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output0_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR0_hash{$pmidsid}=$printSTR0_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	FamilyName	score:".$score0_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#DomainMotif D/I/T
		if($output=~/	[D]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output0_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[IT]$/)
			{
				while($output=~/	[IT]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output0_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR0_hash{$pmidsid}=$printSTR0_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	DomainMotif	score:".$score0_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#Cell C/H/L
		if($output=~/	[C]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output0_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[HL]$/)
			{
				while($output=~/	[HL]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output0_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR0_hash{$pmidsid}=$printSTR0_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Cell	score:".$score0_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#ChromosomeLocation X/Y/Z
		if($output=~/	[X]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output0_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[YZ]$/)
			{
				while($output=~/	[YZ]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output0_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR0_hash{$pmidsid}=$printSTR0_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	ChromosomeLocation	score:".$score0_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}		
	}
	
	my %printSTR1_hash=();
	for(my $i=0;$i<$count;$i++)
	{
		my $output=$output1_hash{$i};
		my $location=$location_hash{$i};
		my $start=100000;
		my $last=0;
		my $pmidsid="";
		#Gene G/E/N
		if($output=~/	[G]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output1_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[EN]$/)
			{
				while($output=~/	[EN]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output1_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR1_hash{$pmidsid}=$printSTR1_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Gene	score:".$score1_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#FamilyName F/A/M
		if($output=~/	[F]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output1_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[AM]$/)
			{
				while($output=~/	[AM]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output1_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR1_hash{$pmidsid}=$printSTR1_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	FamilyName	score:".$score1_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}	
		#DomainMotif D/I/T
		if($output=~/	[D]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output1_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[IT]$/)
			{
				while($output=~/	[IT]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output1_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR1_hash{$pmidsid}=$printSTR1_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	DomainMotif	score:".$score1_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#Cell C/H/L
		if($output=~/	[C]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output1_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[HL]$/)
			{
				while($output=~/	[HL]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output1_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR1_hash{$pmidsid}=$printSTR1_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Cell	score:".$score1_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#ChromosomeLocation X/Y/Z
		if($output=~/	[Z]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output1_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[YZ]$/)
			{
				while($output=~/	[YZ]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output1_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR1_hash{$pmidsid}=$printSTR1_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	ChromosomeLocation	score:".$score1_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
	}

	my %printSTR2_hash=();
	for(my $i=0;$i<$count;$i++)
	{
		my $output=$output2_hash{$i};
		my $location=$location_hash{$i};
		my $start=100000;
		my $last=0;
		my $pmidsid="";
		#Gene G/E/N
		if($output=~/	[G]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output2_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[EN]$/)
			{
				while($output=~/	[EN]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output2_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR2_hash{$pmidsid}=$printSTR2_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Gene	score:".$score2_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#FamilyName F/A/M
		if($output=~/	[F]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output2_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[AM]$/)
			{
				while($output=~/	[AM]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output2_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR2_hash{$pmidsid}=$printSTR2_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	FamilyName	score:".$score2_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}	
		#DomainMotif D/I/T
		if($output=~/	[D]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output2_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[IT]$/)
			{
				while($output=~/	[IT]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output2_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR2_hash{$pmidsid}=$printSTR2_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	DomainMotif	score:".$score2_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#Cell C/H/L
		if($output=~/	[C]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output2_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[HL]$/)
			{
				while($output=~/	[HL]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output2_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR2_hash{$pmidsid}=$printSTR2_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	Cell	score:".$score2_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
		#ChromosomeLocation X/Y/Z
		if($output=~/	[X]$/)
		{
			if($location=~/^([^\t]+	[^\t]+)	[^\t]+	([0-9]+)	([0-9]+)/)
			{
				$pmidsid=$1;
				if($2<$start){$start=$2;}
				if($3>$last){$last=$3;}
			}			
			$i++;
			$output=$output2_hash{$i};
			$location=$location_hash{$i};
			if($output=~/	[YZ]$/)
			{
				while($output=~/	[YZ]$/)
				{
					if($location=~/^[^\t]+	[^\t]+	[^\t]+	([0-9]+)	([0-9]+)/)
					{
						if($1<$start){$start=$1;}
						if($2>$last){$last=$2;}
					}
					$i++;
					$output=$output2_hash{$i};
					$location=$location_hash{$i};
				}
			}
			$printSTR2_hash{$pmidsid}=$printSTR2_hash{$pmidsid}.$pmidsid."	".($start-1)."	".$last."	".substr($article_hash{$pmidsid},$start-1,$last-($start-1))."	ChromosomeLocation	score:".$score2_hash{$i}."\n";
			if($output!~/	O$/)
			{
				$i--;
			}
		}
	}
	close output;
	close location;
	
	my %printSTR_hash=();
	foreach my $pmidsid (@pmidsidARR)
	{
		my %Ann0_start_last=();
		my @Ann = split /\n/, $printSTR0_hash{$pmidsid};
		foreach my $Ann(@Ann) 
		{
			my ($i,$j)= split(/\tscore:/, $Ann);
			if($i=~/^([^\t]+	[^\t]+)	([0-9]+)	([0-9]+)	([^\t]+)	([^\t]+)/)
			{
				my $pmidsid=$1;
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $type=$5;
				$printSTR_hash{$pmidsid} = $printSTR_hash{$pmidsid}.$i."\t".$j."\n";
				$Ann0_start_last{$i}=$type;
			}
			
		}
			
		@Ann = split /\n/, $printSTR1_hash{$pmidsid};
		foreach my $Ann(@Ann) 
		{
			my ($i,$j)= split(/\tscore:/, $Ann);
			if($j>$threshold && $i=~/^([^\t]+	[^\t]+)	([0-9]+)	([0-9]+)	([^\t]+)	([^\t]+)/)
			{
				my $pmidsid=$1;
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $type=$5;
				$mention=~s/[\W\-\_]/\./g;
				my $overlap=0;
				foreach my $STR(keys %Ann0_start_last)
				{
					my $type0=$Ann0_start_last{$STR};
					if($STR=~/^$pmidsid	([0-9]+)	([0-9]+)/)
					{
						my $start0=$1;
						my $last0=$2;
						if($start==$start0 && $last==$last0 && $type0 ne "Gene" && $type eq "Gene" && $j>$threshold_Gene)
						{
							$printSTR_hash{$pmidsid}=~s/($pmidsid\t$start\t$last\t$mention\t)$type0/$1$type/g;
						}
						elsif( ($start>=$start0 && $start<$last0) || ($last>$start0 && $last<=$last0) )
						{
							$overlap=1;#overlap
						}
					}
				}
				if($overlap == 0)
				{
					$printSTR_hash{$pmidsid} = $printSTR_hash{$pmidsid}.$i."\t".$j."\n";
				}
				$Ann0_start_last{$i}=$type;
			}
		}
		
		@Ann = split /\n/, $printSTR2_hash{$pmidsid};
		foreach my $Ann(@Ann) 
		{
			my ($i,$j)= split(/\tscore:/, $Ann);
			if($j>$threshold && $i=~/^([^\t]+	[^\t]+)	([0-9]+)	([0-9]+)	([^\t]+)	([^\t]+)/)
			{
				my $pmidsid=$1;
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $type=$5;
				$mention=~s/[\W\-\_]/\./g;
				my $overlap=0;
				foreach my $STR(keys %Ann0_start_last)
				{
					my $type0=$Ann0_start_last{$STR};
					if($STR=~/^$pmidsid	([0-9]+)	([0-9]+)/)
					{
						my $start0=$1;
						my $last0=$2;
						if($start==$start0 && $last==$last0 && $type0 ne "Gene" && $type eq "Gene" && $j>$threshold_Gene)
						{
							$printSTR_hash{$pmidsid}=~s/($pmidsid\t$start\t$last\t$mention\t)$type0/$1$type/g;
						}
						elsif( ($start>=$start0 && $start<$last0) || ($last>$start0 && $last<=$last0) )
						{
							$overlap=1;#overlap
						}
					}
				}
				if($overlap == 0)
				{
					$printSTR_hash{$pmidsid} = $printSTR_hash{$pmidsid}.$i."\t".$j."\n";
				}
			}
		}
	}
	
	open output,">".$gmention_extraction;
	foreach my $pmidsid (@pmidsidARR)
	{
		foreach my $line(split("\n",$printSTR_hash{$pmidsid}))
		{
			if($line=~/^([^\t]+)	([^\t]+)	([0-9]+)	([0-9]+)	([^\t]+)	([^\t]+)	([^\t]+)/)
			{
				my $pmid=$1;
				my $sid=$2;
				my $start=$3;
				my $last=$4;
				my $mention=$5;
				my $type=$6;
				my $score=$7;
				print output "<$type pmid='$pmid' sid='$sid' start='$start' end='$last' score='$score'>$mention</$type>\n";
			}
			elsif($line=~/^([^\t]+)	([^\t]+)	([0-9]+)	([0-9]+)	([^\t]+)	([^\t]+)/)
			{
				my $pmid=$1;
				my $sid=$2;
				my $start=$3;
				my $last=$4;
				my $mention=$5;
				my $type=$6;
				print output "<$type pmid='$pmid' sid='$sid' start='$start' end='$last' score='0.001'>$mention</$type>\n";
			}
		}
	}
	close output;
	
	#open clear,"|rm tmp/*.sentence.xml.data";	close clear;
	#open clear,"|rm tmp/*.sentence.xml.output";	close clear;
	#open clear,"|rm tmp/*.sentence.xml.location";	close clear;

	return 1;
}

sub PostProcessing
{
	my $filename=$_[0];
	my $sentence_extraction=$_[1];
	my $gmention_extraction=$_[2];
	my $tax_extraction=$_[3];
	my $threshold_Voting=$_[4];
	
	my %Family2Num_hash=(); #pmid&mention -> Num
	my %Domain2Num_hash=(); #pmid&mention -> Num
	my %Cell2Num_hash=(); #pmid&mention -> Num
	my %Mention2Type_hash=(); #pmid&mention -> type
	my %Mention2Num_hash=(); #pmid&mention -> Num
	my %Mention2Score_hash=(); #pmid&mention -> Score
	
	my %article_hash=(); #pmid&sid -> sentence
	my %pmidsid_hash=(); #temporary : For generating @pmidsidARR
	my %pmid_hash=(); #temporary : For generating @pmidsidARR
	my %gene_substr_hash=(); #sub-string of a gene mention (used to infer family name: "TIF" of TIF1)
	
	#Read stop-words
	my %PostProessStopWord_hash=();
	open PostProessStopWord,"<./Library/Regex/PostProessStopWord.txt";	
	while(<PostProessStopWord>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$PostProessStopWord_hash{$tmp}=length($tmp);
	}
	close PostProessStopWord;
	open tax_extraction,"<$tax_extraction";	
	while(<tax_extraction>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/>(.+)</)
		{
			$PostProessStopWord_hash{$1}=length($1);
		}
	}
	close tax_extraction;
	my @PostProessStopWordArr = reverse sort {$PostProessStopWord_hash{$a} <=> $PostProessStopWord_hash{$b}} keys %PostProessStopWord_hash;
	my $PostProessStopWordRegEx = Regex::PreSuf::presuf(@PostProessStopWordArr);
	
	#Read sentence
	my $pmidsid_count=0;
	open input,"<".$sentence_extraction;
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^<TEXT pmid='([^<>]+)' sid='([^<>]+)'>(.+)<\/TEXT>$/)
		{
			$pmid=$1;
			$sid=$2;
			$sentence=$3;
			$sentence=~s/[\n\r]//g;
			$article_hash{$pmid."\t".$sid}=$sentence;
			$pmid_hash{$pmid}=$pmidsid_count;
			$pmidsid_hash{$pmid."\t".$sid}=$pmidsid_count++;
		}
	}
	close input;
	my @pmidsidARR = sort {$pmidsid_hash{$a} <=> $pmidsid_hash{$b}} keys %pmidsid_hash;
	my @pmidARR = sort {$pmid_hash{$a} <=> $pmid_hash{$b}} keys %pmid_hash;
	
	#Mention2Type_hash/Family2Num_hash/Domain2Num_hash/Cell2Num_hash establishment
	{
		my @tmp_arr;
		open input,"<".$gmention_extraction;
		while(<input>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<(.+) pmid='(.+)' sid='(.+)' start='(.+)' end='(.+)' score='(.+)'>(.+)<\/.+>/)
			{
				my $type=$1;
				my $pmid=$2;
				my $sid=$3;
				my $start=$4;
				my $last=$5;
				my $score=$6;
				my $mention=$7;
				
				#Check suffix
				if($mention=~/(cell|cells)$/i){$type="Cell";}
				elsif($mention=~/(disease|diseases|syndrome|syndromes|tumor|tumour|deficiency|dysgenesis|dystrophy|frame|factors|family|families|complex|genes|proteins)$/i){$type="FamilyName";}
				elsif($mention=~/(domain|motif|domains|motifs|sequences)$/i){$type="DomainMotif";}
				
				my $lcmention=lc($mention);
				$lcmention=~s/[\n\r\W\-\_]/ /g;
				
				#filter the StopWords
				# [0-9]+[qp][0-9\.]+	:	chromosome
				# [0-9,]+[km]{0,1}b|[0-9,\.]+	chromosome
				# table.*|fig.*	Table/Figure
				if( (length($lcmention)>=3 && $lcmention!~/^($PostProessStopWordRegEx|[0-9]+[qp][0-9.]+|[0-9,]+[km]{0,1}b|[0-9,.]+|table.*|fig.*|.*virus|.* [mp]ab)$/i) ||
					(length($lcmention)>1 && $lcmention!~/^($PostProessStopWordRegEx|[0-9]+[qp][0-9.]+|[0-9,]+[km]{0,1}b|[0-9,.]+|table.*|fig.*|.*virus|.* [mp]ab)$/i && $score>=0.5)		
					)
				{
					push(@tmp_arr,$type."\t".$pmid."\t".$sid."\t".$start."\t".$last."\t".$score."\t".$mention);
					if($type eq "Gene")
					{
						$gene_substr_hash{$pmid."\t".$mention}=0;
						#sub-string : FamilyName
						foreach my $pmid_compared_mention (keys %gene_substr_hash)
						{
							if($pmid_compared_mention=~/^$pmid	(.+)$/)
							{
								my $compared_mention=$1;
								my $C_compared_mention=lc($compared_mention);
								$C_compared_mention=~s/[\W\-\_]/\./g;
								if($mention=~/^$C_compared_mention[\W\-\_]{0,1}([^\W\-\_]+)/i && lc($1)=~/^(alpha|beta|gamma|kappa|theta|delta|[A-Za-z0-9])$/ )
								{
									$gene_substr_hash{$pmid."\t".lc($compared_mention)}=1;
									delete $gene_substr_hash{$pmid."\t".lc($mention)};
								}
								my $C_mention=lc($mention);
								$C_mention=~s/[\W\-\_]/\./g;
								if($compared_mention=~/^$C_mention[\W\-\_]{0,1}([^\W\-\_]+)/i && lc($1)=~/^(alpha|beta|gamma|kappa|theta|delta|[A-Za-z0-9])$/ )
								{
									$gene_substr_hash{$pmid."\t".lc($mention)}=1;
									delete $gene_substr_hash{lc($pmid."\t".$compared_mention)};
								}
							}
						}
					}
				}
			}
		}
		close input;
		foreach my $pmid(@pmidARR)
		{
			foreach my $t (@tmp_arr)
			{
				if($t=~/^([^\t]+)\t($pmid)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
				{
					my $type=$1;
					my $pmid=$2;
					my $sid=$3;
					my $start=$4;
					my $last=$5;
					my $score=$6;
					my $mention=lc($7);
					
					$Mention2Num_hash{$pmid."\t".$mention}++;
					$Mention2Score_hash{$pmid."\t".$mention}=$Mention2Score_hash{$pmid."\t".$mention}+$score;
					
					my $mention_tmp=$mention;
					if($mention_tmp!~/^([mrc]{0,1}[DR]NA|gene|protein|mouse|mice|human|met|asp|gly|cys|ala|phe|ile||leu|gln|arg|ser|thr|val|trp|xaa|tyr)$/i)
					{
						if($type eq "Gene" && $Mention2Type_hash{$pmid."\t".$mention}=~/(FamilyName|DomainMotif|Cell)/){}
						else{$Mention2Type_hash{$pmid."\t".$mention}=$type;}
						
						if($type eq "FamilyName")
						{
							$Family2Num_hash{$pmid."\t".$mention}++;
						}
						elsif($type eq "DomainMotif")
						{
							$Domain2Num_hash{$pmid."\t".$mention}++;
						}
						elsif($type eq "Cell")
						{
							$Cell2Num_hash{$pmid."\t".$mention}++;
						}
						elsif($gene_substr_hash{$pmid."\t".$mention} == 1) #Assign FamilyName if it is a sub-string of a gene
						{
							$Family2Num_hash{$pmid."\t".$mention}++;
							$Mention2Type_hash{$pmid."\t".$mention}="FamilyName";
						}
						
						#Check the surrounding words
						$sentence=$article_hash{$pmid."\t".$sid};
						if(substr($sentence,$last+1,15)=~/(superfamily|superfamilies|subfamily|subfamilies|family|families|sequences)/ && substr($sentence,$start-1,1) ne "/")
						{
							$Family2Num_hash{$pmid."\t".$mention}=100;
							$Mention2Type_hash{$pmid."\t".$mention}="FamilyName";
						}
						elsif(substr($sentence,$last+1,7)=~/domain/)
						{
							$Domain2Num_hash{$pmid."\t".$mention}=100;
							$Mention2Type_hash{$pmid."\t".$mention}="DomainMotif";
						}
						elsif(substr($sentence,$last+1,7)=~/cell[s]{0,1}[\W\-\_]/)
						{
							$Cell2Num_hash{$pmid."\t".$mention}=100;
							$Mention2Type_hash{$pmid."\t".$mention}="Cell";
						}
					}
				}
			}
		}
	}
	
	#Infer Gene type to family/Domain/Cell mentions based on frequency
	my %AbbreviationMapping_hash=();
	my %pattern_hash=(); #For Pattern Extension 1
	{
		open Abb,"<tmp/".$filename.".Ab3P";
		my $pmid="";
		while(<Abb>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/^([^\|]+)$/)
			{
				$pmid=$1;
				$tmp=<Abb>; #sentences				
			}
			elsif($tmp=~/^  (.+)\|(.+)\|([0-9\.]+)/)
			{
				my $Abb=lc($1);
				my $FN=lc($2);
				$AbbreviationMapping_hash{$pmid."\t".$Abb}=$pmid."\t".$FN;
				$AbbreviationMapping_hash{$pmid."\t".$FN}=$pmid."\t".$Abb;
			}
		}
		close Abb;
		
		foreach my $pmid_mention (keys %Mention2Type_hash)
		{
			#Pattern establishment
			if($Mention2Type_hash{$pmid_mention} eq "Gene" && $pmid_mention=~/^(.+)	(.+)/)
			{
				my $pmid=$1;
				my $mention=$2;
				my $pattern=$mention;
				$pattern=~s/[\W\-\_]/\[\\W\\\-\\\_\]/g;
				$pattern=~s/[0-9]/\[0\-9\]/g;
				$pattern=~s/[IV]/\[IV\]/g;
				$pattern=~s/[A-Z]$/\[A\-Z\]/g;
				$pattern=~s/[a-z]$/\[a\-z\]/g;
				$pattern_hash{$pmid."\t".$pattern}=1;
			}
			#FamilyName -> Gene by Threshold
			elsif($Mention2Type_hash{$pmid_mention}=~/(FamilyName|DomainMotif|Cell)/ && $pmid_mention=~/^(.+)	(.+)/) #FamilyName -> Gene
			{
				my $pmid=$1;
				my $FDC_Num=$Family2Num_hash{$pmid_mention}+
							$Domain2Num_hash{$pmid_mention}+
							$Cell2Num_hash{$pmid_mention}+
							$Family2Num_hash{$AbbreviationMapping_hash{$pmid_mention}}+
							$Domain2Num_hash{$AbbreviationMapping_hash{$pmid_mention}};
				
				my $All_Num=$Mention2Num_hash{$pmid_mention}; #All: GFDC(Gene,Family,Domain,Cell)
				if( ($FDC_Num/$All_Num) <= $threshold_Voting ) #Threshold(+);Recall(+)
				{
					delete $Family2Num_hash{$pmid_mention};
					delete $Domain2Num_hash{$pmid_mention};
					delete $Cell2Num_hash{$pmid_mention};
					$Mention2Type_hash{$pmid_mention}="Gene";
				}
			}
		}
	}
	
	#Pattern Extension 1: extend to exist mentions
	foreach my $pmid_mention (keys %Mention2Type_hash)
	{
		foreach my $pmid_pattern (keys %pattern_hash)
		{
			if($Mention2Type_hash{$pmid_mention} ne "Gene" && $pmid_mention=~/^$pmid_pattern$/)
			{
				$Mention2Type_hash{$pmid_mention}="Gene";
				delete $Family2Num_hash{$pmid_mention};
				delete $Domain2Num_hash{$pmid_mention};
				delete $Cell2Num_hash{$pmid_mention};
			}
		}
	}
	
	#"Abbreviation" mention type follows "Full Name" type
	my $pmid="";
	my %Abbreivation_GeneType_hash=();
	my $AbbThreshold=0.9;
	open Abb,"<tmp/".$filename.".Ab3P";
	while(<Abb>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\|]+)$/)
		{
			$pmid=$1;
			$tmp=<Abb>; #sentences			
		}
		elsif($tmp=~/^  (.+)\|(.+)\|([0-9\.]+)/)
		{
			my $Abb=lc($1);
			my $FN=lc($2);
			my $weight=$3;
			
			if($FN=~/(disease|diseases|syndrome|syndromes|tumor|tumour|deficiency|dysgenesis|atrophy|frame|dystrophy)$/i)
			{
				delete $Mention2Type_hash{$pmid."\t".$FN};
				delete $Mention2Type_hash{$pmid."\t".$Abb};
			}
			elsif($FN=~/(cell|cells)$/i)
			{
				$Mention2Type_hash{$pmid."\t".$FN}="Cell";
				if($weight>=$AbbThreshold){$Mention2Type_hash{$pmid."\t".$Abb}="Cell";}
			}
			elsif($FN=~/(disease|diseases|syndrome|syndromes|tumor|tumour|deficiency|dysgenesis|frame|dystrophy|factors|family|families|complex|genes|proteins)$/i)
			{
				$Mention2Type_hash{$pmid."\t".$FN}="FamilyName";
				if($weight>=$AbbThreshold){$Mention2Type_hash{$pmid."\t".$Abb}="FamilyName";}
			}
			elsif($FN=~/(domain|motif|domains|motifs)$/i)
			{
				$Mention2Type_hash{$pmid."\t".$FN}="DomainMotif";
				if($weight>=$AbbThreshold){$Mention2Type_hash{$pmid."\t".$Abb}="DomainMotif";}
			}
			else
			{
				if((exists $Mention2Type_hash{$pmid."\t".$Abb}) && (exists $Mention2Type_hash{$pmid."\t".$FN}))
				{
					if($weight>=$AbbThreshold)
					{
						$Mention2Type_hash{$pmid."\t".$Abb}=$Mention2Type_hash{$pmid."\t".$FN};
						$Abbreivation_GeneType_hash{$pmid."\t".$FN}=$Mention2Type_hash{$pmid."\t".$FN};
						$Abbreivation_GeneType_hash{$pmid."\t".$Abb}=$Abbreivation_GeneType_hash{$pmid."\t".$FN};
						if($Mention2Type_hash{$pmid."\t".$FN} eq "FamilyName"){$Family2Num_hash{$pmid."\t".$Abb}=1;}
						elsif($Mention2Type_hash{$pmid."\t".$FN} eq "DomainMotif"){$Domain2Num_hash{$pmid."\t".$Abb}=1;}
						elsif($Mention2Type_hash{$pmid."\t".$FN} eq "Cell"){$Cell2Num_hash{$pmid."\t".$Abb}=1;}
					}
					elsif($Mention2Type_hash{$pmid."\t".$FN} eq "Gene")
					{
						$Abbreivation_GeneType_hash{$pmid."\t".$Abb}=$Abbreivation_GeneType_hash{$pmid."\t".$FN};
						delete $Family2Num_hash{$pmid."\t".$Abb};
						delete $Domain2Num_hash{$pmid."\t".$Abb};
						delete $Cell2Num_hash{$pmid."\t".$Abb};
					}
				}
				elsif((exists $Mention2Type_hash{$pmid."\t".$Abb}) && (not exists $Mention2Type_hash{$pmid."\t".$FN}))
				{
					if($weight>=$AbbThreshold){ $Mention2Type_hash{$pmid."\t".$FN}=$Mention2Type_hash{$pmid."\t".$Abb}; }
				}
				elsif((not exists $Mention2Type_hash{$pmid."\t".$Abb}) && (exists $Mention2Type_hash{$pmid."\t".$FN}))
				{
					if($weight>=$AbbThreshold)
					{
						$Mention2Type_hash{$pmid."\t".$Abb}=$Mention2Type_hash{$pmid."\t".$FN};
						if($Mention2Type_hash{$pmid."\t".$FN} eq "Gene")
						{
							$Abbreivation_GeneType_hash{$pmid."\t".$Abb}="Gene";
							$Abbreivation_GeneType_hash{$pmid."\t".$FN}="Gene";
						}
					}
				}
			}
		}
	}
	close Abb;
	
	#Pattern Extension 2: Pattern match in Text, put in %Mention2Type_hash
	{
		my %pattern_text_hash=();
		foreach my $pmid_mention (keys %Mention2Type_hash)
		{
			#Pattern establishment
			if($Mention2Type_hash{$pmid_mention} eq "Gene" && $pmid_mention=~/^(.+)	(.+)/)
			{
				my $pmid=$1;
				my $mention=$2;
				my $pattern=$mention;
				$pattern=~s/[\W\-\_]/\[\\W\\\-\\\_\]/g;
				$pattern=~s/[a-c]$/\[a\-c\]/g;
				$pattern=~s/[0-9]/\[0\-9\]/g;
				$pattern_text_hash{$pmid."\t".$pattern}=length($mention);
			}
		}
		foreach my $pmid_pattern (keys %pattern_text_hash)
		{
			my $len=$pattern_text_hash{$pmid_pattern};
			if($len>3 && $pmid_pattern=~/^(.+)	(.+)/)
			{
				my $pmid=$1;
				my $RegEx=$2;
				foreach my $pmidsid(@pmidsidARR)
				{
					if($pmidsid=~/$pmid	/)
					{
						my $article_tmp=" ".$article_hash{$pmidsid}." ";
						while($article_tmp =~ /^(.*[\W\-\_])($RegEx)([\W\-\_].*)$/i)
						{
							my $str_1=$1;
							my $str_2=$2;
							my $mention=$2;
							my $str_3=$3;
							$str_2=~ s/./@/g;
							$article_tmp = $str_1.$str_2.$str_3;
							if(not exists $Mention2Type_hash{$pmid."\t".lc($mention)})
							{
								$Mention2Type_hash{$pmid."\t".lc($mention)}="Gene";
								$Mention2Num_hash{$pmid."\t".lc($mention)}=$Mention2Num_hash{$pmid."\t".lc($mention)}+0.5;
								if($Mention2Num_hash{$pmid_mention}==0){$Mention2Num_hash{$pmid_mention}=1;}
								$Mention2Score_hash{$pmid."\t".lc($mention)}=($Mention2Score_hash{$pmid_mention}/$Mention2Num_hash{$pmid_mention})*0.5;
							}
						}
					}
				}
			}
		}
	}
	
	#Extend identified mention from %Mention2Type_hash to Text
	my %All_Num_hash=();
	foreach my $pmid (@pmidARR)
	{	
		my $article=" ";
		foreach my $pmidsid(keys %article_hash)
		{
			if($pmidsid=~/^$pmid	/)
			{
				$article=$article.$article_hash{$pmidsid}." ";
			}
		}
		#calculate %All_Num_hash
		foreach my $pmid_men(keys %Mention2Type_hash)
		{
			if($pmid_men=~/^$pmid	(.+)$/)
			{
				my $article_tmp=$article;
				my $men=$1;
				$men=~s/[\W\-\_]/\./g;
				while( $article_tmp =~ /^(.*?[\W\-\_])($men)([\W\-\_].*)$/i)
				{
					$article_tmp = $1.$3;
					$All_Num_hash{$pmid_men}++;
				}
			}
		}
	}
	foreach my $pmid_mention (keys %Mention2Type_hash)
	{
		if(exists $AbbreviationMapping_hash{$pmid_mention})
		{
			$Mention2Score_hash{$pmid_mention}=$Mention2Score_hash{$pmid_mention}+$Mention2Score_hash{$AbbreviationMapping_hash{$pmid_mention}};
			$Mention2Score_hash{$AbbreviationMapping_hash{$pmid_mention}}=$Mention2Score_hash{$pmid_mention};
			$All_Num_hash{$pmid_mention}=$All_Num_hash{$pmid_mention}+$All_Num_hash{$AbbreviationMapping_hash{$pmid_mention}};
			$All_Num_hash{$AbbreviationMapping_hash{$pmid_mention}}=$All_Num_hash{$pmid_mention};
		}
	}
	
	my %annotation_hash=();
	my %annotation_sort_hash=();
	my $count=0;
	#open score,">score.txt";
	foreach my $pmidsid (@pmidsidARR)
	{
		my %tmp_Mention_hash=();
		my %tmp_Type_hash=();
		my %tmp_Sort_hash=();
		my ($pmid,$sid)=($pmidsid=~/^([^\t]+)\t([^\t]+)$/);
		my $article=" ".$article_hash{$pmidsid}." ";
		
		foreach my $pmid_men(keys %Mention2Type_hash)
		{
			if($pmid_men=~/^$pmid	(.+)$/)
			{
				$men=$1;
				if(length($men)>2)
				{
					my $men_org=$pmid."	".$men;						
					$men=lc($men);
					$men=~s/[\W\-\_]/\./g;

					if($Mention2Num_hash{$pmid_men} == 0){$Mention2Num_hash{$pmid_men}=1;}
					if($All_Num_hash{$pmid_men} == 0) {}
					elsif(($Mention2Score_hash{$pmid_men}/$All_Num_hash{$pmid_men})>=$threshold_Gene)
					{
						print score $pmid_men."\t".$Mention2Score_hash{$pmid_men}."\t".$All_Num_hash{$pmid_men}."\n";
						$article_tmp=$article;
						if($men!~/^\.+$/)
						{
							while( $article_tmp =~ /^(.*?[\W\-\_])($men)([\W\-\_].*)$/i)
							{
								my $str_1=$1;
								my $str_2=$2;
								my $mention=$2;
								my $str_3=$3;
								$str_2=~ s/./@/g;
								$article_tmp = $str_1.$str_2.$str_3;
								my $start=length($str_1)-1;
								my $last =length($str_1)+length($str_2)-1;
								$tmp_Mention_hash{$pmidsid."\t".$start."\t".$last}=$mention;
								$tmp_Type_hash{$pmidsid."\t".$start."\t".$last}=$Mention2Type_hash{$pmid_men};
								$tmp_Sort_hash{$pmidsid."\t".$start."\t".$last}=$start;
							}
						}
					}
				}
			}
		}
		foreach my $Ext1(keys %tmp_Mention_hash)
		{
			my ($pmid1,$sid1,$start1,$last1)=split("\t",$Ext1);
			my $type1=$tmp_Type_hash{$Ext1};
			my $mention1=$tmp_Mention_hash{$Ext1};
			foreach my $Ext2(keys %tmp_Mention_hash)
			{
				my ($pmid2,$sid2,$start2,$last2)=split("\t",$Ext2);
				my $type2=$tmp_Type_hash{$Ext2};
				my $mention2=$tmp_Mention_hash{$Ext2};
				
				if($pmid1 eq $pmid2 && $sid1 eq $sid2 && $Ext1 ne $Ext2)
				{
					if(($start1<$start2 && $last1>=$last2)||($start1<=$start2 && $last1>$last2))
					{
						delete $tmp_Mention_hash{$Ext2};
					}
				}
			}
		}
		foreach my $Ext(keys %tmp_Mention_hash)
		{
			my $mention=$tmp_Mention_hash{$Ext};
			my $type=$tmp_Type_hash{$Ext};
			my ($pmid,$sid,$start,$last)=split("\t",$Ext);
			#not display Cell/Chromosome
			if($tmp_Type_hash{$Ext} ne "Cell" && $tmp_Type_hash{$Ext} ne "Chromosone")
			{
				if((exists $Mention2Score_hash{$pmid."\t".lc($mention)}) && (exists $All_Num_hash{$pmid."\t".lc($mention)}))
				{				
					my $weight = round( $Mention2Score_hash{$pmid."\t".lc($mention)}/$All_Num_hash{$pmid."\t".lc($mention)} , 4 );
					$annotation_hash{$pmid."\t".$sid."\t".$Ext}="<$type pmid='$pmid' sid='$sid' start='$start' end='$last' score='$weight'>$mention</$type>\n";
					$annotation_sort_hash{$pmid."\t".$sid."\t".$Ext}=$start;
				}
			}
		}
	}
	#close score;

	my %output_hash=();
	my @annotation_sort = sort {$annotation_sort_hash{$a} <=> $annotation_sort_hash{$b}} keys %annotation_sort_hash;
	foreach my $tmp(@annotation_sort)
	{
		my ($pmid,$sid,$start,$last)=split(/\t/,$tmp);
		$output_hash{$pmid."\t".$sid}=$output_hash{$pmid."\t".$sid}.$annotation_hash{$tmp};
	}
	
	#Extending all family/domain names
	foreach my $pmidsid (@pmidsidARR)
	{
		foreach my $T(keys %Family2Num_hash)
		{
			if(not exists $Abbreivation_GeneType_hash{$T})
			{
				my ($pmid,$family)=split(/\t/,$T);
				$family=~s/[\W\-\_]/\./g;
				$output_hash{$pmidsid}=~s/\t($family)\t(Gene|Cell)/\t$1\tFamilyName/ig;
			}
		}
		foreach my $T(keys %Domain2Num_hash)
		{
			if(not exists $Abbreivation_GeneType_hash{$T})
			{
				my ($pmid,$domain)=split(/\t/,$T);
				$domain=~s/[\W\-\_]/\./g;
				$output_hash{$pmidsid}=~s/\t($domain)\t(Gene|Cell)/\t$1\tDomainMotif/ig;
			}
		}
	}
	
	open output,">".$gmention_extraction;
	foreach my $pmidsid (@pmidsidARR)
	{
		if($output_hash{$pmidsid} ne "")
		{
			print output $output_hash{$pmidsid};
		}
	}
	close output;	

	return 1;
}

sub Gene_Name_Recognition
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	my $TrainorTest="Test";
	
	my $threshold=0.005;
	my $threshold_Gene=0.09;
	my $threshold_Voting=0.5; # (Family+Domain+Cell)/All
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $gmention_extraction="tmp/".$filename.".gmention.xml";
	my $tax_extraction="tmp/".$filename.".tax.xml";

	FeatureExtraction($filename,$sentence_extraction,$TrainorTest);
	if($^O=~/Win/)
	{
		open(CRF,"|CRFmodel/crf_test -n 3 -m CRFmodel/Model.GMR -o ".$sentence_extraction.".output ".$sentence_extraction.".data");
		close CRF;
	}
	else
	{
		open(CRF,"|./CRFmodel/crf_test -n 3 -m CRFmodel/Model.GMR -o ".$sentence_extraction.".output ".$sentence_extraction.".data")|| die ("An Error in CRF module. \nPlease re-install CRF module:\n\n\tsh Installation.sh\n");	close(CRF);
		close CRF;
	}

	TranslateOutput($sentence_extraction,$sentence_extraction.".output",$sentence_extraction.".location",$gmention_extraction,$threshold,$threshold_Gene);
	PostProcessing($filename,$sentence_extraction,$gmention_extraction,$tax_extraction,$threshold_Voting);
	
	return 1;
}

sub ID_Extraction
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $id_extraction="tmp/".$filename.".id.xml";
	
	open sen,"<$sentence_extraction";
	open output,">$id_extraction";
	while(<sen>)
	{
		my $tmp=$_;
		$tmp=~s/[\r\n]//g;
		my @array_tagged;
		if($tmp=~/<TEXT pmid='(.+)' sid='(.+)'>(.+)<\/TEXT>/)
		{
			my $pmid=$1;
			my $sid=$2;
			my $sentence=$3;
			my $sentence_org=$3;
			$sentence=lc($sentence);
			$sentence=~ s/\n/ /g;
			
			$sentence =~s/[,\(\)\{\}\[\]]/ /g;
			$sentence =~s/\_/ZZZ/g;
			#Extraction
			{
				if ( $sentence =~ /^(\S*[0-9]+\S*[A-Za-z]+\S*)([^0-9A-Za-z]+.*)$/ ||  
					 $sentence =~ /^(\S*[A-Za-z]+\S*[0-9]+\S*)([^0-9A-Za-z]+.*)$/)
				{
					my $str1="";
					my $str_temp=$1;
					my $str3=$2;
					my $len_str1=$str1;	$len_str1=~s/ZZZ/\_/g;
					my $len_str_temp=$str_temp;	$len_str_temp=~s/ZZZ/\_/g;
					my $str2=substr($sentence_org,length($len_str1),length($len_str_temp));
					$str_temp=~ s/ZZZ/\@\@\@/g;
					$str_temp=~ s/[^@]/A/g;
					$str_temp=~ s/\@\@\@/ZZZ/g;
					$sentence = $str1.$str_temp.$str3;
					push (@array_tagged, "$str1$str2");
				}		
				while ( $sentence =~ /^(.*[^0-9A-Za-z]+)(\S*[0-9]+\S*[A-Za-z]+\S*)([^0-9A-Za-z]+.*)$/ || 
						$sentence =~ /^(.*[^0-9A-Za-z]+)(\S*[0-9]+\S*[A-Za-z]+\S*)$/ ||  
						$sentence =~ /^(.*[^0-9A-Za-z]+)(\S*[A-Za-z]+\S*[0-9]+\S*)([^0-9A-Za-z]+.*)$/ || 
						$sentence =~ /^(.*[^0-9A-Za-z]+)(\S*[A-Za-z]+\S*[0-9]+\S*)$/ )
				{
					my $str1=$1;
					my $str_temp=$2;
					my $str3=$3;
					my $len_str1=$str1;	$len_str1=~s/ZZZ/\_/g;
					my $len_str_temp=$str_temp;	$len_str_temp=~s/ZZZ/\_/g;
					my $str2=substr($sentence_org,length($len_str1),length($len_str_temp));
					$str_temp=~ s/ZZZ/\@\@\@/g;
					$str_temp=~ s/[^@]/A/g;
					$str_temp=~ s/\@\@\@/ZZZ/g;
					$sentence = $str1.$str_temp.$str3;
					push (@array_tagged, "$str1$str2");
				}
			}
			
			foreach $tagged(@array_tagged)
			{
				$tagged=~s/\_/ZZZ/g;
				if($tagged =~ /^(.*[^0-9A-Za-z]+)(\S*[0-9]+\S*[A-Za-z]+\S*)$/ || 
				   $tagged =~ /^(.*[^0-9A-Za-z]+)(\S*[A-Za-z]+\S*[0-9]+\S*)$/)
				{
					$str_1=$1;
					$str_2=$2;	
					$str_2 =~s/ZZZ/\_/g;
					$start_site=length($str_1);
					$end_site =length($str_1)+length($str_2);
					$str_2 =~s/\.//g;
					
					if($str_2 !~/[;\$\#\&\/\"\',\+\(\)\{\}\[\]]/)
					{
						if($str_2 =~ /^(.+[^0-9])([0-9]+)\-([0-9]+)$/)#numeration
						{
							$tmp1=$1;
							$num1=$2;
							$num2=$3;
							my $prefix="";
							if($num1=~/^([0]+)/)
							{
								$prefix=$1;
							}
							if(($num2-$num1)<=20)
							{
								for($i=$num1;$i<=$num2;$i++)
								{
									if(length($tmp1."".$prefix."".$i)>=5)
									{
										if($i>$num1)
										{
											print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$prefix."".$i."<\/ID>\n";
										}
										else
										{
											print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$i."<\/ID>\n";
										}
									}
								}
							}
							else
							{
								print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$num1."<\/ID>\n";
								print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$num2."<\/ID>\n";
							}
						}
						elsif($str_2 =~ /^(.+[^0-9])([0-9]+)\-(.+[^0-9])([0-9]+)$/)#numeration
						{
							$tmp1=$1;
							$num1=$2;
							$tmp2=$3;
							$num2=$4;
							if(lc($tmp1) eq lc($tmp2))
							{
								my $prefix="";
								if($num1=~/^([0]+)/)
								{
									$prefix=$1;
								}
								if(($num2-$num1)<=20)
								{
									for($i=$num1;$i<=$num2;$i++)
									{
										if(length($tmp1."".$prefix."".$i)>=5)
										{
											if($i>$num1)
											{
												print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$prefix."".$i."<\/ID>\n";
											}
											else
											{
												print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$i."<\/ID>\n";
											}
										}
									}
								}
								else
								{
									print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$num1."<\/ID>\n";
									print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$tmp1."".$num2."<\/ID>\n";
								}
							}
						}
						if(length($str_2)>=5)
						{
							print output "<ID pmid='$pmid' sid='".$sid."' start='".$start_site."' end='".$end_site."'>".$str_2."<\/ID>\n";
						}
					}
				} 	
			}
		}		
	}
	close output;	
	return 1;
}

return 1;