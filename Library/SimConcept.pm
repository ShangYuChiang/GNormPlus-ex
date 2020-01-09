#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

require './Library/Text/English.pm';
	
sub SimConcept
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sa_extraction="./tmp/".$filename.".sa.xml";
	
	my @Array_mention;
	open input,"<$sa_extraction";
	open output,">$sa_extraction.data";
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		push(@Array_mention,$tmp);
		if($tmp=~/<(.+) pmid='(.+)' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.*)'>(.*)<\/.+>/)
		{
			my $mention_type=$1;
			my $pmid=$2;
			my $sid=$3;
			my $start=$4;
			my $last=$5;
			my $tax_id=$6;
			my $mention=$7;
			my $mention_org=$7;
			
			$mention=~s/([0-9])([A-Za-z])/$1 $2/g;
			$mention=~s/([A-Z])([a-z])/$1 $2/g;
			$mention=~s/([a-z])([A-Z])/$1 $2/g;
			$mention=~s/([A-Za-z])([0-9])/$1 $2/g;
			$mention=~s/([\W\-\_])/ $1 /g;
			$mention=~s/[ ]+/ /g;
			
			#tokens
			my @tokens = split(" ",$mention);
			
			#Stemming
			my @stemmed_tokens = Text::English::stem( @tokens ); 
			
			#Stemming
			my @states = split("",$state);
			my %repeat_hash=();
			my $count_token=0;
			foreach my $token(@tokens)
			{
				if($stemmed_tokens[$count_token] eq "")
				{
					$stemmed_tokens[$count_token]=$token;
				}
				my $temmed=$stemmed_tokens[$count_token];
				my $STA=$states[$count_token];
				
				#white space
				my $WhiteSpaceFront="WSF:Y";
				my $WhiteSpaceBack="WSB:Y";
				my $token_rev=$token;
				$token_rev=~s/([\W\-\_])/\\$1/g;
				if($mention_org=~/^[ ]+$token_rev(.+)$/)
				{
					$mention_org=$1;
					$WhiteSpaceFront="WSF:N";
				}
				elsif($mention_org=~/^$token_rev(.+)$/)
				{
					$mention_org=$1;
				}
				elsif($mention_org=~/^[ ]+$token_rev$/)
				{
					$WhiteSpaceFront="WSB:N";
					$WhiteSpaceBack="WSB:N";
				}
				elsif($mention_org=~/^$token_rev$/)
				{
					$WhiteSpaceFront="WSB:Y";
					$WhiteSpaceBack="WSB:N";
				}
				else
				{
					print "error: (WSF) $mention	$mention_org\n";
					exit;
				}
				
				if($mention_org=~/^[ ]+/)
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
				if($tmp=~/^(or|and)$/){$MentionType="-Type_Conjunction-";}
				elsif($tmp=~/^[\/]$/ && (length($tokens[$count_token-1]) == length($tokens[$count_token+1]))){$MentionType="-Type_Continuous-";}
				elsif($tmp=~/^(to|\-)$/ && $tokens[$count_token-1]=~/[0-9\-]+/){$MentionType="-Type_Numerous-";}
				elsif($token!~/[\W\-\_]/ && (exists $repeat_hash{$token}) && (not exists $repeat_hash{","}) && (not exists $repeat_hash{"and"}) && ($tokens[$count_token-1] eq "to" || $tokens[$count_token-1] eq "-")){$MentionType="-Type_Repeat-";}
				elsif($tmp=~/^(ytochrome|cytochrome)$/){$MentionType="-Type_cytochrome-";}
				elsif($tmp=~/(related|regulated|associated|correlated|reactive)/){$MentionType="-Type_relation-";}
				elsif($lc_tmp=~/(polymorphism|mutation|deletion|insertion|duplication|genotype|genotypes)/){$MentionType="-Type_VariationTerms-";}
				elsif($tmp=~/(adrenoceptor|kinase|kinese|receptor|subunit|unit|transporter|transcription|antigen|protein|gene|oxidase|factor|transferase|transferases|member|molecule|regulator|domain)/){$MentionType="-Type_suffix-";}
				elsif($tmp=~/(alpha|beta|gamma|kappa|theta|delta)/){$MentionType="-Type_strain-";}
				elsif($tmp=~/^[abcdefgyr]$/){$MentionType="-Type_strain-";}
				else{$MentionType="__nil__";}
				
				#Protein symbols
				my $uc_tmp=$token;
				my $lc_tmp=lc($token);
				my $ProteinSym="";
				if($lc_tmp=~/(glutamine|glutamic|leucine|valine|isoleucine|lysine|alanine|glycine|aspartate|methionine|threonine|histidine|aspartic|asparticacid|arginine|asparagine|tryptophan|proline|phenylalanine|cysteine|serine|glutamate|tyrosine|stop|frameshift)/){$ProteinSym="-ProteinSymFull-";}
				elsif($lc_tmp=~/^(cys|ile|ser|gln|met|asn|pro|lys|asp|thr|phe|ala|gly|his|leu|arg|trp|val|glu|tyr|fs|fsx)$/){$ProteinSym="-ProteinSymTri-";}
				elsif($lc_tmp=~/^(ys|le|er|ln|et|sn|ro|ys|sp|hr|he|la|ly|is|eu|rg|rp|al|lu|yr)$/ && $last_token=~/^[CISQMNPKDTFAGHLRWVEYX]$/){$ProteinSym="-ProteinSymTriSub-";}
				elsif($uc_tmp=~/^[CISQMNPKDTFAGHLRWVEYX]$/){$ProteinSym="-ProteinSymChar-";}
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

				print output $token." ".$temmed." ".$Num_num." ".$Num_Uc." ".$Num_Lc." ".$Num_All." ".$SpecificC." ".$CHEM." ".$MentionType." ".$ProteinSym." ".$Pattern1." ".$Pattern2." ".$Pattern3." ".$Pattern4." ".$prefix." ".$suffix." ".$WhiteSpaceFront." ".$WhiteSpaceBack." ".$STA."\n";
				$count_token++;
				$repeat_hash{$token}=1;
			}
			print output "                         \n";
		}
	}
	close output;
	close input;
	
	if($^O=~/Win/)
	{
		open(CRF,"|CRFmodel/crf_test -m CRFmodel/Model.SimConcept -o $sa_extraction.output $sa_extraction.data");
		close CRF;
	}
	else
	{
		open(CRF,"|./CRFmodel/crf_test -m CRFmodel/Model.SimConcept -o $sa_extraction.output $sa_extraction.data")|| die ("\nAn Error in CRF module. \nPlease re-install CRF module:\n\n\tsh Installation.sh\n\n");	
		close(CRF);
	}
	
=pod
	my $count=0;
	my %POFT_hash=();
	my %POFT_CN_hash=();
	open input,"<$sa_extraction.output";
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		#IDH	idh	N:0	U:3	L:0	A:0	__nil__	__nil__	__nil__	__nil__	P1:AAA	P2:aaa	P3:A	P4:a	I	ID	IDH	__nil__	__nil__	H	DH	IDH	__nil__	__nil__	A
		if($tmp=~/^([^\t]+).*([ASCNO])$/)
		{
			$POFT_hash{$count}=$POFT_hash{$count}.$1." ".$2."|";
			if($2 eq "C"){$POFT_CN_hash{$count}="C";}
			elsif($2 eq "N"){$POFT_CN_hash{$count}="N";}
		}
		else
		{
			$POFT_hash{$count}=~s/\|$//g;
			$count++;
		}
	}
	close input;
	open output,">$sa_extraction";
	for(my $i=0;$i<$count;$i++)
	{
		my $pre="";my $mention="";my $post="";
		
		if($Array_mention[$i]=~/(<.+ pmid='.+' sid='.+' start='.+' end='.+' tax_id='.*')>(.*)(<\/.+>)/){$pre=$1;$mention=$2;$post=$3;}
		
		if(not exists $POFT_CN_hash{$i})
		{
			print output $Array_mention[$i]."\n";
		}
		elsif($POFT_CN_hash{$i}=~/[CN]/)
		{
			my @strains;
			my $strain="";my $CNO_continous=0;
			my $A="";
			my @split_tmp=split(/\|/,$POFT_hash{$i});
			foreach my $tmp (@split_tmp)
			{
				if($tmp=~/^(.+) (.+)$/)
				{
					my $token=$1;
					my $state=$2;
					if($state=~/[A]/)
					{
						$A=$A.$token." ";
					}
					elsif($state=~/[S]/ && $count_S==0)
					{
						$A=~s/s $/ /g;
						$A=$A."STRAIN";
						$strain=$strain.$token." ";
						$CNO_continous=0;
					}
					elsif($state=~/[CN]/ && $CNO_continous==0)
					{
						push(@strains,$strain);
						$strain="";
						$S_count=0;
						$CNO_continous++;
					}
				}
			}
			if($strain ne ""){push(@strains,$strain);}
			$A=~s/(STRAIN){1,}/STRAIN/g;
			$A=~s/s $/ /g;
			
			if($POFT_CN_hash{$i} eq "C")
			{
				foreach $strain (@strains)
				{
					my $tmp=$A;
					$tmp=~s/STRAIN/$strain/g;
					$tmp=~s/  / /g;$tmp=~s/ $//g;
					print output $pre." org='".$mention."'>".$tmp.$post."\n";
				}
			}
			elsif($POFT_CN_hash{$i} eq "N")
			{
				my $strain1=$strains[0];
				my $strain2=$strains[1];
				if(($strain2-$strain1)<=20)
				{
					for(my $i=$strain1;$i<=$strain2;$i++)
					{
						my $tmp=$A;
						$tmp=~s/STRAIN/$i/g;
						$tmp=~s/  / /g;$tmp=~s/ $//g;
						print output $pre." org='".$mention."'>".$tmp.$post."\n";
					}
				}
				else
				{
					my $tmp=$A;
					$tmp=~s/STRAIN/$strain1/g;
					$tmp=~s/  / /g;$tmp=~s/ $//g;
					print output $pre." org='".$mention."'>".$tmp.$post."\n";
				}
			}
		}
	}
	close output;
=cut
	my %Abb_hash=();
	open Abb,"<tmp/".$filename.".Ab3P";
	while(<Abb>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$tmp=<Abb>;
		while(<Abb>)
		{
			$tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/^  (.+)\|(.+)\|(.+)$/)
			{
				my $SF=$1;my $LF=$2;
				$SF=~s/([0-9])([A-Za-z])/$1 $2/g;
				$SF=~s/([A-Z])([a-z])/$1 $2/g;
				$SF=~s/([a-z])([A-Z])/$1 $2/g;
				$SF=~s/([A-Za-z])([0-9])/$1 $2/g;
				$SF=~s/([\W\-\_])/ $1 /g;
				$SF =~ s/[ ]+/ /g;
				$LF=~s/([0-9])([A-Za-z])/$1 $2/g;
				$LF=~s/([A-Z])([a-z])/$1 $2/g;
				$LF=~s/([a-z])([A-Z])/$1 $2/g;
				$LF=~s/([A-Za-z])([0-9])/$1 $2/g;
				$LF=~s/([\W\-\_])/ $1 /g;
				$LF =~ s/[ ]+/ /g;
				$Abb_hash{$SF}=$LF;
				$Abb_hash{$LF}=$SF;
			}
			else
			{
				last;
			}
		}
	}
	close Abb;
	
	my $MentionCount=0;
	my %Tokens_hash=();
	my %States_hash=();
	open input,"<$sa_extraction.output";
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^([^\t]+).*([ASCJBEFNOP])$/)
		{
			my $tkn=$1;
			my $sta=$2;
			#C:Continuous
			#J:conJunction
			#B:abbreviation-Begin
			#E:abbreviation-End
			#F:abbreviation-Follow with suffix/strain
			#N:Numerous
			#O:others	
			#P:its
			$Tokens_hash{$MentionCount}=$Tokens_hash{$MentionCount}.$tkn." ";
			$States_hash{$MentionCount}=$States_hash{$MentionCount}.$sta;
		}
		else
		{
			if($Tokens_hash{$MentionCount}=~/[\(\[]/)
			{
				foreach my $Abb (keys %Abb_hash)
				{
					my $Abb1_tmp=$Abb;
					my $Abb2_tmp=$Abb_hash{$Abb};
					$Abb1_tmp=~s/[\(\)\[\]]/\./g;
					$Abb2_tmp=~s/[\(\)\[\]]/\./g;
					if($Tokens_hash{$MentionCount}=~/^(.*)$Abb1_tmp [\(\[] $Abb2_tmp [\)\]] (.+)$/)
					{
						my @Num1=split(" ",$1);
						my @NumAbbL=split(" ",$Abb);
						my @NumAbbS=split(" ",$Abb_hash{$Abb});
						my @Num2=split(" ",$2);
						#print @Num1."\t".@NumAbbL."\t".@NumAbbS."\t".@Num2."\n";
						my $N1=substr($States_hash{$MentionCount},0,@Num1);
						my $N2=substr($States_hash{$MentionCount},@Num1,@NumAbbL);
						my $N3=substr($States_hash{$MentionCount},@Num1+@NumAbbL+1,@NumAbbS);
						my $N4=substr($States_hash{$MentionCount},@Num1+@NumAbbL+@NumAbbS+2,@Num2);
						$N2=~s/./W/g;$N3=~s/./W/g;
						#print $States_hash{$MentionCount}."\n";
						$States_hash{$MentionCount}=$N1.$N2."X".$N3."O".$N4;
						#print $States_hash{$MentionCount}."\n";
					}
				}
			}
			$MentionCount++;
		}
	}
	close input;
	
	my $TargetTitle="";
	my $TargetAbstract="";
	my %TargetMention_hash=();
	my $TargetCount=0;
	open sa,"<$sa_extraction";
	while(<sa>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/<(.+) pmid='(.+)' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.*)'>(.*)<\/.+>/)
		{
			my $mention_type=$1;
			my $pmid=$2;
			my $sid=$3;
			my $start=$4;
			my $last=$5;
			my $tax_id=$6;
			my $mention=$7;
			my $mention_org=$7;
			
			if($mention=~/^.+\((.+)\)$/ && length($1)<=2)
			{
				$TargetMention_hash{$TargetCount}=$pmid."\t".$sid."\t".$start."\t".$last."\t".$mention."\t".$tax_id."\t".$mention_type;
				if($States_hash{$TargetCount}=~/(.*)(B[AS]+E)$/)
				{
					my $pre=$1;
					my $BE=$2;$BE=~s/./S/g;
					$States_hash{$TargetCount}=$pre.$BE;
				}
				elsif($States_hash{$TargetCount}=~/(.*)(B[AS]+)$/)
				{
					my $pre=$1;
					my $BE=$2;$BE=~s/./S/g;
					$States_hash{$TargetCount}=$pre.$BE;
				}
				elsif($States_hash{$TargetCount}=~/(.*)([AS]+E)$/)
				{
					my $pre=$1;
					my $BE=$2;$BE=~s/./S/g;
					$States_hash{$TargetCount}=$pre.$BE;
				}
			}
			else
			{
				$TargetMention_hash{$TargetCount}=$pmid."\t".$sid."\t".$start."\t".$last."\t".$mention."\t".$tax_id."\t".$mention_type;
			}
			$TargetCount++;
		}
	}
	close sa;
	
	open output,">".$sa_extraction;
	my @Output_Split_mention;
	for(my $i=0;$i<$MentionCount;$i++) # i mention
	{
		my @Split_mention=();
		my @Split_state=();
		my @mention=($TargetMention_hash{$i}=~/^[^\t]+	[^\t]+	[^\t]+	[^\t]+	([^\t]+)/);
		$mention[0]=~s/([0-9])([A-Za-z])/$1 $2/g;
		$mention[0]=~s/([A-Z])([a-z])/$1 $2/g;
		$mention[0]=~s/([a-z])([A-Z])/$1 $2/g;
		$mention[0]=~s/([A-Za-z])([0-9])/$1 $2/g;
		$mention[0]=~s/([\W\-\_])/ $1 /g;
		$mention[0] =~ s/[ ]+/ /g;
		my @tokens = split(" ",$mention[0]);
		my @states = split("",$States_hash{$i});
		my $tmp_mention="";my $tmp_state="";
		
		#Split BEJ
		for(my $j=0;$j<@states;$j++)
		{
			if($states[$j]=~/^[BEJ]$/)
			{
				push(@Split_mention,$tmp_mention);
				push(@Split_state,$tmp_state);
				$tmp_mention="";
				$tmp_state="";
			}
			else #ASCNO
			{
				$tmp_mention=$tmp_mention.$tokens[$j]." ";
				$tmp_state=$tmp_state.$states[$j];
			}
		}
		if($tmp_mention ne ""){push(@Split_mention,$tmp_mention);}
		if($tmp_state ne ""){push(@Split_state,$tmp_state);}
		
		#Split G/F
		for(my $j=0;$j<@Split_mention;$j++) # j split mention
		{
			my $A="";my $STAA="";
			my $strainX="";
			my $STAstrainX="";
			my $X_continous=0;
			my @strainsX;
			my @STAstrainsX;
			my $X="";
			my @each_token = split(" ",$Split_mention[$j]);
			my @each_state = split("",$Split_state[$j]);
			for(my $k=0;$k<@each_state;$k++) # k token
			{
				if($each_state[$k]=~/[ACNS]/) #BEJ have been split
				{
					$A=$A.$each_token[$k]." ";
					$STAA=$STAA.$each_state[$k];
				}
				elsif($each_state[$k]=~/[W]/)
				{
					$A=$A."STRAIN";
					$STAA=$STAA."@";
					$strainX=$strainX.$each_token[$k]." ";
					$STAstrainX=$STAstrainX.$each_state[$k];
					$X_continous=0;
				}
				elsif($each_state[$k]=~/[X]/ && $X_continous==0)
				{
					$X=$each_state[$k];
					push(@strainsX,$strainX);
					push(@STAstrainsX,$STAstrainX);
					$strainX="";
					$STAstrainX="";
					$X_continous++;
				}
			}
			if($strainX ne ""){push(@strainsX,$strainX);}
			if($STAstrainX ne ""){push(@STAstrainsX,$STAstrainX);}
			$A=~s/(STRAIN){1,}/STRAIN/g;
			$STAA=~s/(@){1,}/@/g;

			if($X eq "X")
			{
				for(my $l=0;$l<@strainsX;$l++) # k token
				{
					my $strain=$strainsX[$l];
					my $state=$STAstrainsX[$l];
					my $tkns=$A;$tkns=~s/STRAIN/$strain/g;$tkns=~s/  / /g;$tkns=~s/ $//g;
					my $STAAs=$STAA;$STAAs=~s/@/$state/g;
					push(@Split_mention,$tkns);
					push(@Split_state,$STAAs);
				}
				delete $Split_mention[$j];
				delete $Split_state[$j];
			}
		}
		
		foreach my $state(@Split_state)
		{
			$state=~s/W/A/g;
		}
		
		#Split CN
		for(my $j=0;$j<@Split_mention;$j++) # j split mention
		{
			my $A="";
			my $strainCN="";
			my $CNO_continous=0;
			my @strainsCN;
			my $CorN="";
			
			#refinement
			if($Split_state[$j]=~/^([A]+)(S[NC][O]*)SS$/){$Split_state[$j]=$1.$2."OS";}
			elsif($Split_state[$j]=~/^([A]+)SS(S[NC][O]*S)$/){$Split_state[$j]=$1."AA".$2;}
			elsif($Split_state[$j]=~/^([A]+)S(S[NC][O]*S)$/){$Split_state[$j]=$1."A".$2;}
			
			my @each_token = split(" ",$Split_mention[$j]);
			my @each_state = split("",$Split_state[$j]);
			for(my $k=0;$k<@each_state;$k++) # k token
			{
				if($each_state[$k]=~/[A]/)
				{
					$A=$A.$each_token[$k]." ";
				}
				elsif($each_state[$k]=~/[S]/)
				{
					if(length($A)>=4)
					{
						$A=~s/s $/ /g;
					}
					$A=$A."STRAIN";
					$strainCN=$strainCN.$each_token[$k]." ";
					$CNO_continous=0;
				}
				elsif($each_state[$k]=~/[CN]/ && $CNO_continous==0)
				{
					$CorN=$each_state[$k];
					push(@strainsCN,$strainCN);
					$strainCN="";
					$CNO_continous++;
				}
			}
			if($strainCN ne ""){push(@strainsCN,$strainCN);}
		
			$A=~s/(STRAIN){1,}/STRAIN/g;
			if($A=~/^(.+)s / && length($1)>=3)
			{
				$A=~s/s $/ /g;
			}
			if($CorN eq "C")
			{
				foreach $strain (@strainsCN)
				{
					my $tmp=$A;
					$tmp=~s/STRAIN/$strain/g;
					$tmp=~s/  / /g;$tmp=~s/ $//g;
					$Output_Split_mention[$i]=$Output_Split_mention[$i].$tmp."|";
				}
			}
			elsif($CorN eq "N")
			{
				my $strain1=$strainsCN[0];
				my $strain2=$strainsCN[1];
				if($strain1=~/[0-9]+/)
				{
					if(($strain2-$strain1)<=20)
					{
						for(my $s=$strain1;$s<=$strain2;$s++)
						{
							my $tmp=$A;
							$tmp=~s/STRAIN/$s/g;
							$tmp=~s/  / /g;$tmp=~s/ $//g;
							$Output_Split_mention[$i]=$Output_Split_mention[$i].$tmp."|";
						}
					}
					else
					{
						$Output_Split_mention[$i]=$Output_Split_mention[$i].$Split_mention[$j]."|";
					}
				}
				elsif($strain1=~/[A-Z]/)
				{
					if((ord($strain2)-ord($strain1))<=20)
					{
						for(my $s=ord($strain1);$s<=ord($strain2);$s++)
						{
							my $tmp=$A;
							my $chr_s=chr($s);
							$tmp=~s/STRAIN/$chr_s/g;
							$tmp=~s/  / /g;$tmp=~s/ $//g;
							$Output_Split_mention[$i]=$Output_Split_mention[$i].$tmp."|";
						}
					}
					else
					{
						$Output_Split_mention[$i]=$Output_Split_mention[$i].$Split_mention[$j]."|";
					}
				}
				else
				{
					$Output_Split_mention[$i]=$Output_Split_mention[$i].$Split_mention[$j]."|";
				}
			}
			else
			{
				$Output_Split_mention[$i]=$Output_Split_mention[$i].$Split_mention[$j]."|";
			}
		}
		
		$Output_Split_mention[$i]=~s/[ ]*\|[ ]*/\|/g;
		$Output_Split_mention[$i]=~s/[ ]*\-[ ]*/\-/g;
		$Output_Split_mention[$i]=~s/\|$//g;$Output_Split_mention[$i]=~s/^\|//g;$Output_Split_mention[$i]=~s/[\|]+/\|/g;
		if(exists $TargetMention_hash{$i})
		{
			$States_hash{$i}=~s/WO/AF/g;
			$States_hash{$i}=~s/W/A/g;
			$States_hash{$i}=~s/X/B/g;
			#<DomainMotif pmid='16380819' sid='Figure_26' start='53' end='77' tax_id='10245'>MV-GFP (MOI 3) or EV-GFP</DomainMotif>
			#16380819	Figure_26	53	77	MV-GFP (MOI 3) or EV-GFP	DomainMotif|MV-GFP|MOI 3|EV-GFP|AAABASEJAAA
			my ($pmid,$sid,$start,$last,$mention,$tax_id,$type)=($TargetMention_hash{$i}=~/^([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)/);
			my @mentions=split(/\|/,$Output_Split_mention[$i]);
			if(@mentions>1 && $type eq "Gene")
			{
				foreach my $splited_mention(@mentions)
				{
					print output "<$type pmid='$pmid' sid='$sid' start='$start' end='$last' tax_id='$tax_id' org='$mention'>$splited_mention</$type>\n";
				}
			}
			else
			{
				print output "<$type pmid='$pmid' sid='$sid' start='$start' end='$last' tax_id='$tax_id'>$mention</$type>\n";
			}
		}
	}
	close output;
	
	#open clear,"|rm tmp/*.sa.xml.data";	close clear;
	#open clear,"|rm tmp/*.sa.xml.output";	close clear;
	
	return 1;
}
return 1;