#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

my %Roman2Arabian_hash=();
my %TargetTax_hash=();
my %IDF_Entity_hash=();
my %IDF_Token_hash=();

my @Pmidlist;
my @Token_array;
my @Entity_array;
my %Entity_hash=();
my %Family_hash=(); #FamilyNameTokenization, PartialMatch

my %CandidateEntity_token_hash=();
my %CandidateEntity_weight_hash=();
my %CandidateEntity_tax_hash=();
my %Final_CandidateEntity_token_hash=();
my %Final_CandidateEntity_weight_hash=();
my %Final_CandidateEntity_tax_hash=();

my $finalSSThreshold;
my $simiarlarscorethreshold;
my %commonword_hash=();

my $specialScore;
my $offset;
my $normal;
my $base;

sub log10 
{
	my $n = shift;
	return log($n)/log(10);
}

sub round
{
  my($value,$rank) = @_;
  if($value>0){
    return int($value * 10**$rank + 0.5) / 10**$rank;
  }else{
    return int($value * 10**$rank - 0.4 )/ 10**$rank;
  }
}

sub Tokenization
{
	my ($pmid) = @_[0];
	my ($sa_extraction) = @_[1];
	
	my %two_Entity_hash=();
	%Entity_hash=();
	
	open input,"<$sa_extraction";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/<Gene pmid='$pmid' sid='(.+)' start='.+' end='.+' tax_id='(.+)' org='(.+?)' org='.+'>(.+)<\/Gene>/ ||
			$tmp=~/<Gene pmid='$pmid' sid='(.+)' start='.+' end='.+' tax_id='(.+)' org='(.+)'>(.+)<\/Gene>/ ||
			$tmp=~/<Gene pmid='$pmid' sid='(.+)' start='.+' end='.+' tax_id='(.+)'>(.+)<\/Gene>/
		)
		{
			my @Entitys;
			my $sid=$1;
			my $Tax_id=$2;
			my $mention=$3;
			my $mention2=$4;
			push(@Entitys,$mention);
			if($mention2 ne "")
			{
				push(@Entitys,$mention2);
			}
			
			foreach my $Entity_org(@Entitys)
			{
				if( length($Entity_org)==2)
				{
					$two_Entity_hash{$Entity_org} = $Tax_id;
				}
				elsif( length($Entity_org)>=3)
				{
				
					$Entity = $Entity_org;
					$Entity =~ s/([A-Za-z]+)([0-9]+)/$1 $2/g;
					$Entity =~ s/([0-9]+)([A-Za-z]+)/$1 $2/g;
					$Entity =~ s/[ ]+/ /g;
					$Entity =~ s/^[ ]+//g;
					
					if($Entity =~/^(.+) (binding|binding protein[s]*|homeoprotein[s]*|interacting|interacting protein[s]*|protein[s]*|gene[s]*|homologue|homolog|promoter)$/) #-protein(s)
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}	
					
					if($Entity =~/^(.+) molecule$/) #molecule[s]
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}	
					elsif($Entity =~/^(\S+) (and) (\S+) molecules$/)
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}	
									
					if($Entity =~/^(.+) subunit (.*)$/) #-subunit(s)
					{
						$temp1=$1;$temp2=$2;
						$tmp_symbol = $temp1."".$temp2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $temp1." ".$temp2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $temp1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}	
					elsif($Entity =~/^(.+) subunits/)
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}

					if($Entity =~/.+(protein|gene|ase)[s]* ([A-Za-z0-9 ]+)$/) #-(protein|gene|ase)[s]*
					{
						$tmp_symbol = $2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if( $tmp_symbol=~/[0-9]/)
						{
							if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						}
					}
					
					if($Entity =~/(.+) associated(.+)/)
					{
						$tmp_symbol = $1."".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						$tmp_symbol =~ s/[ ]+/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}	
					
					if($Entity =~/^(.+)[\W\-\_]*([IX]{2,})$/)
					{
						my $tmp_symbol=$1." ".$Roman2Arabian_hash{$2};$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					if($Entity =~/^(.+)[\W\-\_]*([IX]{2,})([a-cA-C])$/)
					{
						my $tmp_symbol=$1." ".$Roman2Arabian_hash{$2}." ".$3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)alpha$/)
					{
						$temp1=$1;
						if(length($temp1) >= 3)
						{
							my $tmp_symbol=$temp1."a";$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
							if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						}
					}
					elsif($Entity =~/^(.+)beta$/)	
					{
						$temp1=$1;
						if(length($temp1) >= 3)
						{
							my $tmp_symbol=$temp1."b";$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
							if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						}
					}
					elsif($Entity =~/^(.+)alpha(.+)$/)	
					{
						$tmp_symbol=$1."alpha ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)beta(.+)$/)	
					{
						$tmp_symbol=$1."beta ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)gamma(.+)$/)	
					{
						$tmp_symbol=$1."gamma ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
							
					if($Entity =~/(KIAA [0-9]+)/)
					{
						$tmp_symbol=$1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}

					if($Entity =~/^(.+[0-9A-Z]) [p]$/)
					{
						$temp1=$1;
						$tmp_symbol=$temp1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol=$temp1." protein";$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+[A-Z])[A]$/)
					{
						$tmp_symbol=$1." A";$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					
					if($Entity =~/^(.+)nu$/)
					{
						$tmp_symbol=$$1." nu";$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					
					if($Entity =~/^(.+)s ([A-Za-b1-9]), ([A-Za-b1-9]), ([A-Za-b1-9]){0,1}, (and) ([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1." ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$4;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$6;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)s ([A-Za-b1-9]), ([A-Za-b1-9]) (and) ([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1." ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$5;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)s ([A-Za-b1-9]) (and) ([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1." ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$4;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+) ([A-Za-b1-9]) (and) ([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1." ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$4;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+) (and) (.+)$/)
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)s ([A-Za-b1-9])[,\/]([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1." ".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1." ".$3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)([A-Za-b1-9])[,\/]([A-Za-b1-9])$/)
					{
						$tmp_symbol = $1."".$2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $1."".$3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					
					if($Entity =~/^(.+) ([0-9]) (to) (.+) ([0-9])$/)
					{
						$tmp1 = $1;
						$num1 = $2;
						$tmp2 = $4;
						$num2 = $5;			
						if($tmp1 eq $tmp2)
						{
							for($i=$num1;$i<=$num2;$i++)
							{
								$tmp_symbol =$tmp1." ".$i;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
								if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
							}
						}
					}
					elsif($Entity =~/^(.+) (or) (.+)$/)
					{
						$tmp_symbol = $1;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol = $3;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					
					if($Entity =~/^(.+)[(](.+)[)](.+)$/)
					{
						$tmp_symbol1=$1;
						$tmp_symbol2=$2;
						$tmp_symbol3=$3;
						$tmp_symbol=$tmp_symbol1." ".$tmp_symbol3;
						$tmp_symbol =~ s/[ ]+/ /g;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						
						if($tmp_symbol3 eq " protein" || $tmp_symbol3 eq " gene")
						{
							$tmp_symbol=$tmp_symbol2;
							$tmp_symbol =~ s/[ ]+/ /g;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
							if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						}
					}
					if($Entity =~/^(.+) [(]([abABIX]+)[)]$/)
					{
						$tmp1 = $1;
						$tmp2 = $2;
						$tmp_symbol=$tmp1." ".$tmp2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						$tmp_symbol=$tmp1."".$tmp2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
					}
					elsif($Entity =~/^(.+)[(](.+?)[)]$/)
					{
						$tmp_symbol = $1;
						$tmp_symbol_2 = $2;
						if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						if(length($tmp_symbol_2)>2)
						{
							$tmp_symbol=$tmp_symbol_2;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
							if(exists $Entity_hash{$tmp_symbol}){if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/){}else{$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";}}else{$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";}
						}
					}
					$Entity =~ s/[()]/ /g;
					$Entity =~ s/[ ]+/ /g;
					
					$tmp_symbol=$Entity;$tmp_symbol =~ s/[-\+\(\)\[\]\'!@\#\$\&\%:;.\=*><\"]/ /g;
					if(exists $Entity_hash{$tmp_symbol})
					{
						if($Entity_hash{$tmp_symbol}=~/\@$Tax_id\@/)
						{}
						else
						{
							$Entity_hash{$tmp_symbol} = $Entity_hash{$tmp_symbol}.$Tax_id."@";
						}
					}
					else
					{
						$Entity_hash{$tmp_symbol} = "@".$Tax_id."@";
					}
				}
			}
		}
	}
	close input;
	
	@Token_array=@empty_array;
	@Entity_array=@empty_array;	
	
	foreach  my $ent (keys %Entity_hash)
	{
		$Entity_hash{$ent}=~s/^\@(.+)\@$/$1/g;
		if( length($ent)>=3)
		{
				my @tokens = split(" ",$ent);
				push @Token_array, @tokens;
				push @Entity_array, $ent;				
		}
	}
	return 1;
}

sub Tokenization_FamilyDomain
{
	my ($pmid) = @_[0];
	my ($sa_extraction) = @_[1];
	
	#my %Tmp_Entity_hash=();
	foreach my $ent (keys %Entity_hash){$Entity_hash{$ent}="@".$Entity_hash{$ent}."@";}
	
	open input,"<$sa_extraction";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^<(FamilyName|DomainMotif) pmid='$pmid' sid='[^ ]+' start='[^ ]+' end='[^ ]+' tax_id='([^ ]+)'>(.+)<\/(FamilyName|DomainMotif)>/)
		{
			my $Tax_id = $2;
			my $Entity = $3;
			$Entity =~ s/([A-Za-z]+)([0-9]+)/$1 $2/g;
			$Entity =~ s/([0-9]+)([A-Za-z]+)/$1 $2/g;
			$Entity=~ s/[\W\-\_]/ /g;
			$Entity =~ s/[ ]+/ /g;
			$Entity =~ s/^[ ]+//g;
			$Entity =~ s/[ ]+$//g;
			
			if(exists $Entity_hash{$Entity})
			{
				if($Entity_hash{$Entity}=~/\@$Tax_id\@/){}
				else{
					$Entity_hash{$Entity} = $Entity_hash{$Entity}.$Tax_id."@";
				}
			}
			else
			{
				$Entity_hash{$Entity} = "@".$Tax_id."@";
			}
		}
	}
	close input;
	
	foreach my $Entity (keys %Entity_hash)
	{
		$Entity_hash{$Entity}=~s/^\@(.+)\@$/$1/g;
		if( length($Entity)>=3)
		{
			my @tokens = split(" ",$Entity);
			push (@Token_array, @tokens);
			push (@Entity_array, $Entity);
			$Entity=~s/[\W\-\_]//g;
			$Family_hash{lc($Entity)}=1;
		}
	}
	return 1;
}	

sub similarity
{
	my ($sim_entrez_id) = @_[0];
	my ($exact_or_partial) = @_[1];
	my ($suitable_tax) = @_[2];
	my @temp_Token_array;
	my %tem_pattern_hash=();
	my %tmp_IDF_hash=();
	my $score=0;
	my $score2=0;
	
	my $dic_term="";
	if( $exact_or_partial == 0 )
	{
		$dic_term="Token";
		%tmp_IDF_hash=%IDF_Token_hash;
	}
	else
	{
		$dic_term="Entity";
		%tmp_IDF_hash=%IDF_Entity_hash;
	}
	
	open output,"<".$dictionary."/TF_".$dic_term."/Tax_".$suitable_tax.".txt";
	while(<output>)
	{
		my $tmp=$_;
		if($tmp=~/^(.+)	(.+)	(.+)	(.+)	(.+)$/ && ($1 == $sim_entrez_id))
		{
			$pattern=$2;
			$Cj=$3;
			$Norm=$3;
			if($Norm==0){$Norm=1;}
			$maximum=$5;
		}
	}
	close output;
	
	@temp_Token_array = @Token_array;
	@split_pattern = split(",",$pattern);
	foreach $each_pattern(@split_pattern)
	{
		($token,$freq) = split("-",$each_pattern);
		$token =~s/ $//;
		$tem_pattern_hash{"\L$token"}=$freq;
	}
	$IDFi=1;
	foreach $entity(@temp_Token_array)
	{
		$entity="\L$entity";
		if(exists $tem_pattern_hash{$entity})
		{
			$IDFi=1;
			if(exists $commonword_hash{$entity})
			{
				$IDFi=1;
			}
			elsif(exists $tmp_IDF_hash{$entity})
			{
				$IDFi=$tmp_IDF_hash{$entity};
			}
		}
		if($tem_pattern_hash{$entity}==$maximum )
		{
			$score = $score + $IDFi;
		}
		else
		{
			$Fij = $tem_pattern_hash{"\L$entity"} / ($maximum+1);
			$score=$score + $Fij * $IDFi * (1/(1-$Fij));
		}
	}
	if($Norm==0){$Norm=1;}
	$score = $Cj * (1/$Norm) * $score;
	return $score;
}

sub OfficialName 
{
	my ($dictionary) = @_[0]; 
	my ($entity) = @_[1];
	my ($TAX) =	@_[2];
	
	my $entity_nospace = $entity;
	$entity_nospace =~ s/[\W\-\_]+//g;	
	my %entrezid_hash=();
	open input,"<".$dictionary."/OfficialName/OfficialName_".$TAX.".txt";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^(.+)::(.+)$/)
		{
			my $Extra_entity=$1;
			my $Extra_entrezids=$2;
			if(lc($entity_nospace) eq lc($Extra_entity))
			{
				my @split_Extra_entrezids=split(/\|/,$Extra_entrezids);
				foreach my $Extra_entrezid(@split_Extra_entrezids)
				{
					$Extra_entrezid=~s/[\n\r\t\0]//g;
					$entrezid_hash{$Extra_entrezid}=1;
				}
			}
		}
	}
	close input;

	my $return_ids="";
	foreach my $entrezid(keys %entrezid_hash)
	{
		$return_ids=$return_ids.$entrezid."|";
		if(exists $CandidateEntity_token_hash{$entrezid})
		{		
			$CandidateEntity_token_hash{$entrezid}=$CandidateEntity_token_hash{$entrezid}."\|\L$entity*";
			$CandidateEntity_weight_hash{$entrezid}= $CandidateEntity_weight_hash{$entrezid} + ($specialScore*1.5) + similarity($entrezid,1,$TAX);
			$if_null++;
		}
		else
		{
			$CandidateEntity_token_hash{$entrezid} = "\L$entity*";
			$CandidateEntity_weight_hash{$entrezid} = ($specialScore*1.5) + similarity($entrezid,1,$TAX);
			$CandidateEntity_tax_hash{$entrezid} = $TAX;
			$if_null++;
		}	
	}
	$return_ids=~s/\|$//g;
	return $return_ids;
}

sub ModifiedName 
{
	my ($dictionary) = @_[0]; 
	my ($entity) = @_[1];
	my ($TAX) = @_[2];
	
	my $entity_nospace = $entity;
	$entity_nospace =~ s/[\W\-\_]+//g;
	my %entrezid_hash=();
	open input,"<".$dictionary."/ModifiedName/ModifiedName_".$TAX.".txt";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^(.+)::(.+)$/)
		{
			my $Extra_entity=$1;
			my $Extra_entrezids=$2;
			if(lc($entity_nospace) eq lc($Extra_entity))
			{
				my @split_Extra_entrezids=split(/\|/,$Extra_entrezids);
				foreach my $Extra_entrezid(@split_Extra_entrezids)
				{
					$Extra_entrezid=~s/[\n\r\t\0]//g;
					$entrezid_hash{$Extra_entrezid}=1;
				}
			}
		}
	}
	close input;

	foreach my $entrezid(keys %entrezid_hash)
	{
		if(exists $CandidateEntity_token_hash{$entrezid})
		{		
			$CandidateEntity_token_hash{$entrezid}=$CandidateEntity_token_hash{$entrezid}."\|\L$entity+";
			$CandidateEntity_weight_hash{$entrezid}= $CandidateEntity_weight_hash{$entrezid} + $specialScore + similarity($entrezid,1,$TAX);
			$if_null++;
		}
		else
		{
			$CandidateEntity_token_hash{$entrezid} = "\L$entity+";
			$CandidateEntity_weight_hash{$entrezid} = $specialScore + similarity($entrezid,1,$TAX);
			$CandidateEntity_tax_hash{$entrezid} = $TAX;
			$if_null++;
		}
	}
	return 1;
}

sub PartialMatch
{
	my ($dictionary) = @_[0]; 
	my ($entity) = @_[1];
	my ($TAX) =	@_[2];
	
	my %partial_match_weight_hash=();
	my %partial_match_entrez_hash=();
	my @empty_array;
	my @split_entities_Roman2Arabian;
	my %tokens_hash=();			
	#=====
	#split tokens
	my $entity_org=$entity;
	$entity_org=~s/([0-9])([A-Za-z])/$1 $2/g;
	$entity_org=~s/([A-Za-z])([0-9])/$1 $2/g;
	my @tokens = split(" ",$entity);
	push(@tokens,split(" ",$entity_org));
	foreach my $token(@tokens)
	{
		$tokens_hash{lc($token)}=1;
	}
	
	#=====
	#Partial Match: Find the partial match candidate EntrezID
	open input,"<".$dictionary."/PartialName/PartialName_".$TAX.".txt";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/^(.+)\:\:(.+)$/)
		{
			my $tok=$1;
			if($tok!~/^[0-9iv]+$/ && exists $tokens_hash{$tok})
			{
				my $entrezids=$2;
				my @split_entrezid=split(/\|/,$entrezids);
				foreach my $partial_entrez_id (@split_entrezid)
				{
					$partial_entrez_id=~s/[\n\r\t\0]//g;
					$partial_match_entrez_hash{$partial_entrez_id}=1;		
				}
			}
		}
	}
	close input;

		
	foreach my $each_entrez_id (keys %partial_match_entrez_hash )
	{
		if(exists $CandidateEntity_token_hash{$each_entrez_id})
		{
			$entity=lc($entity);
			my $entity_tmp=$entity;$entity_tmp=~s/[\W\-\_]/\./g;
			my $CandidateEntity_token_tmp=$CandidateEntity_token_hash{$each_entrez_id};
			$CandidateEntity_token_tmp=~s/[\+\*]/ /g;
			my @split_entity=split(/\|/,$CandidateEntity_token_tmp);
			my $match=0;
			foreach my $ent(@split_entity)
			{
				if($ent eq $entity_tmp)
				{
					$match=1;
					last;
				}
			}
			if($match==0)
			{
				$entity_tmp=~s/[\W\-\_]//g;
				if(exists $Family_hash{$entity_tmp})
				{
					$CandidateEntity_token_hash{$each_entrez_id}=~s/\+/\*/g;
				}
				$CandidateEntity_token_hash{$each_entrez_id} = $CandidateEntity_token_hash{$each_entrez_id}."\|$entity";
			}
			$CandidateEntity_weight_hash{$each_entrez_id} = $CandidateEntity_weight_hash{$each_entrez_id} + similarity($each_entrez_id,0,$TAX);
		}
	}
	return 1;
}	

sub ExactMatch
{
	my ($dictionary) = @_[0];
	my ($entity) = @_[1];
	my ($TAX) = @_[2]; 
	
	#Match
	my $entity_tmp=$entity;
	$entity_tmp=~s/ //g;
	
	if(length($entity_tmp)>=3)
	{
		my $FindIDs=OfficialName($dictionary,$entity,$TAX);
		if($FindIDs eq "")
		{
			ModifiedName($dictionary,$entity,$TAX);
		}
	}
	return 1;
}

sub ID_Match
{
	my ($dictionary) = @_[0];
	my ($pmid) = @_[1];
	my ($table_extraction) = @_[2];
	my ($id_extraction) = @_[3];
	my ($gn_extraction) = @_[4];
	
	my %CID_hash=();
	my %CIDindex_hash=();
 	my %dictionary_hash=();
 	
 	open input,"<$id_extraction";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/<ID pmid='$pmid' .*>(.+)<\/ID>/)
		{
			my $symbol=lc($1);
			$symbol =~s/\W//g;
			$symbol =~s/\_//g;
			$symbol =~s/ //g;
			if(length($symbol)>=6 && length($symbol)<=30)
			{
				$CID_hash{$symbol}=1;
				if(lc(substr($symbol,-1,1))=~/[a-z]/)
				{
					$CIDindex_hash{"letter"}=1;
				}
				else
				{
					$CIDindex_hash{substr($symbol,-1,1)}=1;
				} 
			}
		}
	}
	close input;
	
	open input,"<$table_extraction";
	while(<input>)
	{
		my $tmp=$_;
		if($tmp=~/<TABLE pmid='$pmid'>(.+)<\/TABLE>/)
		{
			my $symbol=lc($1);
			$symbol =~s/[\W\-\_]//g;
			if(length($symbol)>=6 && length($symbol)<=30)
			{
				$CID_hash{$symbol}=1;
				if(lc(substr($symbol,-1,1))=~/[a-z]/)
				{
					$CIDindex_hash{"letter"}=1;
				}
				else
				{
					$CIDindex_hash{substr($symbol,-1,1)}=1;
				} 
			}
		}
	}
	close input;
	
	foreach my $symbol (keys %CIDindex_hash)
	{
		open input,"<".$dictionary."/CID/CID_".$symbol.".txt";
		while(<input>)
		{
			$tmp=$_;
			if($tmp=~/^(.+)	(.+)	(.+)$/)
			{
				my $tax_id=$1;
				my $entrez_id=$2;
				my $Cid=lc($3);
				$Cid=~s/[\n\r\0\W\-\_]//g;
				$dictionary_hash{$Cid}=$tax_id."|".$entrez_id;
			}
		}
		close input;
	}
	
	open output,">>$gn_extraction";
	foreach my $symbol (keys %CID_hash)
	{
		if (exists $dictionary_hash{lc($symbol)})
		{
			$dictionary_hash{lc($symbol)}=~/^(.+)\|(.+)$/;
			my $tax_id=$1;
			my $entrez_id=$2;
			my $weight=(($specialScore*1.5)-$offset)/$normal+$base;
			printf output "<GN pmid='$pmid' tax_id='".$tax_id."' entrezID='".$entrez_id."' score='%5.4f'>".$symbol."*<\/GN>\n",$weight;
		}			
	}
	close output;
	return 1;
}

sub GNPostProcessing
{
	my ($filename) = @_[0];
	my ($TargetPmid) = @_[1];
	my ($sa_extraction) = @_[2];
	my ($gn_extraction) = @_[3];
	
	my %ShortForm2LongForm_hash=();
	my %ShortForm2Score_hash=();
	
	my $tmpGN="";
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
		}
		elsif($tmp=~/^  (.+)\|(.+)\|(.+)$/)#$abbr_hash{$pmid}=1;
		{
			
			if($pmid eq $TargetPmid)
			{
				my $ShortForm=lc($1);
				my $LongForm=lc($2);
				my $Score=$3;
				$ShortForm=~s/[\W\-\_]//g;
				$ShortForm2LongForm_hash{$ShortForm}=$LongForm;
				$ShortForm2Score_hash{$ShortForm}=$Score;
			}
		}
		else
		{
			if($pmid eq $TargetPmid)
			{
				my %FamilyDomain_hash=();
				open sa,"<$sa_extraction";
				while(<sa>)
				{
					my $tmp=$_;
					$tmp=~s/[\n\r]//g;
					if($tmp=~/^<(FamilyName|DomainMotif) pmid='$pmid' sid='[^ ]+' start='[^ ]+' end='[^ ]+' tax_id='[^ ]+'>(.+)<\/(FamilyName|DomainMotif)>/)
					{
						my $mention=lc($2);
						$mention=~s/[\W\-\_]//g;
						$FamilyDomain_hash{$mention}=1;
					}
				}
				close sa;
				
				my %gn_hash=();
				my %gn_len_hash=();
				open gn,"<$gn_extraction";
				while(<gn>)
				{
					my $tmp=$_;
					$tmp=~s/[\n\r]//g;
					if($tmp=~/<GN pmid='$pmid' tax_id='.+' entrezID='.+' score='.+'>(.+)<\/GN>/)
					{
						my $mentions=$1;
						$gn_hash{$mentions}=$tmp."\n";
						$gn_len_hash{$mentions}=length($mentions);
					}
					else
					{
						$tmpGN=$tmpGN.$tmp."\n";
					}
				}
				close gn;
				
				my %AbbMatched_hash=();
				my @mention_arr = reverse sort {$gn_len_hash{$a} <=> $gn_len_hash{$b}} keys %gn_len_hash;
				foreach my $mentions(@mention_arr)
				{
					my @split_mention=split(/\|/,$mentions);
					my %mention_hash=();
					foreach my $mention (@split_mention)
					{
						if($mention=~/^[ ]*(.+)[\+\*]$/)
						{
							$mention=$1;
							my $mention_org=$mention;
							$mention=~s/[\W\-\_]//g;
							if($FamilyDomain_hash{$mention}) # Family/Domain name
							{
								$gn_hash{$mentions}=~s/([>\|])($mention_org)[\+\*]([<\|])/$1$2$3/g;
							}
							elsif($ShortForm2LongForm_hash{$mention}) # Abbreviation
							{
								my $LForm=$ShortForm2LongForm_hash{$mention};
								my $match=0;
								foreach  my $mention_tmp (@split_mention)
								{
									$mention_tmp=~s/[\+\*]$//g;
									$mention_tmp=" ".lc($mention_tmp)." ";
									$mention_tmp=~s/[\W\-\_]+/ /g;
									$mention_tmp=~s/([0-9])([A-Za-z])/$1 $2/g;
									$mention_tmp=~s/([A-Za-z])([0-9])/$1 $2/g;
									$mention_tmp=~s/([A-Z])([a-z])/$1 $2/g;
									$mention_tmp=~s/([a-z])([A-Z])/$1 $2/g;
									$LForm=" ".lc($LForm)." ";
									$LForm=~s/[\W\-\_]+/ /g;
									$LForm=~s/([0-9])([A-Za-z])/$1 $2/g;
									$LForm=~s/([A-Za-z])([0-9])/$1 $2/g;
									$LForm=~s/([A-Z])([a-z])/$1 $2/g;
									$LForm=~s/([a-z])([A-Z])/$1 $2/g;
									$LForm=~s/[ ]+/ /g;
									my $LForm_formatch=$LForm;
									$LForm_formatch=~s/ /\[\\W\\\-\\\_\]\*/g;
									if($LForm=~/$mention_tmp/ || $mention_tmp=~/$LForm_formatch/)
									{
										$match=1;
									}
								}
								if($match==0 && (not exists $AbbMatched_hash{$mention}))
								{
									#print "Filtered by GNPostProcessing: $gn_hash{$mentions}\n\n"; #output filter mentions
									my $shortForm=$mention;
									$shortForm=" ".lc($shortForm)." ";
									$shortForm=~s/[\W\-\_]+/ /g;
									$shortForm=~s/([A-Za-z])(a )$/$1 $2/g;
									$shortForm=~s/([0-9])([A-Za-z])/$1 $2/g;
									$shortForm=~s/([A-Za-z])([0-9])/$1 $2/g;
									$shortForm=~s/([A-Z])([a-z])/$1 $2/g;
									$shortForm=~s/([a-z])([A-Z])/$1 $2/g;
									$shortForm=~s/[ ]+/ /g;
									$shortForm=~s/ /\[\\W\\\-\\\_\]\*/g;
									$gn_hash{$mentions}=~s/([>\|])($shortForm)[\+\*]([<\|])/$1$2$3/g;
									$gn_hash{$mentions}=~s/([>\|])($shortForm)[\+\*]([<\|])/$1$2$3/g;
								}
								else
								{
									$AbbMatched_hash{$mention}=1;
									$LForm=~s/ /\[\\W\\\-\\\_\]\*/g;
									my $found=0;
									foreach my $tmp_mention(keys %gn_hash)
									{
										if($gn_hash{$tmp_mention}=~/([>\|])($LForm)[\+\*]([<\|])/)
										{
											$found=1;
										}
									}
									if($found == 0)
									{
										$gn_hash{$mentions}=~s/([>\|])($LForm)([<\|])/$1$2\+$3/g;
									}
								}
							}
						}
					}
					
				}
				
				open gn,">$gn_extraction";
				print gn $tmpGN;
				foreach my $mentions (keys %gn_hash)
				{
					print gn $gn_hash{$mentions};
				}
				close gn;
			}
			
			%ShortForm2LongForm_hash=();
			%ShortForm2Score_hash=();
		}
	}
	close Abb;
	return 1;
}

sub Disambiguation_GN
{
	my ($finalCNT) = @_[0];
	my %Tmp_CandidateEntity_hash=();
	my %redundance_hash=();
	foreach my $each_entrez_id1(keys %CandidateEntity_token_hash)
	{
		my $tax_1=$CandidateEntity_tax_hash{$each_entrez_id1};
		my $termSet1=$CandidateEntity_token_hash{$each_entrez_id1};
		my @split_termSet1 = split(/\|/,$termSet1);
		my $redundanceFlag =0; #0=false ; 1=true
		if(not exists $redundance_hash{$termSet1})
		{
			foreach my $each_entrez_id2(keys %CandidateEntity_token_hash)
			{
				my $tax_2=$CandidateEntity_tax_hash{$each_entrez_id2};
				if($each_entrez_id1 != $each_entrez_id2 && $tax_1 eq $tax_2)
				{
					my $termSet2=$CandidateEntity_token_hash{$each_entrez_id2};
					if( $termSet1 ne $termSet2 )
					{
						if(not exists $redundance_hash{$termSet2})
						{	
							my @split_termSet2 = split(/\|/,$termSet2);
							my $count=0;
							my %split_termSet1_hash=();
							my %split_termSet2_hash=();
							my $count_split_termSet1=0;
							my $count_split_termSet2=0;
							foreach my $t(@split_termSet1){	$split_termSet1_hash{$t}=1;	}
							foreach my $t(@split_termSet2){	$split_termSet2_hash{$t}=1;	}
							foreach my $t(keys %split_termSet1_hash){$count_split_termSet1++;}
							foreach my $t(keys %split_termSet2_hash){$count_split_termSet2++;}
							foreach my $temp1(keys %split_termSet1_hash)
							{
								foreach my $temp2(keys %split_termSet2_hash)
								{
									if($temp1 eq $temp2)
									{
										$count++;
										last;
									}
								}
							}
							if($count==$count_split_termSet1)
							{
								$redundance_hash{$termSet1}=0;
								$redundanceFlag=1;
								last;
							}
							elsif($count==$count_split_termSet2)
							{
								$redundance_hash{$termSet2}=0;
							}
						}
					}
					
				}
			}
			if($redundanceFlag==0)
			{
				if(exists $Tmp_CandidateEntity_hash{$termSet1})
				{
					$Tmp_CandidateEntity_hash{$termSet1}=$Tmp_CandidateEntity_hash{$termSet1}."\t".$each_entrez_id1;
				}
				else
				{
					$Tmp_CandidateEntity_hash{$termSet1}=$each_entrez_id1;
				}
			}
		}
	}	

	foreach my $termSet (keys %Tmp_CandidateEntity_hash) 
	{
		my @split_entrez_id = split ("\t",$Tmp_CandidateEntity_hash{$termSet});
		my %tmp_CandidateEntity_token_hash=();
		my %tmp_CandidateEntity_weight_hash=();
		my %tmp_CandidateEntity_rank_hash=();
		foreach my $entrez_id_tmp ( @split_entrez_id )
		{
			$tmp_CandidateEntity_token_hash{$entrez_id_tmp}=$CandidateEntity_token_hash{$entrez_id_tmp};
			$tmp_CandidateEntity_weight_hash{$entrez_id_tmp}=$CandidateEntity_weight_hash{$entrez_id_tmp};
			$tmp_CandidateEntity_rank_hash{$entrez_id_tmp}=$CandidateEntity_weight_hash{$entrez_id_tmp}*10000000000000-$entrez_id_tmp;
			$tmp_CandidateEntity_tax_hash{$entrez_id_tmp}=$CandidateEntity_tax_hash{$entrez_id_tmp};
		}
		my @rank = reverse sort {$tmp_CandidateEntity_rank_hash{$a} <=> $tmp_CandidateEntity_rank_hash{$b}} keys %tmp_CandidateEntity_rank_hash;
		for($i=0;$i<$finalCNT;$i++)
		{
			$Final_CandidateEntity_token_hash{$rank[$i]}=$tmp_CandidateEntity_token_hash{$rank[$i]};
			$Final_CandidateEntity_weight_hash{$rank[$i]}=$tmp_CandidateEntity_weight_hash{$rank[$i]};
			$Final_CandidateEntity_tax_hash{$rank[$i]}=$tmp_CandidateEntity_tax_hash{$rank[$i]};
		}
	}
	return 1;
}

sub Gene_Normalization
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $id_extraction="tmp/".$filename.".id.xml";
	my $sa_extraction="tmp/".$filename.".sa.xml";
	my $table_extraction="tmp/".$filename.".table.xml";
	my $gn_extraction="tmp/".$filename.".gn.xml";
	$finalSSThreshold=0;
	$simiarlarscorethreshold=0;
	$specialScore=1000;
	$offset=500;
	$normal=1350;
	$base=0.2;
	my $finalCNT=1;
		
	#Initial
	#:Roman2Arabian_hash
	#:IDF_Entity_hash
	#:IDF_Token_hash
	{
		@Pmidlist=();
		my $tmppmid="";
		open sentence,"<$sentence_extraction";
		while(<sentence>)
		{
			my $tmp=$_;
			if($tmp=~/^<TEXT pmid='(.+)' sid='.+'>.+<\/TEXT>/)
			{
				if($tmppmid ne $1)
				{
					push(@Pmidlist,$1);
					$tmppmid=$1;
				}
			}
		}
		close sentence;
		
		%Roman2Arabian_hash=('ii' => 2, 'iii' => 3, 'iv' => 4, 'vi' => 6, 'vii' => 7, 'viii' => 8, 'ix' => 9, 'xi' => 11, 'xii' => 12, 'xiii' => 13, 'II' => 2, 'III' => 3, 'IV' => 4, 'VI' => 6, 'VII' => 7, 'VIII' => 8, 'IX' => 9, 'XI' => 11, 'XII' => 12, 'XIII' => 13);
	
		open sa,"<$sa_extraction";#IDF extraction	
		while(<sa>)
		{
			my $tmp=$_;
			if($tmp=~/tax_id=\'([0-9]+)\'/)
			{
				$TargetTax_hash{$1}=1;
			}
		}
		close sa;
		foreach my $tax(keys %TargetTax_hash) # IDF_Entity_hash / IDF_Token_hash
		{
			open input,"<".$dictionary."/DF_Entity/Tax_".$tax.".txt";
			my $sum=<input>;
			$sum=~s/[\n\t\r]//g;
			while(<input>)
			{
				my $tmp=$_;
				if($tmp=~/^(.+)\|([0-9]+)$/)
				{
					my $tmp1=$1;
					my $tmp2=$2;
					$tmp1=~s/[\n\t\r]//g;
					$tmp2=~s/[\n\t\r]//g;
					$IDF_Entity_hash{$tmp1}=log10($sum/$tmp2);
				}
			} 
			close input;
			open input,"<".$dictionary."/DF_Token/Tax_".$tax.".txt";
			$sum=<input>;
			$sum=~s/[\n\t\r]//g;
			while(<input>)
			{
				my $tmp=$_;
				if($tmp=~/^(.+)\|([0-9]+)$/)
				{
					my $tmp1=$1;
					my $tmp2=$2;
					$tmp1=~s/[\n\t\r]//g;
					$tmp2=~s/[\n\t\r]//g;
					$IDF_Token_hash{$tmp1}=log10($sum/$tmp2);
				}
			} 
			close input;
		}	
	}
	open output,">$gn_extraction";close output;
	foreach my $TargetPmid (@Pmidlist)
	{
		%CandidateEntity_token_hash=();
		%CandidateEntity_weight_hash=();
		%CandidateEntity_tax_hash=();
		%Final_CandidateEntity_token_hash=();
		%Final_CandidateEntity_weight_hash=();
		%Final_CandidateEntity_tax_hash=();
		
		%Entity_hash=(); 
		@Token_array=(); 
		@Entity_array=(); 
		Tokenization($TargetPmid,$sa_extraction);
		Tokenization_FamilyDomain($TargetPmid,$sa_extraction);
		
		#ExactMatch::OfficialName&ModifiedName
		foreach my $entity (keys %Entity_hash)
		{
			my @TAXs = split("@",$Entity_hash{$entity});
			foreach my $TAX (@TAXs)
			{ 
				ExactMatch($dictionary,$entity,$TAX);
			}
		}
		
		#Partial Match
		foreach my $entity (keys %Entity_hash)
		{
			my @TAXs = split("@",$Entity_hash{$entity});
			foreach my $TAX (@TAXs)
			{ 
				PartialMatch($dictionary,$entity,$TAX);
			}
		}
		
		#Disambiguation
		Disambiguation_GN($finalCNT);
	
		#Output
		open output,">>$gn_extraction";
		foreach $each_entrez_id(keys %Final_CandidateEntity_token_hash)
		{
			my $weight=($Final_CandidateEntity_weight_hash{$each_entrez_id}-$offset)/$normal+$base;
			if($weight>1){$weight=1;}
			my $entrez_id=$each_entrez_id;
			$each_entrez_id=~s/[\s\W\t\n\r]//g;
			printf output "<GN ".
							"pmid='".$TargetPmid."' ".
							"tax_id='".$Final_CandidateEntity_tax_hash{$entrez_id}."' ".
							"entrezID='".$each_entrez_id."' ".
							"score='%5.4f'".
							">".$Final_CandidateEntity_token_hash{$entrez_id}."<\/GN>\n",$weight;
		}
		close output;

		#ID Match
		ID_Match($dictionary,$TargetPmid,$table_extraction,$id_extraction,$gn_extraction);
		
		#GNPostProcessing
		GNPostProcessing($filename,$TargetPmid,$sa_extraction,$gn_extraction);
	}
	
	return 1;
}

sub Gene_Mapping
{
	my ($output)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $gmention_extraction="tmp/".$filename.".gmention.xml";
	my $sa_extraction="tmp/".$filename.".sa.xml";
	my $gn_extraction="tmp/".$filename.".gn.xml";
	my $ga_extraction=$output."/".$filename.".ga.xml";
	
	my $family_RegEx="(?:1(?:2otetradecanoylphorbol13acetate|4kda|acylglycerol3phosphateoacyltransferases)|2(?:hscdc2|mdabrca2containing)|5ht(?:serotonin)?|8bromocamp|90kda|a(?:a(?:rs[12]|tp)|b(?:c[abcdefg]|h(?:ydrolasedomaincontaining|d))|c(?:e(?:tyl(?:atedhistoneh4|c(?:holinereceptorrelated|oacarboxylase))|r)|kr|ot|ti(?:nrelatedprotein23|vatedgrowth)|ylcoa(?:synthetase|thioesterases)|s)|d(?:am(?:metallopeptidase(?:domaincontaining|swiththrombospondintype1)|ts)?|cy|en(?:ovirale1a|yl(?:ate(?:cycl|kin)ases|ylcyclase))|hesion|ipo(?:cytespecifictranscription|r)|o(?:ii|ra)|pribosylation(?:factor(?:gtpaseactivating|slike))?|r[abp]|vancedglycosylationendproductsof|[hp])|g(?:es|pat|tr|[eo])|if|k(?:ap|inaseanchor|rs|r)|l(?:coholdehydrogenases|d(?:ehydedehydrogenases|oketoreductases?|h)|k(?:alineceramidase|ylation(?:repairhomologs)?|b)|ox|pha(?:2aa?drenoceptor|classglutathionetransferase)s|uy)|m(?:ino(?:a(?:cyltrnasynthetases|ndcarboxyterminal)|terminal)|pdependentproteinkinase|yloidprecursor)|n(?:ap(?:hasepromoting|c)|k(?:rd|yrinrepeatdomaincontaining)|nexins|p32(?:acidicnuclearphospho)?|xa|o)|p(?:o(?:bec|lipo(?:proteinbmrnaeditingenzymes)?)|[1p])|q(?:uaporinwaterchannel|p)|r(?:achidonatelipoxygenases|fgap|gonautepiwi|hg(?:ap|ef)|m(?:adillorepeatcontaining|c)|pc|ylsulfatase|[fls])|sic|t(?:axins|p(?:ase(?:alphachain|s)?|bindingcassettetransporters|dependent(?:dnaunwindingenzym|proteas)e)|xn|[fp])|uto(?:crineandparacrinegrowth|somaldominantretinitispigmentosa)|vpr|xonemaldyneins|k)|b(?:3ga?t|4gt|65|a(?:culoviraliaprepeatcontaining|s(?:eexcisionrepair|ic(?:helixloophelix|leucinezipper|regionleucinezipper)))|c(?:g(?:1tobcg88|72|84)|kdh)|dkr|e(?:nd(?:omaincontaining)?|st|ta(?:1(?:3glucuronyl|6glcnac)transferases|3glycosyltransferases|4glycosyltransferases|a(?:drenergic|rrestin(?:greenfluorescent)?|r)|subunitofg)|r)|hlh|i(?:ogenesisoflysosomalorganellescomplex1|rc)|lo(?:c1s|odgroup(?:antigens)?)|mps?|onemorphogenetic|pif(?:oldcontaining)?|r(?:anchedchainalphaketoaciddehydrogenase|c(?:ax|t)|ic(?:hosdomaincontaining|d)|s)|tb(?:pozdomaincontaining|d)|zip(?:like)?|2)|c(?:1set|2(?:chtypezincfinger|h2zincfinger|set)|a(?:2(?:binding|calmodulin(?:dependentproteinkin|stimulatedproteinphosphat)ase)|at|cn1?|dherins|l(?:c(?:i(?:neurin|um(?:activatedpotassiumchannels?|bindings100|channel))|r)|modulincalmodulinlike|poninhomology2)|mkii|n(?:didatetumorsuppressor|nabinoid)|r(?:bo(?:nicanhydrases|xy(?:lesterases|terminal))|cinoembryonicantigencearelated)|s(?:pases|scaffolding|[prs])|t(?:alytic|hepsins|sper)|n)|bc|c(?:chemokines?|knr|[lr])|d(?:34|43intracellular|c(?:likekinases|2)|gshironsulfurdomaincontaining|hr|k(?:inhibitors?|s)|molecules|[3k])|e(?:d3|ll(?:deathinducingdff45likeeffector|surface)|ntralandperipheralcannabinoid|rs|s)|ftr|gmpgatedcationchannels|h(?:a(?:rgedmultivesicularbody|p)|chd|emo(?:attractant|kine(?:beta|ligands)?)|mp|ol(?:esterolesterase|inergic)|r(?:omatinmodifyingenzymes|[mn])|r)|i(?:de|sd)|kis?|l(?:audins|cn|dn|ec|ic|mf|k)|n(?:ar|[gr])|o(?:iledcoilhelixcoiledcoilhelixdomaincontaining|l(?:ec|l(?:agen(?:oustrimerictypeiimembrane|vi|s)?|ectins)|pg|6)|m(?:i(?:ii|[iv])|p(?:lement(?:system)?|onentsofoligomericgolgi)|i)|ro(?:nins)?|g)|r(?:eatinekinase|hr)|spg|t(?:ypelectin(?:domaincontaining)?|[dfs])|ut|x(?:3cr|c(?:chemokines?|r))|y(?:b5ap|cli(?:cnucleotidegatedchannel|ndependent(?:celldivision)?kinase)s|p4(?:50)?|sticfibrosis|to(?:chrome(?:p450s?|b)|kine|plasmic(?:dyneins)?|toxiclymphocytematuration)|[bp])|[9adf])|d(?:bdh|caf|d(?:b1andcul4associated|x)|e(?:a[dh]boxes|f(?:ensins(?:alph|bet)a|[ab])|nn(?:madddomaincontaining|d)|xamethasone|x)|hx|i(?:hydropyridine|s(?:coidin|hevelledhomologs))|n(?:a(?:b(?:inding|ound(?:tfiiatbp|vp16activation))|dependentatpase|ligases|polymerase(?:ii|s)?|j)|lz)|o(?:l(?:ichyldmannosylphosphatedependentmannosyltransferases|pm)|ublezincfinger)|rd|u(?:alspecificityphosphatase|sp[acmpqst])|vl|y(?:slexiacandidate|n)|[1n])|e(?:6binding|c(?:mpg|d)|dnr|f(?:calciumbinding|hand(?:domaincontaining|motif)?|n)|l(?:mo|ongatoracetyltransferase|p)|mid(?:omaincontaining)?|n(?:ac|do(?:genousligands|lig|plasmicreticulum(?:golgi)?)|gulfmentandcellmotility|hancedrnaithreeprimemrnaexonucleases)|p(?:hr(?:elated|ins)|id(?:ermal(?:cellenvelope|growth(?:factorlike)?)|idymis)|h)|ri|st(?:rogen(?:responsive)?)?|t[fs]|x(?:ostosinglycosyltransferase|portins|tracellularcytokinebinding|t)|[236cr])|f(?:1(?:atp|f0atpsynthase)|2r|a(?:bp|ds|nc(?:onianemiacomplementationgroups)?|t(?:hd|tyacid(?:binding|desaturases|hydroxylasedomaincontaining)|p)|k)|b(?:ln|oxes(?:leucinerichrepeats|other|wd40)|x[low])|erm(?:itins|t)?|far|ib(?:r(?:i(?:llarcollagentypev|nogencdomaincontaining)|onectintypeiii(?:domaincontaining|like))|ulins|c)|k506|lywch|mo|n3|o(?:caladhesionkinase|r(?:khead(?:boxes)?|skolin)|x)|pr|reac|u(?:cosyltransferases|t)|zd|c)|g(?:1cyclin|a(?:b(?:br|ra|[ar])|l(?:4dnabinding|beta13glcnacbeta16galnac|nactransferase|r)|mmaglutamyltransferases|ta(?:zincfingerdomaincontaining|d)|p)|beta|cgr|eneraltranscription|fp|gt|hsr|i(?:coupled|map)|l(?:ra|t[1268]|uta(?:m(?:ate|inerich)|thionestransferases?)|y(?:c(?:erolkinases|o(?:genphosphorylases|syltransferase)|o)|rs|r))|nr(?:hr|ps|p)|olgi|p(?:atch(?:domaincontaining)?|cr(?:ao|bo|co|gproteincoupledreceptorkinase|s)?|nloopgtpases|roteincoupled(?:chemoattractant)?|[cn])|qcoupled|r(?:eenfluorescent|owingiceced3|[im])|s(?:coupled|giandgqproteincoupled|tk|t)|t(?:p(?:ase(?:activating|regulator|simap)?|binding)|shr|[fp])|uan(?:inenucleotidereleasing|osine(?:nucleotidebindingproteincoupled|triphosphategtpbinding))|[jkr])|h(?:aus(?:augminlike)?|c(?:rt|a)r|e(?:atshock|lixloophelix|terogeneousnuclearribonucleo)|fshrsextracellular|hr23|i(?:ghmobilitygroup|st(?:o(?:compatibility|nes))?)|l(?:a(?:classii|dr|g)|[ah])|mg(?:a1p|b(?:ox|p)|coasynth(?:ase|esis)|np|x)?|n(?:rnps?|f)|o(?:meobox(?:es)?|xl)|rh|s(?:as|d17b|p(?:70|90|[bc]))|tr(?:gp|3)?|v(?:cn|mat)|y(?:al(?:ectan|uronatebinding)|dro(?:phobic|xyaminoaciddehydratases)))|i(?:ap|ce|f(?:f[123456o]|[nt])|g(?:ho|ko|l(?:ike|o)|[adeghjkl])|k(?:a(?:ppab(?:kinases?)?|roszincfingers)|ks|zf|k)|l1[2r]|m(?:muno(?:deficiencyvirustype1|globulin(?:super|s)?)|p(?:dh|ortins))|n(?:hibitor(?:ofapoptosis|yglycine)|itiation|o80|t(?:e(?:grin(?:binding|s)|r(?:ferons|leukin(?:12|il12|sandinterleukin)|mediatefilament(?:familyorphans|stype(?:i(?:ii|[iv])|vi|[iv]))))|ra(?:cellularmannosespecificlectin|flagellartransporthomologs)))|on(?:channels|otropicglutamate)|po|set|t(?:pr|g)|[gl])|j(?:nks?|unkinases)|k(?:a(?:inateselectivechannel|llikreins|nk|t)|c(?:hannel(?:beta|s)|n[jk]|[an])|dm|e(?:lchlike|ratin(?:associated)?)|i(?:llercell(?:immunoglobul|lect)inlike|n(?:esin(?:likecoiledcoil|s)|in)|[fr])|l(?:hl|[fkr])|m(?:sks|t)|nmotifandankyrinrepeatdomaincontaining|r(?:ab(?:a(?:andkrabbboxes|box)|b)?|ppelassociatedbox|tap|uppelliketranscription)|u(?:autoantigen|nitz(?:typeserineproteaseinhibitor)?)|v(?:beta|channel)|[uv])|l(?:a(?:ir|minins|r(?:ibonucleoproteindomaincontaining|p)|te(?:cornifiedenvelopes|nttransforminggrowthfactorbetabinding)|m)|c[en]|d(?:associated|lr|l)|e(?:ctinsgalactosidebinding|rks|u(?:cinezippertranscription|kocyte(?:associatedig|immunoglobulin)like))|gals|i(?:gand(?:gatedionchannel|softheephrelatedkinase)s|lr|p(?:iddropletassociated|o(?:calin|xygenase)s)|ssencephaly|m)|ncrna|o(?:n(?:gnoncodingrnas|protease)|wdensitylipo)|par|ribosomal|t(?:bp|nr)|y(?:6antigen|mpho(?:cyteantigen|idnuclear)|rm(?:otifcontaining)?|so(?:phospholip(?:ases?|ids)|somalcysteineproteinase))|1)|m(?:a(?:drelated|estroheatlikerepeatcontaining|nnosespecificmembranelectin|p(?:2k|3k|4k|k(?:inase|s)|k)|r(?:ch(?:membraneassociatedringfingers)?|sardnabinding)|xdimerization|[dp])|c(?:dh|hr|nr|m)|e(?:f2|tallothioneins|x3(?:homologs)?)|gst|i(?:s(?:crna|matchrepair)|to(?:af|chondrialr(?:espiratorychain(?:complexassembly)?|ibosomal)|genactivatedproteinkinase(?:cascade|s)?)|r)|lnr|mr|o(?:bkinaseactivators|lecularweight|tilitystimulatingexophosphodiesterase|b)|r(?:oh|p[los])|sg|t(?:nr|trna)|u(?:cin(?:typeolinkedproteinglycosylation|s)|rinehomeoboxcontaining|c)|xd|y(?:bp|elocytomatosisviralonco|hii|o(?:cyteenhancer(?:factor2)?|i(?:ii|x)|sin(?:binding|s)|vii?|x(?:ix|vi(?:ii)?|v)|[ivx])|[cl])|[rt])|n(?:a(?:d(?:ependentinorganicpicotransporters|phoxidase|h)|gammaaminobutyricacidtransporter|ktransportingatpaseinteracting|l(?:cn|phaacetyltransferase)|p(?:hosphate|o4)cotransporter|a)|bpf|c(?:rnas|k)|e(?:cab|rvegrowth|trins|uro(?:blastomabreakpoint|nspecificsynapticvesicleassociatedphospho|tr(?:ansmittertransporter|ophic)))|f(?:at|e2(?:binding|relatedtranscription|transcription)|il2luciferase|kappab)|gf|k(?:ain|l)|lr|mur|o(?:n(?:glutamateleucinearginine|receptortyrosinekinases?)|p2sundomaincontaining|radrenalinetransportinhibitors)|p(?:bwr|ffr|sr|xy|yr|t)|sun|t(?:erminalefhandcalciumbinding|sr|n)|u(?:cle(?:ar(?:factorofactivatedtcells|hormone)|otidebindingdomainandleucinerichrepeatcontaining)|d(?:ixmotifcontaining|t))|r)|o(?:7tm|glycan(?:core2|s)|megahydroxylase|p[nr]|r(?:1[01234]|5[1256]|ai(?:calciumreleaseactivatedcalciummodulators)?|phangproteincoupled|[123456789])|sbp|tud(?:omaincontaining)?|x(?:idationofcoproporphyrinogeniiitoprotoporphyrinogenix|ysterolbinding))|p(?:2(?:r[xy]|[ry])|450|a(?:di|iredbox(?:es)?|nx|paintypecysteineproteinase|r(?:a(?:neoplasticmaantigen|oxonase)s|kinsondisease|vins|[12kpsv])|t(?:atinlikephospholipasedomaincontaining|[ep])|[krx])|c(?:dh[cn]|gf)|d[eiz]|e(?:l(?:linohomologs|i)|ptidylargininedeiminases|rilipins)|h(?:actr|osph(?:at(?:ase(?:andactinregulators)?|idylinositol(?:3kinase|glycananchorbiosynthesis))|o(?:diesterase(?:inucleotidepyrophosphatase|s)|inositide3kinases?))|f)|i(?:rc|g)|k[cn]|l(?:asmamembraneassociatedgtpbinding|e(?:ckstrinhomologyphdomaincontaining|kh|xins)|in|xn)|n(?:pl|m)a|o(?:l(?:y(?:adpribosepolymerases?|combgroupringfingers|peptidenacetylgalactosaminyltransferases|unsaturatedfattyacids)|r)|mbecelldivisioncycle2cdc2likekinases|t(?:assiumchannels|eankyrindomaincontaining|e)|[lnu])|p(?:125|p[12346]r|[mp])|r(?:enyltransferasealphasubunitrepeatcontaining|mt|o(?:gesterone|kr|linerich(?:transmembrane)?|t(?:e(?:as(?:eactivated(?:kinase)?|omeprosomemacropain)|in(?:argininemethyltransferases|disulfideisomerases|kinasec|tyrosine(?:kinase|phosphatases?))|oglycans)|oonco)|x)|rt|ss|d)|s(?:betag|eudoautosomalregions|m)|t(?:af?r|gr|hnr|p(?:ase|[23enr])|s(?:1containing|[12])|k)|u(?:rinergic|tative(?:calciumbinding|transmembrane|zincfinge))|v92|yg|r)|qpgy|r(?:a(?:b(?:gdi|memberrasonco|relatedgtpbinding)|cs|d23|mp|pidlygrowingfamilyofeukaryotictranscriptionalregulators|srelated(?:c3botulinumtoxinsubstrate)?|[bcrs])|bm|e(?:ceptor(?:accessory|gproteincoupledactivitymodifying|proteintyrosinekinases?|t(?:ransporter|yrosinekinase))|ep|gulatorsofgproteinsignaling|plicationfactorc|tin(?:alspecific|oicacid))|f(?:apr|c)|gs|ho(?:g(?:tpaseactivating|uaninenucleotideexchange)|a)?|i(?:bonucleasesrnasea|ng(?:finger|h2zincfinger|typec3hc4zincfingers)?|h)|n(?:a(?:bindingmotifrrmcontaining|helicase|p(?:olymerase(?:ii(?:holoenzyme|poliiholoenzyme|transcripts))?|seudouridylatesynthasedomaincontaining)|se)|[fy])|p(?:tks?|usd|[ls])|rna|tp|vnr|xfp|yr)|s(?:1(?:00(?:calciumbinding)?|pr)|a(?:ll|m(?:andsh3domaincontaining|d)|pks?|sh)|c(?:a(?:mp|rna|d)|gb|nn|rna|n)|d(?:r(?:c[123]|[ae])|c)|e(?:cret(?:ed(?:frizzledrelated|photo)|o(?:globins|ry(?:carriermembran|tumornecrosisfactorinducibl)e))|ma(?:phorins)?|pt(?:ins)?|r(?:ine(?:a(?:ndtyrosinephosphatase|rgininerichsplicing)|orcysteinepeptidaseinhibitors|peptidase(?:inhibitorskazaltype|s)|threoninep(?:hosphatase|roteinkina)s)|pin)|ventransmembranespanning)|f(?:rp|xn)|gs[mt]|h(?:2d(?:omaincontaining)?|3binding|fm|isa(?:homologs)?|mt|ortchain(?:alcoholdehydrogenase|dehydrogenasereductasesuper)|[23])|i(?:al(?:icacidbindingiglikelectins|omucin|yltransferases)|deroflexins|g(?:lec|nal(?:regulatory|transducersandactivatorsoftranscription))|rp|x)|k(?:channels?|itranscriptionalcorepressors|or)|l(?:rr|c)|m(?:a(?:ds|llgproteinsignalingmodulators|d)|c)|n(?:a(?:ilhomologs|re|i)|or[ad]|rna|[fhx])|o(?:cs|diumchannels|lutecarriers|rtingnexins|[cx])|p(?:dy|e(?:cific(?:ityprotein)?transcription|edyhomologs|rmantigen)|ink|lithandsplitfootmalformation|t)|r(?:chomology(?:domain3|sh[23]|[23])?|ibosomal|sf|y(?:boxcontaining|relatedtranscription|sexdeterminingregionyboxes)|c)|s(?:tr|x)|t(?:3g|a(?:r(?:relatedlipidtransferstartdomaincontaining|d)|t)|er(?:ilealphamotifsamdomaincontaining|oid(?:thyroidhormonevitamin)?)|oreoperatedchannels|r(?:essactivatedproteinkinases?|ucturalmaintenanceofchromosomes))|u(?:bmaxillarygland|l(?:fotransferases(?:cytosolic|membranebound)|tm|t)|p(?:eroxidedismutase|pressorsofcytokinesignaling))|y(?:n(?:ap(?:sin|to(?:brevin|tagmins?))|taxin1)|t)|[cp])|t(?:47d|a(?:ar|cr|fs|le|s(?:1r|2r|nr|tereceptorstype[12])|ta|[ft])|b(?:oxes|p(?:associated|i)|x)|c(?:ell|tn|r)|drd|e(?:ashirtzincfingers|ctonic|tra(?:spanins|tricopeptidettcrepeatdomaincontaining))|fii[dh]|g(?:fbeta(?:inducedheteromericsignaling)?|m)|h(?:ap(?:domaincontaining)?|oc|rombin|o)|m(?:cc|prss)|n(?:f(?:alpha|binding|like|r(?:elatedligands?|sf|s)|sf|r)|rc|f)|p(?:a25|cn|ks|[km])|r(?:a(?:ffickingproteinparticle|ns(?:acting(?:cellular)?|criptionalactivationofclassii|forminggrowthfactorbeta|glutaminases|ientreceptorpotential|membrane(?:andcoiledcoildomaincontaining|serineprotease)?)|ppc|s)|bo|i(?:mp|nucleotidecagrepeatcontaining|p(?:aritemotifcontaining|lehelical)|m)|na(?:splicingendonuclease)?|opomyosins|pl|[bdgp])|s(?:en|hz|pan)|t(?:ll|c)|u(?:b(?:eroussclerosis|ulin(?:tyrosineligaselike|s))|dordomaincontaining|mor(?:necrosis|suppressor)|b)|y(?:pe(?:1tnf|i(?:i(?:c(?:a2dependentlectin|ollagenopath(?:ies|y))|keratins|ntegralmembrane)|kinase)|svandxiprocollagen|vcollagen)|rosine(?:kinases?|p(?:hosphatases?|roteinkinase))))|u(?:b(?:e[12]|iqui(?:lin|tin(?:conjugatingenzyme(?:se2)?|like(?:modifieractivatingenzymes)?|proteinligase(?:e3componentnrecognins)?|specificpeptidases))|ox(?:domaincontaining)?|qln|x(?:domaincontaining|n)|r)|dp(?:g(?:alnac|l(?:cnacpyrophosphorylase|uc(?:ose|uronosyltransferases)))|nacetyl(?:alphadgalactosamine|glucosaminepyrophosphorylases))|gt|pstreamstimulatoryactivityusaderivedco|ridine5diphosphoglucose|sp)|v(?:a(?:cuolaratpase|mp|nins|tp|ult(?:rnas)?)|dac|esic(?:leassociatedmembrane|ularmonoaminetransporter)|h1|ippacr|mat|n(?:1ru?|2r|n)|o(?:ltage(?:dependentan|gated)ionchannels|meronasal)|set|trnas)|w(?:a(?:ardenburgshahsyndrome|pfourdisulfidecoredomaincontaining|s(?:proteinhomologs|h))|dr(?:epeatdomaincontaining)?|fdc|inglesstypemmtvintegrationsites|nt|wc|w)|x(?:cr|linkedagammaglobulinemia|po)|y(?:ip[1f]|rnasroassociated)|z(?:acn|bed|c(?:2hc|3hc?|4h2|chc)|d(?:bf|hhc)|f(?:and|c3h1|hx|yve)|incfingers?|m(?:at|iz|y(?:nd|m))|n(?:fecontainingactive|hit|f)|onapellucidaglyco|ranb|swim|yg11(?:cellcycleregulator)?|zz|p)|g)";
	my $suffix_RegEx="(?:complex(?:es)?|do(?:amins|main)|fa(?:ctors?|mil(?:ies|y))|genes?|proteins?|subfamil(?:ies|y))";
	
	if(@Pmidlist eq 0)
	{
		my $tmppmid="";
		open sentence,"<$sentence_extraction";
		while(<sentence>)
		{
			my $tmp=$_;
			if($tmp=~/^<TEXT pmid='(.+)' sid='.+'>.+<\/TEXT>/)
			{
				if($tmppmid ne $1)
				{
					push(@Pmidlist,$1);
					$tmppmid=$1;
				}
			}
		}
		close sentence;
	}
	
	open ga,">$ga_extraction";close ga;
	foreach my $TargetPmid (@Pmidlist)
	{
		my $pmid="";	
		my %abb_hash=();
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
			elsif($tmp=~/^  (.+)\|(.+)\|(.+)$/)
			{
				if($pmid eq $TargetPmid)
				{
					if($3>0.9)
					{
						my $ShortForm=lc($1);
						my $LongForm=lc($2);
						$LongForm=~s/[\W\-\_]//g;
						$ShortForm=~s/[\W\-\_]//g;
						$abb_hash{$ShortForm}=$LongForm;
						$abb_hash{$LongForm}=$ShortForm;
					}
				}
			}
		}
		close Abb;
		
		#=====
		#domain_mapping
		my %Domain_hash=();
		my %Family_hash=();
		my %Domain_GM_hash=();
		my %Domain_tkn_hash=();
		my %Family_GM_hash=();
		open DomainFamily,"<$gmention_extraction";
		while(<DomainFamily>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/^<DomainMotif pmid='$TargetPmid' sid='(.+)' start='(.+)' end='(.+)' score='.+'>(.+)<\/DomainMotif>$/)
			{
				my $sid=$1;
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $nospace_mention=lc($mention);
				$nospace_mention=~s/[\W\-\_]//g;
				$Domain_GM_hash{$nospace_mention}=$Domain_GM_hash{$nospace_mention}.$tmp."\n";	
				$Domain_tkn_hash{$nospace_mention}=lc($mention);		
			}
			elsif($tmp=~/^<FamilyName pmid='$TargetPmid' sid='(.+)' start='(.+)' end='(.+)' score='.+'>(.+)<\/FamilyName>$/)
			{
				my $sid=$1;
				my $start=$2;
				my $last=$3;
				my $mention=$4;
				my $nospace_mention=lc($mention);
				$nospace_mention=~s/[\W\-\_]//g;
				$Family_GM_hash{$nospace_mention}=$Family_GM_hash{$nospace_mention}.$tmp."\n";
			}
		}
		close DomainFamily;
		
		open gn,"<$gn_extraction";
		my %gn_hash=();
		my %gn_length_hash=();
		my %gn_entrez_hash=();
		my %entrez_hash=();
		while(<gn>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<GN pmid='$TargetPmid' tax_id='(.+)' entrezID='(.+)' score='(.+)'>(.+)<\/GN>/)
			{
				my $tax=$1;
				my $entrezid=$2;
				my $weight_org=$3;
				$weight_org = round($weight_org , 4);
				my $mentions=$4;
				$entrez_hash{$entrezid}=1;
				if($mentions=~/[\*\+]/) #Family/Domain exact matches haven't */+
				{
					my @split_mentions=split(/\|/,$mentions);
					foreach $entity (@split_mentions)
					{
						my $weight=$weight_org;
						if($entity=~/^(.+)[\*\+]$/) #Entity inference
						{
							my $mention=lc($1);
							my $mention_tmp=$mention;
							$mention_tmp=~s/[\W\-\_]//g;
							if($mention_tmp!~/^($family_RegEx)($suffix_RegEx|)$/)
							{
								if((not exists $gn_hash{$mention."@".$tax}) || $gn_hash{$mention."@".$tax}<$weight)
								{
									$gn_length_hash{$mention."@".$tax}=length($mention);
									$gn_hash{$mention."@".$tax}=$weight;
									$gn_entrez_hash{$mention."@".$tax}=$entrezid;
								}
							}
						}
						
						my $mention=$entity;
						$mention=~s/[\W\-\_]//g;
						if(exists $Domain_GM_hash{$mention})
						{
							if($Domain_hash{$mention}!~/(^|,)$entrezid,/)
							{
								$Domain_hash{$mention}=$Domain_hash{$mention}.$entrezid.",";
								$Domain_hash{$abb_hash{$mention}}=$Domain_hash{$abb_hash{$mention}}.$entrezid.",";
							}
						}
						elsif(exists $Family_GM_hash{$mention})
						{
							if($Family_hash{$mention}!~/(^|,)$entrezid,/)
							{
								$Family_hash{$mention}=$Family_hash{$mention}.$entrezid.",";
								$Family_hash{$abb_hash{$mention}}=$Family_hash{$abb_hash{$mention}}.$entrezid.",";
							}
						}
					}
				}
			}
		}
		close gn;
		
		my %Domain_pssmid_hash=();
		open input,"<".$dictionary."/GeneDomain2DomainName.txt";
		while(<input>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/^([0-9]+)\t([0-9]+)\t(.+)$/)
			{
				if(exists $entrez_hash{$1})
				{
					my $geneid=$1;
					my $pssmid=$2;
					my $tokens=$3;
					my @split_token=split(" ",$tokens);
					foreach my $token(@split_token)
					{
						if($token!~/(domain|binding)/)
						{
							$Domain_pssmid_hash{$token}=$Domain_pssmid_hash{$token}.$pssmid."|";
						}
					}
				}
			}
		}
		close input;
		
		my %mention2pssmids_hash=();
		foreach my $nospace_mention(keys %Domain_tkn_hash)
		{
			my %tmp_hash=();
			my @tokens=split(/[\W\-\_]/,$Domain_tkn_hash{$nospace_mention});
			foreach my $tkn(@tokens)
			{
				if(exists $Domain_pssmid_hash{$tkn})
				{
					$Domain_pssmid_hash{$tkn}=~s/\|$//g;
					my @pssmids=split(/\|/,$Domain_pssmid_hash{$tkn});
					foreach my $pmssid(@pssmids)
					{
						$tmp_hash{$pmssid}++;
					}
				}
			}
			my @tmp = reverse sort {$tmp_hash{$a} <=> $tmp_hash{$b}} keys %tmp_hash;
			$mention2pssmids_hash{$nospace_mention}=$tmp[0];
		}
		
		my %annotation_rank_hash=();
		my %annotation_hash=();
		my %tmp_hash=();
		open sa,"<$sa_extraction";
		while(<sa>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/^<Gene pmid='$TargetPmid' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.+)' org='(.+)'>(.+)<\/Gene>$/
			 || $tmp=~/^<Gene pmid='$TargetPmid' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.+)'>((.+))<\/Gene>$/ 
			)
			{
				my $sid=$1;
				my $start_org=$2;
				my $last_org=$3;
				my $tax=$4;
				my $mention_org=$5;
				my $mention_mod=$6;
				my $sid_num=0;
				if($sid=~/([0-9]+)$/)
				{
					$sid_num=$1;
				}
				my $suffixp="N";
				if($mention_mod=~/^(.+)p$/){$suffixp="Y";}
				$mention=lc($mention_mod);
				
				my @gn_length_rank = reverse sort {$gn_length_hash{$a} <=> $gn_length_hash{$b}} keys %gn_length_hash;
				foreach my $tmp (@gn_length_rank)
				{
					my $start=$start_org;
					my $last=$last_org;
					if($tmp=~/^(.+)\@(.+)$/)
					{
						my $mention_gn=$1;
						my $tax_gn=$2;
						if($tax eq $tax_gn)
						{
							my $mention_gn_forcompare=$mention_gn;
							my $mention_gn_forcompare2=$mention_gn;
							$mention_gn_forcompare2=~s/[\W\_\-]/\.\*/g;
							$mention_gn_forcompare=~s/[\W\_\-]//g;
							my $mention_tmp=$mention;
							my $conj=0;
							if($mention_tmp=~/(,| and | or | to |\/)/){$conj=1;}
							$mention_tmp=~s/[\W\_\-]//g;
							$mention_gn_forcompare=~s/1/\(1\|i\)/g;
							$mention_gn_forcompare=~s/2/\(2\|ii\)/g;
							$mention_gn_forcompare=~s/3/\(3\|ii\)/g;
							$mention_gn_forcompare=~s/4/\(4\|iv\)/g;
							$mention_gn_forcompare=~s/a/\(a\|alpha\)/g;
							$mention_gn_forcompare=~s/b/\(b\|beta\)/g;
							my $paran=0;
							if($mention_gn_forcompare=~/\(.*\)/){$paran=1;}
							if ($mention_tmp=~/^($mention_gn_forcompare([p]{0,1}))$/)
							{
								my $str1=$1;
								my $sufp=$2;if($paran==1){$sufp=$3;}
								if($sufp eq "p" && $suffixp eq "N")	
								{
									#print "$mention_gn_forcompare\n";
								}
								else
								{
									my $weight=$gn_hash{$mention_gn."@".$tax};
									my $entrezid=$gn_entrez_hash{$mention_gn."@".$tax};
									if(exists $annotation_hash{$sid."|".$start."|".$last})
									{
										$annotation_hash{$sid."|".$start."|".$last}=$annotation_hash{$sid."|".$start."|".$last}."\n<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_org<\/Gene>";
									}
									else
									{
										$annotation_hash{$sid."|".$start."|".$last}="<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_org<\/Gene>";
									}
									$annotation_rank_hash{$sid."|".$start."|".$last}=$sid_num*1000+$start;
								}
							}
							elsif ($mention_tmp=~/^$mention_gn_forcompare2/ && $conj==1)
							{
								my $weight=$gn_hash{$mention_gn."@".$tax};
								my $entrezid=$gn_entrez_hash{$mention_gn."@".$tax};
								if(exists $annotation_hash{$sid."|".$start."|".$last})
								{
									$annotation_hash{$sid."|".$start."|".$last}=$annotation_hash{$sid."|".$start."|".$last}."\n<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_org<\/Gene>";
								}
								else
								{
									$annotation_hash{$sid."|".$start."|".$last}="<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_org<\/Gene>";
								}
								$annotation_rank_hash{$sid."|".$start."|".$last}=$sid_num*1000+$start;
							}
							else
							{
								my $mention_nospace=$mention_gn;
								$mention_nospace=~s/[\W\-\_]//g;
								$mention_gn_forcompare=$mention_gn;
								$mention_gn_forcompare=~s/[\W\_\-]/\[\\W\\\_\\\-\]\*/g;
								$mention_gn_forcompare=~s/1/\(1\|i\)/g;
								$mention_gn_forcompare=~s/2/\(2\|ii\)/g;
								$mention_gn_forcompare=~s/3/\(3\|iii\)/g;
								$mention_gn_forcompare=~s/4/\(4\|iv\)/g;
								$mention_gn_forcompare=~s/a/\(a\|alpha\)/g;
								$mention_gn_forcompare=~s/b/\(b\|beta\)/g;
								if ((not exists $Family_GM_hash{$mention_nospace}) &&
									(not exists $Domain_GM_hash{$mention_nospace}) &&
									$mention_org=~/^(.*)($mention_gn_forcompare[p]{0,1})(.*)$/i)
								{
									my $str1=$1;
									my $str2=$2;
									my $startL=length($str1)+length($str2);
									my $lengL=length($mention_org)-$startL;
									my $str3=substr($mention_org,$startL,$lengL);
									if($str2=~/p$/ && $suffixp eq "N")	{}
									else
									{
										if(($str1 eq "" || $str1=~/[\W\-\_]$/) && ($str3 eq "" || $str3=~/^[\W\-\_]/))
										{
											$start=$start+length($str1);
											my $mention_tmp=substr($mention_org,length($str1),length($str2));
											$last=$start+length($mention_tmp);
											my $weight=$gn_hash{$mention_gn."@".$tax};
											my $entrezid=$gn_entrez_hash{$mention_gn."@".$tax};
											if(exists $annotation_hash{$sid."|".$start."|".$last})
											{
												$annotation_hash{$sid."|".$start."|".$last}=$annotation_hash{$sid."|".$start."|".$last}."\n<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_tmp<\/Gene>";
											}
											else
											{
												$annotation_hash{$sid."|".$start."|".$last}="<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax' entrezID='$entrezid' score='$weight'>$mention_tmp<\/Gene>";
											}
											$annotation_rank_hash{$sid."|".$start."|".$last}=$sid_num*1000+$start;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		close sa;
		
		foreach my $annot1(keys %annotation_rank_hash)
		{
			my @S_annot1=split(/\|/,$annot1);
			foreach my $annot2(keys %annotation_rank_hash)
			{
				my @S_annot2=split(/\|/,$annot2);
				if(	$S_annot2[0] eq $S_annot1[0] && ((($S_annot2[1]<$S_annot1[1]) && ($S_annot2[2]>=$S_annot1[2])) || (($S_annot2[1]<=$S_annot1[1]) && ($S_annot2[2]>$S_annot1[2]))))
				{
					delete $annotation_rank_hash{$annot1};
				}
			}
		}
		
		foreach my $tmp (keys %annotation_hash)
		{
			if($annotation_hash{$tmp}=~/\n/)
			{
				my %entrezID_set=();
				my @split_annotation=split(/\n/,$annotation_hash{$tmp});
				my $sid="";
				my $start="";
				my $last="";
				my $tax_id="";
				my $weight=1;
				my $mention="";
				foreach my $anno(@split_annotation)
				{
					if($anno=~/<Gene pmid='$TargetPmid' sid='(.+)' start='(.+)' end='(.+)' tax_id='(.+)' entrezID='(.+)' score='(.+)'>(.+)<\/Gene>/)
					{
						$sid=$1;
						$start=$2;
						$last=$3;
						$tax_id=$4;
						my $entrezID=$5;
						$mention=$7;
						if($weight>$6){$weight=$6;}
						if($entrezID_set{$sid."\t".$start."\t".$last}=~/,$entrezID,/){}
						elsif(exists $entrezID_set{$sid."\t".$start."\t".$last})
						{
							$entrezID_set{$sid."\t".$start."\t".$last}=$entrezID_set{$sid."\t".$start."\t".$last}.$entrezID.",";
						}
						else
						{
							$entrezID_set{$sid."\t".$start."\t".$last}=",".$entrezID.",";
						}
					}
				}
				$entrezID_set{$sid."\t".$start."\t".$last}=~s/^,//g;$entrezID_set{$sid."\t".$start."\t".$last}=~s/,$//g;
				$annotation_hash{$tmp}="<Gene pmid='$TargetPmid' sid='$sid' start='$start' end='$last' tax_id='$tax_id' entrezID='".$entrezID_set{$sid."\t".$start."\t".$last}."' score='$weight'>$mention<\/Gene>";
			}
		}

		open ga,">>$ga_extraction";
		my @rank = sort {$annotation_rank_hash{$a} <=> $annotation_rank_hash{$b}} keys %annotation_rank_hash;
		foreach $tmp(@rank)
		{
			print ga $annotation_hash{$tmp}."\n";
		}
		foreach my $mention (keys %Domain_hash)
		{
			$Domain_hash{$mention}=~s/,$//g;
			if(exists $Domain_GM_hash{$mention})
			{
				$Domain_GM_hash{$mention}=~s/end='([0-9]+)'/end='$1' entrezID='$Domain_hash{$mention}' PSSMID='$mention2pssmids_hash{$mention}'/g;
				print ga $Domain_GM_hash{$mention};
			}
		}
		foreach my $mention (keys %Family_hash)
		{
			$Family_hash{$mention}=~s/,$//g;
			if(exists $Family_hash{$mention})
			{
				$Family_GM_hash{$mention}=~s/end='([0-9]+)'/end='$1' entrezID='$Family_hash{$mention}'/g;
				print ga $Family_GM_hash{$mention};
			}
		}
		close ga;
	}
	return 1;
}
return 1;