#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

sub Output_Result_PubTator
{
	my ($input)=@_[0];
	my ($output)=@_[1];
	my ($filename)=@_[2];
	my ($dictionary)=@_[3];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $tax_extraction="tmp/".$filename.".tax.xml";
	my $ga_extraction=$output."/".$filename.".ga.xml";
	my $result=$output."/".$filename.".PubTator";
	
	@Pmidlist=();
	my $pmid="";
	open sentence,"<".$sentence_extraction;
	while(<sentence>)
	{
		my $tmp=$_;
		if($tmp=~/^<TEXT pmid='(.+)' sid='.+'>.+<\/TEXT>/)
		{
			if($pmid ne $1)
			{
				push(@Pmidlist,$1);
				$pmid=$1;
			}
		}
	}
	close sentence;
	
	open output,">".$result;
	foreach my $pmid(@Pmidlist)
	{
		my $startpos=0;
		my $title="";
		my $abstract="";
		my %sentence2start=();
		open sentence,"<".$sentence_extraction;
		while(<sentence>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<TEXT pmid='$pmid' sid='((TITLE|ABSTRACT)\_([0-9]+))'>(.+)<\/TEXT>/i)
			{
				my $tiabssid=$1;
				my $tiabs=$2;
				my $sid=$3;
				my $sentence=$4;

				if($tiabs eq "TITLE")
				{
					$title=$sentence;
				}
				else #ABSTRACT
				{
					if($abstract eq "")
					{
						$abstract=$sentence;
					}
					else
					{
						$abstract=$abstract." ".$sentence;
					}
				}
				$sentence2start{$tiabssid}=$startpos;
				$startpos=$startpos+length($sentence)+1;
			}
		}
		close sentence;
		#$abstract=~s/^ //g;
		
		my %Annotation_hash=();
		#tax
		open tax,"<".$tax_extraction;
		while(<tax>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<Tax pmid='$pmid' sid='((TITLE|ABSTRACT)\_([0-9]+))' start='(.+)' end='(.+)' tax_id='(.+)'>(.+)<\/Tax>/i)
			{
				my $tiabssid=$1;
				my $tiabs=$2;
				my $sid=$3;
				my $start=$4+$sentence2start{$tiabssid};
				my $last=$5+$sentence2start{$tiabssid};
				my $tax_id=$6;
				my $tax_mention=$7;
				$Annotation_hash{$pmid}{$pmid."\t".$start."\t".$last."\t".$tax_mention."\tSpecies\t".$tax_id}=$start;
			}
		}
		close tax;
		
		#ga
		open ga,"<".$ga_extraction;
		while(<ga>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if(	$tmp=~/<([^ ]+) pmid='$pmid' sid='((TITLE|ABSTRACT)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' tax_id='[^ ]+' entrezID='([^ ]+)' score='([^ ]+)'>(.+)<\/.+>/i ||
				$tmp=~/<(FamilyName) pmid='$pmid' sid='((TITLE|ABSTRACT)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' entrezID='([^ ]+)' score='([^ ]+)'>(.+)<\/FamilyName>/i ||
				$tmp=~/<(DomainMotif) pmid='$pmid' sid='((TITLE|ABSTRACT)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' entrezID='([^ ]+)' PSSMID='([^ ]*)' score='[^ ]+'>(.+)<\/DomainMotif>/i
			)
			{
				my $type=$1;
				my $tiabssid=$2;
				my $tiabs=$3;
				my $sid=$4;
				my $start=$5+$sentence2start{$tiabssid};
				my $last=$6+$sentence2start{$tiabssid};
				my $entrezID=$7;
				my $PSSMID=$8;
				my $mention=$9;
				if($type eq "DomainMotif")
				{
					$entrezID="GeneID:".$entrezID;
					if($PSSMID ne "")
					{
						$entrezID=$entrezID."/PSSMID:".$PSSMID;
					}
				}
				$Annotation_hash{$pmid}{$pmid."\t".$start."\t".$last."\t".$mention."\t".$type."\t".$entrezID}=$start;
			}
		}
		close ga;
		
		print output $pmid."|t|".$title."\n";
		print output $pmid."|a|".$abstract."\n";
		my @Anno = sort {$Annotation_hash{$pmid}{$a} <=> $Annotation_hash{$pmid}{$b}} keys %{$Annotation_hash{$pmid}};
		foreach my $ann (@Anno) 
		{
			print output $ann."\n";
		}
		print output "\n";
	}
	close output;
	return 1;
}

sub Output_Result_BioC
{
	my ($input)=@_[0];
	my ($output)=@_[1];
	my ($filename)=@_[2];
	my ($dictionary)=@_[3];
	
	my $inputfile=$input."/".$filename;
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $tax_extraction="tmp/".$filename.".tax.xml";
	my $ga_extraction=$output."/".$filename.".ga.xml";
	my $result=$output."/".$filename.".BioC.XML";
	
	@Pmidlist=();
	my $pmid="";
	open sentence,"<".$sentence_extraction;
	while(<sentence>)
	{
		my $tmp=$_;
		if($tmp=~/^<TEXT pmid='(.+)' sid='.+'>.+<\/TEXT>/)
		{
			if($pmid ne $1)
			{
				push(@Pmidlist,$1);
				$pmid=$1;
			}
		}
	}
	close sentence;
	
	my %Annotation_hash=();
	foreach my $pmid(@Pmidlist)
	{
		my $start=0;
		my %sid2setnence=();
		my %sentence2start=();
		open sentence,"<".$sentence_extraction;
		while(<sentence>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<TEXT pmid='$pmid' sid='(([0-9]+)\_([0-9]+))'>(.+)<\/TEXT>/i)
			{
				my $Global_sid=$1;
				my $start_position_of_passage=$2;
				my $sid=$3;
				my $sentence=$4;
				if($sid==1)
				{
					$start=0;
				}
				$sentence2start{$Global_sid}=$start_position_of_passage+$start;
				
				my $preSen=$sentence;
				my $length_preSen=length($preSen);
				$preSen=~s/\&[\#A-Za-z0-9]+\;/\./g;
				my $gap_preSen=$length_preSen-length($preSen);
				
				$start=$start+length($sentence)+1-$gap_preSen;
				$sid2setnence{$Global_sid}=$sentence;
			}
		}
		close sentence;
		#$abstract=~s/^ //g;
		
		#tax
		open tax,"<".$tax_extraction;
		while(<tax>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if($tmp=~/<Tax pmid='$pmid' sid='(([0-9]+)\_([0-9]+))' start='(.+)' end='(.+)' tax_id='(.+)'>(.+)<\/Tax>/i)
			{
				my $Global_sid=$1;
				my $start_position_of_passage=$2;
				my $sid=$3;
				my $start=$4+$sentence2start{$Global_sid};
				my $tax_id=$6;
				my $tax_mention=$7;
				
				my $preSen=substr($sid2setnence{$Global_sid},0,($start-$sentence2start{$Global_sid}));
				my $length_preSen=length($preSen);
				$preSen=~s/\&[\#A-Za-z0-9]+\;/\./g;
				my $gap_preSen=$length_preSen-length($preSen);
								
				$mention_tmp=$tax_mention;
				$mention_tmp=~s/\&[\#A-Za-z0-9]+\;/\./g;
				$Annotation_hash{$pmid}{$start_position_of_passage}{"<infon key='type'>Species</infon><infon key='NCBI Taxonomy'>$tax_id</infon><location offset='".($start-$gap_preSen)."' length='".length($mention_tmp)."' />\n<text>$tax_mention</text></annotation>"}=$start;
			}
		}
		close tax;
		
		#ga
		open ga,"<".$ga_extraction;
		while(<ga>)
		{
			my $tmp=$_;
			$tmp=~s/[\n\r]//g;
			if(	$tmp=~/<([^ ]+) pmid='$pmid' sid='((TITLE|ABSTRACT|[0-9]+)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' tax_id='[^ ]+' entrezID='([^ ]+)' score='([^ ]+)'>(.+)<\/.+>/i ||
				$tmp=~/<(FamilyName) pmid='$pmid' sid='((TITLE|ABSTRACT|[0-9]+)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' entrezID='([^ ]+)' score='([^ ]+)'>(.+)<\/FamilyName>/i ||
				$tmp=~/<(DomainMotif) pmid='$pmid' sid='((TITLE|ABSTRACT|[0-9]+)\_([0-9]+))' start='([^ ]+)' end='([^ ]+)' entrezID='([^ ]+)' PSSMID='([^ ]*)' score='[^ ]+'>(.+)<\/DomainMotif>/i
			)
			{
				my $type=$1;
				my $Global_sid=$2;
				my $start_position_of_passage=$3;
				my $sid=$4;
				my $start=$5+$sentence2start{$Global_sid};
				my $entrezID=$7;
				my $PSSMID=$8;
				my $mention=$9;
				if($type eq "DomainMotif"){$entrezID="GeneID:".$entrezID."/PSSMID:".$PSSMID;}
				
				my $preSen=substr($sid2setnence{$Global_sid},0,($start-$sentence2start{$Global_sid}));
				my $length_preSen=length($preSen);
				$preSen=~s/\&[\#A-Za-z0-9]+\;/\./g;
				my $gap_preSen=$length_preSen-length($preSen);
				
				$mention=~s/\.$//g;
				if(($mention=~/\&/ && $mention!~/\;/) || length($mention)>=200){}else
				{
					$mention_tmp=$mention;
					$mention_tmp=~s/\&[\#A-Za-z0-9]+\;/\./g;
					$Annotation_hash{$pmid}{$start_position_of_passage}{"<infon key='type'>$type</infon><infon key='NCBI Gene'>$entrezID</infon><location offset='".($start-$gap_preSen)."' length='".length($mention_tmp)."' />\n<text>$mention</text></annotation>"}=$start;
				}
			}
		}
		close ga;
	}
	
	my $BioC_text="";;
	open inputfile,"<".$inputfile;
	while(<inputfile>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$BioC_text=$BioC_text.$tmp." ";
	}
	close inputfile;

	my $STR_output="";
	my $pmid="";
	my $Anncount=1;
	#$BioC_text=~s/    [ ]+/ /g;
	if($BioC_text=~/^(.+?)(<document>.+)$/){$STR_output=$STR_output.$1;$BioC_text=$2;} #document-Start
	while($BioC_text=~/^(.*?)<document>(.+?)<\/document>(.*)$/)
	{
		my $annotation_id=0;
		my $pre_doc=$1;
		my $document=$2;
		$BioC_text=$3;
		$STR_output=$STR_output.$pre_doc;	#document-Start
		$STR_output=$STR_output."<document>";
		if($document=~/^[\W\-\_]*<id>(.+?)<\/id>/){$pmid=$1;}
		
		my $tmp_passage=$document;
		while($tmp_passage=~/^(.*?)<passage>(.+?)<\/passage>(.*)$/)
		{
			my $pre_passage=$1;
			my $passage=$2;
			$tmp_passage=$3;
			my $sentence="";
			my $offset=0;
			$STR_output=$STR_output.$pre_passage;	#passage-Start
			$STR_output=$STR_output."<passage>";
			if($passage=~/^(.*?)<infon key="type">(.+?)<\/infon>(.*?)<offset>([0-9]+)<\/offset>(.*?)<text>(.*?)<\/text>(.*?)$/)
			{
				my $pre=$1;
				my $type=$2;
				my $XX1=$3;
				$offset=$4;
				my $XX2=$5;
				$sentence=$6;
				my $post=$7;
				
				$STR_output=$STR_output.$pre;
				$STR_output=$STR_output."<infon key=\"type\">$type<\/infon>";
				$STR_output=$STR_output.$XX1;
				$STR_output=$STR_output."<offset>$offset<\/offset>";
				$STR_output=$STR_output.$XX2;
				$STR_output=$STR_output."<text>$sentence<\/text>";
				$STR_output=$STR_output.$post;
			}
			else
			{
				$STR_output=$STR_output.$passage;
			}
			
			#Annotation insert
			my @AnnotationArr = sort {$Annotation_hash{$pmid}{$offset}{$a} <=> $Annotation_hash{$pmid}{$offset}{$b}} keys %{ $Annotation_hash{$pmid}{$offset} };
			my $sentence_tmp=$sentence;
			foreach my $AnnotationArr (@AnnotationArr)
			{
				$STR_output=$STR_output."<annotation id='".$Anncount++."'>\n".$AnnotationArr;
			}
			$STR_output=$STR_output."<\/passage>";
		}
		$STR_output=$STR_output.$tmp_passage;	#passage-End
		$STR_output=$STR_output."<\/document>";
	}
	$STR_output=$STR_output.$BioC_text;	#document-End
	
	$STR_output=~s/(<(collection|document|passage|annotation)>)/$1\n/g;
	$STR_output=~s/(<\/([A-Za-z0-9]+)>)/$1\n/g;
	open output,">".$result;	
	print output $STR_output."\n";
	close output;
	return 1;
}

return 1;