#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

BEGIN {
        push(@INC, '/home/weic4/BioC');
}
use BioC_full;

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
				$Annotation_hash{$pmid}{$start_position_of_passage}{($start-$gap_preSen)."\t".length($mention_tmp)."\t".$tax_mention."\tSpecies\tNCBI Taxonomy\t".$tax_id}=$start;
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
					$Annotation_hash{$pmid}{$start_position_of_passage}{($start-$gap_preSen)."\t".length($mention_tmp)."\t".$mention."\t".$type."\tNCBI Gene\t".$entrezID}=$start;
				}
			}
		}
		close ga;
	}
	
	my $STR_output="";
	my $pmid="";
	my $Anncount=1;
	
	#Read BioC
	my $input_collection = new BioC_full::Collection();
	my $input_xml = new BioC_full::Connector_libxml();
	$input_xml->start_read($inputfile, $input_collection);
	
	#Write BioC
	my $annotation_collection = new BioC_full::Collection();	# create annotation collection object
	my $annotate_xml = new BioC_full::Connector_libxml();	# create libxml Connector for output
	$annotation_collection->{source}=$input_collection->{source};
	$annotation_collection->{date}=$input_collection->{date};
	$annotation_collection->{key}=$input_collection->{key};
	$annotate_xml->start_write( $result, $annotation_collection );	# start to write contents of annotation collection to stdout
	
	my $document = new BioC_full::Document();
	while ( $input_xml->read_next($document) ) # documents
	{
		my $document_output = new BioC_full::Document();
		my $pmid = $document->{id};
		$document_output ->{id}=$pmid;
		
		my $count=0;my $annotation_id=0;
		my $passage_start=0;my $passage_text="";my $passage_last=0;my $passage_startlen=0;my $passage_startlen_next=0;
		for(my $i=0;$i<$document->{passages}->size();$i++) # passages
		{
			my $psg = new BioC_full::Passage();	# create Passage object
			
			if($document->{passages}->get($i)->{sentences}->size()>0)
			{
				for(my $j=0;$j<$document->{passages}->get($i)->{sentences}->size();$j++) # sentence
				{
					my $Sen = new BioC_full::Sentence();	# create Passage object
					$Sen->{offset} = $document->{passages}->get($i)->{sentences}->get($j)->{offset};	# copy offset element contents
					$Sen->{text} = $document->{passages}->get($i)->{sentences}->get($j)->{text};	# copy text
					
					my @AnnotationArr = sort {$Annotation_hash{$pmid}{$Sen->{offset}}{$a} <=> $Annotation_hash{$pmid}{$Sen->{offset}}{$b}} keys %{ $Annotation_hash{$pmid}{$Sen->{offset}} };
					my $sentence_tmp=$sentence;
					foreach my $AnnotationArr (@AnnotationArr)
					{
						if($AnnotationArr=~/^([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)/)
						{
							my $startA=$1;
							my $lengthA=$2;
							my $mentionA=$3;
							my $typeA=$4;
							my $dbA=$5;
							my $identifierA=$6;
							
							my $annotation = new BioC_full::Annotation(); #define annotation
							$annotation->{id} = "".$annotation_id;
							$annotation->{infons}->set("type", $typeA);
							$annotation->{infons}->set($dbA, $identifierA);
							$annotation->add_location( $startA, $lengthA );
							$annotation->{text} = $mentionA;
							$Sen->{annotations}->push( $annotation );	# can I use shift?
							$annotation_id++;
						}				
					}
					
					$psg->{sentences}->push( $Sen );
				}
			}
			else
			{			
				$psg->{infons}->set("type", $document->{passages}->get($i)->{infons}->get("type"));	# copy type information
				$psg->{offset} = $document->{passages}->get($i)->{offset};	# copy offset element contents
				$psg->{text} = $document->{passages}->get($i)->{text};	# copy text
					
				##====
				##Added existed annotation from input to io
				#for(my $j=0;$j<$document->{passages}->get($i)->{annotations}->size();$j++)	
				#{
				#	$psg->{annotations}->push( $document->{passages}->get($i)->{annotations}->get($j) );
				#	$count++;
				#}
				#for(my $j=0;$j<$document->{passages}->get($i)->{relations}->size();$j++)
				#{
				#	$psg->{relations}->push( $document->{passages}->get($i)->{relations}->get($j) );
				#}
				
				#====
				#insert entities into BioC pasaages
				my @AnnotationArr = sort {$Annotation_hash{$pmid}{$psg->{offset}}{$a} <=> $Annotation_hash{$pmid}{$psg->{offset}}{$b}} keys %{ $Annotation_hash{$pmid}{$psg->{offset}} };
				my $sentence_tmp=$sentence;
				foreach my $AnnotationArr (@AnnotationArr)
				{
					if($AnnotationArr=~/^([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)	([^\t]+)/)
					{
						my $startA=$1;
						my $lengthA=$2;
						my $mentionA=$3;
						my $typeA=$4;
						my $dbA=$5;
						my $identifierA=$6;
						
						my $annotation = new BioC_full::Annotation(); #define annotation
						$annotation->{id} = "".$annotation_id;
						$annotation->{infons}->set("type", $typeA);
						$annotation->{infons}->set($dbA, $identifierA);
						$annotation->add_location( $startA, $lengthA );
						$annotation->{text} = $mentionA;
						$psg->{annotations}->push( $annotation );	# can I use shift?
						$annotation_id++;
					}				
				}
			}
			$document_output->{passages}->push( $psg );
		}
		
		##====
		##Added existed annotation from input to io
		#for(my $j=0;$j<$document->{relations}->size();$j++)
		#{
		#	$document_output->{relations}->push( $document->{relations}->get($j) );
		#}
		
		$annotate_xml->write_next( $document_output );
	}
	
	$annotate_xml->end_write();
	
	return 1;
}

return 1;