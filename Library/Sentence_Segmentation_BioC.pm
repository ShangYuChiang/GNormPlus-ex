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

sub Sentence_Segmentation_BioC
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	
	#read BioC
	my $input_collection = new BioC_full::Collection();
	my $input_xml = new BioC_full::Connector_libxml();
	$input_xml->start_read($input."/".$filename, $input_collection);
	
	open sentence,">$sentence_extraction";
	#output article for Abbreviation resolution
	my $document = new BioC_full::Document();
	while ( $input_xml->read_next($document) ) # documents
	{
		my $pmid = $document->{id};
		my $article="";
		for(my $i=0;$i<$document->{passages}->size();$i++) # passages
		{
			if($document->{passages}->get($i)->{sentences}->size()>0)
			{
				for(my $j=0;$j<$document->{passages}->get($i)->{sentences}->size();$j++) # sentence
				{
					my $type = "sentence"; #infon
					my $offset = $document->{passages}->get($i)->{sentences}->get($j)->{offset}; #offset
					my $text = $document->{passages}->get($i)->{sentences}->get($j)->{text}; #text
					$text=~s/[\n\r]//g;
					
					while($text=~/^(.*)(<.+?>)(.*)$/)
					{
						my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
						$text=$pre.$tag.$post
					}
					$text =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
					$text =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
					@split_text=split("<SPLIT>",$text);
					my $sid_count=1;
					foreach $sentence ( @split_text )
					{
						if($sid_count>1){$sentence =~s/^ //g;}
						if($sentence!~/\.$/){$sentence=$sentence.".";}
						#if($sid_count>1){$sentence=~s/    [ ]+/ /g;}
						print sentence "<TEXT pmid='".$pmid."' sid='".$offset."_".$sid_count."'>$sentence<\/TEXT>\n";
						$sid_count++;
					}
				}
			}
			else
			{
				my $type = $document->{passages}->get($i)->{infons}->get("type"); #infon
				my $offset = $document->{passages}->get($i)->{offset}; #offset
				my $text = $document->{passages}->get($i)->{text}; #text
				$text=~s/[\n\r]//g;
				
				while($text=~/^(.*)(<.+?>)(.*)$/)
				{
					my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
					$text=$pre.$tag.$post
				}
				$text =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
				$text =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
				@split_text=split("<SPLIT>",$text);
				my $sid_count=1;
				foreach $sentence ( @split_text )
				{
					if($sid_count>1){$sentence =~s/^ //g;}
					if($sentence!~/\.$/){$sentence=$sentence.".";}
					#if($sid_count>1){$sentence=~s/    [ ]+/ /g;}
					print sentence "<TEXT pmid='".$pmid."' sid='".$offset."_".$sid_count."'>$sentence<\/TEXT>\n";
					$sid_count++;
				}
			}
		}
	}
	close sentence;
	
	return 1;
}

return 1;