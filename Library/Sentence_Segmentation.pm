#!perl
#===================================================
# Software: GNormPlus
# Date: 2014/12/25
#===================================================

package GNormPlus;

sub Sentence_Segmentation #XML or Free Text
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	my $table_extraction="tmp/".$filename.".table.xml";
	
	my $context="";
	open inputfile,"<".$input."/".$filename;
	while(<inputfile>)
	{
		my $tmp=$_;
		$context=$context.$tmp;
	}
	close inputfile;
	$context=~s/[\n\r\t]//g;
	
#=pod	
	my %ent_hash=();
	open ent,"<".$dictionary."/ent.rev.txt";
	while(<ent>)
	{
		my $tmp=lc($_);
		if($tmp=~/^(.+)	(.+)$/)
		{
			my $tmp=$1;
			my $ent=$2;
			$ent=~s/[\t\r\n]//g;
			$tmp=~s/[\t\r\n]//g;
			$ent_hash{$ent}=$tmp;
			
		}
	}
	foreach my $ent (keys %ent_hash)
	{
		$context =~s/$ent/ $ent_hash{$ent} /g;
	}
#=cut
		
	my $sid_count=1;
	open sentence,">".$sentence_extraction;
	
	#Extract PMID
	my $pmid="";
	if($context=~/<pub-id pub-id-type="pmid">([^<>]+)<\/pub-id>/)
	{
		$pmid=$1;
	}
	
	#Extract TITLE caption content
	if($context =~/<title\-group>[ ]*<article\-title>(.+?)<\/article\-title>/  || $context =~/<[Aa]rticleTitle.*?>(.+?)<\/[Aa]rticleTitle>/ || $context=~/<BookTitle.*?>(.+?)<\/BookTitle>/)
	{
		$title=$1;
		while($title=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$title=$pre.$tag.$post
		}
		$title =~s/^ //g;
		if($title!~/\.$/){$title=$title.".";}
		print sentence "<TEXT pmid='$pmid' sid='TITLE_$sid_count'>$title<\/TEXT>\n";
		$sid_count++;
	}
	
	#Extract ABSTRACT caption content
	if($context =~/<[Aa]bstract.*?>(.*?)<\/[Aa]bstract>/)
	{
		$abstract=$1;
		while($abstract=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$abstract=$pre.$tag.$post
		}
		$abstract =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
		$abstract =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
		@split_abstract=split("<SPLIT>",$abstract);
		foreach $sentence ( @split_abstract )
		{
			$sentence =~s/^ //g;
			if($sentence!~/\.$/){$sentence=$sentence.".";}
			print sentence "<TEXT pmid='$pmid' sid='ABSTRACT_$sid_count'>$sentence<\/TEXT>\n";
			$sid_count++;
		}
	}
	elsif($context =~/<[Aa]bstractText.*?>(.*?)<\/[Aa]bstractText>/ )
	{
		$abstract=$1;
		while($abstract=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$abstract=$pre.$tag.$post
		}
		$abstract =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
		$abstract =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
		@split_abstract=split("<SPLIT>",$abstract);
		foreach $sentence ( @split_abstract )
		{
			$sentence =~s/^ //g;
			if($sentence!~/\.$/){$sentence=$sentence.".";}
			print sentence "<TEXT pmid='$pmid' sid='ABSTRACT_$sid_count'>$sentence<\/TEXT>\n";
			$sid_count++;
		}
	}
	else
	{
		$abstract=$context;
		while($abstract=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$abstract=$pre.$tag.$post
		}
		$abstract =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
		$abstract =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
		@split_abstract=split("<SPLIT>",$abstract);
		foreach $sentence ( @split_abstract )
		{
			$sentence =~s/^ //g;
			if($sentence!~/\.$/){$sentence=$sentence.".";}
			print sentence "<TEXT pmid='$pmid' sid='ABSTRACT_$sid_count'>$sentence<\/TEXT>\n";
			$sid_count++;
		}
	}
	
	#Extract Author Summary caption content
	if($context =~/<abstract abstract.+?><title>Author Summary<\/title>(.*?)<\/abstract>/)
	{
		$Author_Summary=$1;
		while($Author_Summary=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$Author_Summary=$pre.$tag.$post
		}
		$Author_Summary =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
		$Author_Summary =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
		@split_Author_Summary=split("<SPLIT>",$Author_Summary);
		foreach $sentence ( @split_Author_Summary )
		{
			$sentence =~s/^ //g;
			if($sentence!~/\.$/){$sentence=$sentence.".";}
			print sentence "<TEXT pmid='$pmid' sid='Author_Summary_$sid_count'>$sentence<\/TEXT>\n";
			$sid_count++;
		}
	}
	
	#Extract abbreviations caption content
	if($context =~/<glossary><title>Abbreviations<\/title>(.+?)<\/glossary>/)
	{
		$Abbreviations=$1;
		#$Abbreviations =~s/<.+?>//g;
		#$Abbreviations =~s/[ ]+/ /g;
		print sentence "<TEXT sid='ABBREVIATIONS_$sid_count'>$Abbreviations<\/TEXT>\n";
		$sid_count++;
	}
	$context =~s/<glossary><title>Abbreviations<\/title>.+?<\/glossary>//g;
	
	#Extract figure caption content
	@figure=($context =~/<label>(Figure.+?<\/label><caption>.+?)<\/caption>/g);
	foreach $caption (@figure)
	{
		if($caption =~/(Figure.+)<\/label><caption>(.+)/)
		{
			$figure_num=$1;
			$caption=$2;
			while($caption=~/^(.*)(<.+?>)(.*)$/)
			{
				my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
				$caption=$pre.$tag.$post
			}
			$caption =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
			$caption =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
			@split_caption=split("<SPLIT>",$caption);
			foreach $sentence ( @split_caption )
			{
				$sentence =~s/^ //g;
				if($sentence!~/\.$/){$sentence=$sentence.".";}
				print sentence "<TEXT pmid='$pmid' sid='Figure_$sid_count'>$sentence<\/TEXT>\n";
				$sid_count++;
			}
		}
		
	}
	$context =~s/<label>Figure.+?<\/label><caption>.+?<\/caption>//g;
	
	#Extract table caption content
	@table=($context =~/<label>(Table.+?<\/label><caption>.+?)<\/caption>/g);
	foreach $caption (@table)
	{
		if($caption =~/(Table.+)<\/label><caption>(.+)/)
		{
			$table_num=$1;
			$caption=$2;
			while($caption=~/^(.*)(<.+?>)(.*)$/)
			{
				my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
				$caption=$pre.$tag.$post
			}
			$caption =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
			$caption =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
			@split_caption=split("<SPLIT>",$caption);
			foreach $sentence ( @split_caption )
			{
				$sentence =~s/^ //g;
				if($sentence!~/\.$/){$sentence=$sentence.".";}
				print sentence "<TEXT pmid='$pmid' sid='Table_$sid_count'>$sentence<\/TEXT>\n";
				$sid_count++;
			}
		}
		
	}
	$context =~s/<label>Table.+?<\/label><caption>.+?<\/caption>//g;
	
	#Extract body content
	if($context =~/<body>(.+?)<\/body>/)
	{
		$body=$1;
		while($body=~/^(.*)(<.+?>)(.*)$/)
		{
			my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
			$body=$pre.$tag.$post
		}
		$body =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
		$body =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
		@split_body=split("<SPLIT>",$body);
		foreach $sentence ( @split_body )
		{
			$sentence =~s/^ //g;
			if($sentence!~/\.$/){$sentence=$sentence.".";}
			print sentence "<TEXT pmid='$pmid' sid='BODY_$sid_count'>$sentence<\/TEXT>\n";
			$sid_count++;
		}
	}
	close sentence;

	#Extract table content
	if($context =~/<table.+?(>.+?<)\/table>/i)
	{
		open tablefile,">$table_extraction";
		@table=($context =~/<table.+?(>.+?<)\/table>/ig);
		foreach $table_info (@table)
		{
			@records=($table_info =~/>(.+?)</g);
			foreach $each_record (@records)
			{
				$each_record =~s/^\s+//g;
				$each_record =~s/\s+$//g;
				$each_record =~s/\&\#x[a-z0-9][a-z0-9][a-z0-9][a-z0-9][a-z0-9];//g;
				if($each_record !~/[\s\'\"]/ && $each_record =~/[0-9]+/  && $each_record =~/[A-Za-z]+/ )
				{
					if($each_record =~ /^(.+[^0-9])([0-9]+)\-([0-9]+)$/)#numeration
					{
						$tmp1=$1;
						$num1=$2;
						$num2=$3;
						if(($num2-$num1)<=10)
						{
							for($i=$num1;$i<=$num2;$i++)
							{
								print tablefile "<TABLE pmid='$pmid'>".$tmp1."".$i."<\/TABLE>\n"
							}
						}
					}
					print tablefile "<TABLE pmid='$pmid'>$each_record<\/TABLE>\n"
				}
			}
		}
		close tablefile;
	}
	
	return 1;	
}

sub Sentence_Segmentation_PubTator
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	
	my $context="";
	my $sid_count=1;
	open sentence,">".$sentence_extraction;
	open inputfile,"<".$input."/".$filename;
	while(<inputfile>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if ($tmp=~/^(.+)\|t\|(.*)$/)
		{
			my $pmid=$1;
			my $title=$2;
			$sid_count=1;
			print sentence "<TEXT pmid='$pmid' sid='TITLE_$sid_count'>$title<\/TEXT>\n";
		}
		elsif ($tmp=~/^(.+)\|a\|(.*)$/)
		{
			my $pmid=$1;
			my $abstract=$2;
			while($abstract=~/^(.*)(<.+?>)(.*)$/)
			{
				my $pre=$1;	my $tag=$2; my $post=$3; $tag=~s/./ /g;
				$abstract=$pre.$tag.$post
			}
			$abstract =~ s/\.([ ]+)([A-Z\(])/\.<SPLIT>$1$2/g;
			$abstract =~ s/(\([^\(\)]+?)\.<SPLIT>([^\(\)]+?\))/$1\.$2/g;
			@split_abstract=split("<SPLIT>",$abstract);
			foreach $sentence ( @split_abstract )
			{
				$sentence =~s/^ //g;
				if($sentence!~/\.$/){$sentence=$sentence.".";}
				$sid_count++;
				print sentence "<TEXT pmid='$pmid' sid='ABSTRACT_$sid_count'>$sentence<\/TEXT>\n";
			}
		}
	}
	close inputfile;
	close sentence;
	return 1;
}

sub Sentence_Segmentation_BioC
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	
	my $context="";
	open inputfile,"<".$input."/".$filename;
	while(<inputfile>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		$context=$context.$tmp." ";
	}
	close inputfile;
	
	
	open sentence,">".$sentence_extraction;
	my @documents=split("</document>",$context);
	foreach my $document (@documents)
	{
		my @tmp=($document=~/<id>(.+?)<\/id>/);my $pmid=$tmp[0];
		my @passages=split("</passage>",$document);
		foreach my $passage (@passages)
		{
			my @tmp=($passage=~/<offset>(.+?)<\/offset>/);my $offset=$tmp[0];
			my @tmp=($passage=~/<text>(.+?)<\/text>/);my $text=$tmp[0];
			
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
	close sentence;
	return 1;
}

#### Abbreviation resolution: AB3P ####
sub Ab3PExt
{
	my ($input)=@_[0];
	my ($filename)=@_[1];
	my ($dictionary)=@_[2];
	
	my $sentence_extraction="tmp/".$filename.".sentence.xml";
	
	my $pmid="";
	my $context="";
	open input,"<".$sentence_extraction;
	open output,">tmp/".$filename.".Abb";
	while(<input>)
	{
		my $tmp=$_;
		$tmp=~s/[\n\r]//g;
		if($tmp=~/^<TEXT pmid='([^<>]+)' sid='[^<>]+'>(.+)<\/TEXT>$/)
		{
			my $t1=$1;
			my $t2=$2;
			$t2=~s/<[^<>]+>/ /g;
			if($pmid ne "" && $pmid ne $t1)
			{
				print output $pmid."\n".$context."\n\n";
				$context="";
				$pmid=$t1;
			}
			else
			{
				$pmid=$t1;
				$context=$context.$t2." ";
			}
		}
	}
	print output $pmid."\n".$context."\n\n";
	close output;
	close input;
	
	if($^O=~/Win/)
	{
		print "Worming: Operation system is ".$^O.". Ab3P is not workable.              \n";
	}
	else
	{
		open execabb,"|./Ab3P tmp/".$filename.".Abb tmp/".$filename.".Ab3P";
		close execabb;
	}
	
	return 1;
}

return 1;