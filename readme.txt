[Directory]
A. Introduction of folders
B. Installation
C. Instruction
D. FULL Usage
E. Retrieve articles via E-utilities (For PubMed abstracts & PMC fulltext articles)
F. reference

#======================================================================================#

A. [Introduction of folders]
	
	Corpus:
		[NCBIGeneCorpus]
	Library & Module:
		[Library]
		[CRFmodel]
	Dictionary:
		[dictionary]
	Input data folder:
		[input]
	Output data folder:
		[output]
	TMP folder:
		[tmp]
	
B. [Installation] 

	*Users don't need to install CRF++ module if their operation system is windows.
	
	*Installation
		*Windows Environment
			Users need to install Perl. User can download ActivePerl via this link: http://www.activestate.com/activeperl/downloads
		
		*Linux Environment
			Users need to execute "Installation.sh" to install crf_test first.
	
			$ sh Installation.sh
					
C. [Instruction]

	GNormPlus is an updated version of GenNorm. GNormPlus already integrates SR4GN for species identification, SimConcept for handling composite entity names, and Ab3P for abbreviation resolution. Thus, the installation of these separate tools are not necessary. 
	
	$ perl GNormPlus.pl -i [Input] -o [Output] -s [Setup]
	
	Example: perl GNormPlus.pl -i input -o output -s setup.txt
	
	Input: User can provide the input data folder route ("input" or "C:\input"). 
	Output: User can provide the output data folder route ("output" or "C:\output"). 
	Setup: User can choose focus species and other functions that they would like to enable/disable. Default is "setup.txt". ("setup.txt" or "C:\setup.txt"). 


D. [FULL Usage] 

	GNormPlus is developed to find gene/proteins, its family names and domain motifs.	

	INPUT:
	
		Input file folder. Each input file should follow the PubTator format below or the BioC format (http://bioc.sourceforge.net/). 

	RUN: 

		perl GNormPlus.pl -i input -o output -s setup.txt

	OUTPUT:

		Output file folder. Each output file should follow the PubTator format below or BioC format(http://bioc.sourceforge.net/). 
				
	Input File Format (PubTator format):
	
		<PMID1>|t|<TITLE>
		<PMID1>|a|<ABSTRACT>
		<PMID1>	<START>	<LAST>	<MENTION>	<TYPE>	<Normalized Form>
		
		<PMID2>|t|<TITLE>
		<PMID2>|a|<ABSTRACT>
		<PMID2>	<START>	<LAST>	<MENTION>	<TYPE>	<Normalized Form>
		...
		
		Example: 
		10022127|t|TIF1gamma, a novel member of the transcriptional intermediary factor 1 family.
		10022127|a|We report the cloning and characterization of a novel member of the Transcriptional Intermediary Factor 1 (TIF1) gene family, human TIF1gamma. Similar to TIF1alpha and TIF1beta, the structure of TIF1beta is characterized by multiple domains: RING finger, B boxes, Coiled coil, PHD/TTC, and bromodomain. Although structurally related to TIF1alpha and TIF1beta, TIF1gamma presents several functional differences. In contrast to TIF1alpha, but like TIF1beta, TIF1 does not interact with nuclear receptors in yeast two-hybrid or GST pull-down assays and does not interfere with retinoic acid response in transfected mammalian cells. Whereas TIF1alpha and TIF1beta were previously found to interact with the KRAB silencing domain of KOX1 and with the HP1alpha, MODI (HP1beta) and MOD2 (HP1gamma) heterochromatinic proteins, suggesting that they may participate in a complex involved in heterochromatin-induced gene repression, TIF1gamma does not interact with either the KRAB domain of KOX1 or the HP1 proteins. Nevertheless, TIF1gamma, like TIF1alpha and TIF1beta, exhibits a strong silencing activity when tethered to a promoter. Since deletion of a novel motif unique to the three TIF1 proteins, called TIF1 signature sequence (TSS), abrogates transcriptional repression by TIF1gamma, this motif likely participates in TIF1 dependent repression.

		10072587|t|Cloning of a novel gene (ING1L) homologous to ING1, a candidate tumor suppressor.
		10072587|a|The ING1 gene encodes p33(ING1), a putative tumor suppressor for neuroblastomas and breast cancers, which has been shown to cooperate with p53 in controlling cell proliferation. We have isolated a novel human gene, ING1L, that potentially encodes a PHD-type zinc-finger protein highly homologous to p33(ING1). Fluorescence in situ hybridization and radiation-hybrid analyses assigned ING1L to human chromosome 4. Both ING1 and ING1L are expressed in a variety of human tissues, but we found ING1L expression to be significantly more pronounced in tumors from several colon-cancer patients than in normal colon tissues excised at the same surgical sites. Although the significance of this observation with respect to carcinogenesis remains to be established, the data suggest that ING1L might be involved in colon cancers through interference with signal(s) transmitted through p53 and p33(ING1).

		10072769|t|Selection of cDNAs encoding putative type II membrane proteins on the cell surface from a human full-length cDNA bank.
		10072769|a|We have developed a simple method to test whether a hydrophobic segment near the N-terminus of a protein functions as a type II signal anchor (SA) in which the N-terminus faces the cytoplasm. A cDNA fragment containing the putative SA sequence of a target clone was fused in-frame to the 5' end of a cDNA fragment encoding the protease domain of urokinase-type plasminogen activator (u-PA). The resulting fused gene was expressed in COS7 cells. Fibrinolytic activity on the cell surface was measured by placing a fibrin sheet in contact with the transfected COS7 cells after removing the medium. When the cDNA fragment encoded a SA, the fibrin sheet was lysed by the u-PA expressed on the cell surface. The fibrinolytic activity was not detected in the culture medium, suggesting that the u-PA remains on the cell surface anchored via the SA in the membrane without being cleaved by signal peptidase. This fibrin sheet method was successfully applied to select five novel cDNA clones encoding putative type II membrane proteins from a human full-length cDNA bank.
		
	Output File Format:

		<PMID1>|t|<TITLE>
		<PMID1>|a|<ABSTRACT>
		<PMID1>	<START>	<LAST>	<MENTION>	<TYPE>	<Normalized Form>
		
		<PMID2>|t|<TITLE>
		<PMID2>|a|<ABSTRACT>
		<PMID2>	<START>	<LAST>	<MENTION>	<TYPE>	<Normalized Form>
		...
		
		Example: 
		10022127|t|TIF1gamma, a novel member of the transcriptional intermediary factor 1 family.
		10022127|a|We report the cloning and characterization of a novel member of the Transcriptional Intermediary Factor 1 (TIF1) gene family, human TIF1gamma. Similar to TIF1alpha and TIF1beta, the structure of TIF1beta is characterized by multiple domains: RING finger, B boxes, Coiled coil, PHD/TTC, and bromodomain. Although structurally related to TIF1alpha and TIF1beta, TIF1gamma presents several functional differences. In contrast to TIF1alpha, but like TIF1beta, TIF1 does not interact with nuclear receptors in yeast two-hybrid or GST pull-down assays and does not interfere with retinoic acid response in transfected mammalian cells. Whereas TIF1alpha and TIF1beta were previously found to interact with the KRAB silencing domain of KOX1 and with the HP1alpha, MODI (HP1beta) and MOD2 (HP1gamma) heterochromatinic proteins, suggesting that they may participate in a complex involved in heterochromatin-induced gene repression, TIF1gamma does not interact with either the KRAB domain of KOX1 or the HP1 proteins. Nevertheless, TIF1gamma, like TIF1alpha and TIF1beta, exhibits a strong silencing activity when tethered to a promoter. Since deletion of a novel motif unique to the three TIF1 proteins, called TIF1 signature sequence (TSS), abrogates transcriptional repression by TIF1gamma, this motif likely participates in TIF1 dependent repression.
		10022127	0	9	TIF1gamma	Gene	51592
		10022127	33	70	transcriptional intermediary factor 1	FamilyName	51592,8805
		10022127	147	184	Transcriptional Intermediary Factor 1	FamilyName	51592,8805
		10022127	186	190	TIF1	FamilyName	51592,8805
		10022127	211	220	TIF1gamma	Gene	51592
		10022127	233	242	TIF1alpha	Gene	8805
		...

E. [Retrieve PubMed or PMC articles for pre-processing]
	
	$ perl PreProcessing.pl -t [type] -i [input] -o [output]
	
	Example: perl PreProcessing.pl -t PMID -i 22016685 -o 22016685.txt
	
	type: User can choose to one of the input type : PMID | PMCID | PMIDlist | PMCIDlist (PMIDlist and PMCIDlist are files which store the pmid or pmcid list.)
	Input: User can provide the input parameter, such as 22016685(for PMID type), PMIDlist.txt(for PMIDlist type). 
	Input: User can provide the output file name (input.txt). 
	
F. [Reference]

	Chih-Hsuan Wei, Hung-Yu Kao, Zhiyong Lu (2015) "GNormPlus: An Integrative Approach for Tagging Genes, Gene Families and Protein Domains”, BioMed Res Int., in press.
	Chih-Hsuan Wei, Robert Leaman, Zhiyong Lu (2015) "SimConcept: A Hybrid Approach for Simplifying Composite Named Entities in Biomedical Text", IEEE Journal of Biomedical and Health Informatics, DOI:10.1109/JBHI.2015.2422651
	Chih-Hsuan Wei, Hung-Yu Kao, Zhiyong Lu (2012) "SR4GN: a species recognition software tool for gene normalization", PLoS ONE,7(6):e38460 