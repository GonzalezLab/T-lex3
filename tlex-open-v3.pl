#!/usr/bin/perl

#######################################################################################
# To check if an error message returns:                                               #
# ------------------------------------                                                #
# -check the full path of the data                                                    #
# -check the concordance between the TE list, TE annotation, reference sequence files #
# -check the format of the input illumina data                                        #
#######################################################################################


use warnings;
use strict;
use Getopt::Long;
use Time::Local;
use List::Util qw[min max];
use Data::Dumper;
use File::Basename;

&startup;

sub startup { 

    my $strains;
    my $refgenome; 
    my $TE_list; 
    my $newTElist=0;
    my $TE_map;
    my $maxReadLength=100;
    my $output;
    my $species="drosophila";
    my $noFilterTE;
    my $minflankcov=0.5;
    my $bwaonly; 
    my $junction=1000;
    my $buffer=60;
    my $shrimponly;
    my $flank=125;
    my $var=20;
    my $lima=10;
    my $limp=15;
    my $id=95;
    my $combine;
    my $combine_all;
    my $combine_strains;
    my $freqestim;
    my $align;
    my $noclean;
    my $help;
    my $binref; 
    my $binreads; 
    my $pooleddata; 
    my $tsd;
    my $reads=0; 
    my $minqual=30;
    my $processes=0;
    my $PE="no"; # require the reads of the same pair to be in separated files
    my $minreads=3;
    my $maxreads=90;
    my $minpop=1;

    my $startdirectory=`pwd`; 
    chomp $startdirectory;
    
    GetOptions ('R=s' => \$strains,
		'G=s' => \$refgenome,	
		'T=s' => \$TE_list,
		'M=s' => \$TE_map,
		'A=i' => \$maxReadLength,
		'O=s' => \$output,
		's=s' => \$species,
		'noFilterTE' => \$noFilterTE,
		'd=f' => \$minflankcov,		
		'q' => \$bwaonly,
		'j=i' => \$junction,
		'b=i' => \$buffer,
		'minQ=i' => \$minqual,
		'p' => \$shrimponly,
		'f=i' => \$flank, 
		'v=i' => \$var,	
		'lima=i' => \$lima,
		'limp=i' => \$limp,
		'id=i' => \$id,
		'combRes' => \$combine,
		'combAll' => \$combine_all,
		'combData' => \$combine_strains,
		'freq' => \$freqestim, 
		'align' => \$align,
		'tsd' => \$tsd, 
		'noclean' => \$noclean,
		'h|help' => \$help,                                                                                                                                                                                                                                                                                                                                                              
                'binref' => \$binref,
		'binreads' => \$binreads,
		'pooled' => \$pooleddata,
		'pairends=s' => \$PE,
		'processes=i' => \$processes,
        'minR=i' => \$minreads, # V2.5 for setting minimum number of reads for calculating frequencies in pools
        'maxR=i' => \$maxreads, # V2.5 for setting maximum number of reads for calculating frequencies in pools
        'minP=i' => \$minpop); # V2.5 for setting minimum number of individuals for calculating frequencies in individual strains

   
   
    if ($help) {
	&help();
	die "\n";
    }

    my $outputdir;
    if($output) {
		if (index($output, "_") != -1) {
			print "WARNING: strain name cannot contain character '_'!!!!";
			exit;
		} else {
			my $tmp="tlex_$output";
			$outputdir=$tmp;
		}
    }
    else{
	$outputdir="tlex_output"; 
    }

    if($combine_strains){
	print "\n\n **************************************************************************************************************************************\n";
	print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\tCombine the presence/absence results from different strain(s) *\n\t\t\t\t\t\t\n";
	print "\t\t\t * only the parameters for the combination of the presence and absence result may be necessary *\n";
	&combine_results("$startdirectory");
	die "\n";
    }
    
    if($combine_all){
	print "\n\n **************************************************************************************************************************************\n";
	print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\tCombine the frequency estimates, the analysis of the TE flanking regions and the TSD detection *\n\t\t\t\t\t\t\n";
	if ($pooleddata){
	    &CombineAll("$outputdir","$pooleddata");
	}
	else{
	    &CombineAll("$outputdir","");  
	}
	die "\n";
    }

    my $new_TE_list;
    if($combine){
	print "\n\n **************************************************************************************************************************************\n";
	print "\n\t\t\t\t\t\t\t  * T-lex release 2 \n\t\t\t\tCombine the presence/absence results from one strain *\n\t\t\t\t\t\t\n";
	print "* Specify the project name if necessary\n";
	unless ($TE_list) {
	    die "* Must specify the TE list after the TE filtering step\n";
	}
	if ($noFilterTE) {
	    $new_TE_list=$TE_list;
	}
	else{
	    $new_TE_list="$TE_list\_filtered";
	}
	print "* list of TE used : $new_TE_list\n";
	print "$outputdir - $new_TE_list - $flank\n";
	&FinalResults("$outputdir","$new_TE_list", "$TE_map", "$flank");
	die "\n";
    }
       
    if ($freqestim){

    print "$startdirectory\/$outputdir\n";
	
	if ($pooleddata) {
		if (-e "$outputdir\/Tresults"){
			print "\n\n **************************************************************************************************************************************\n";
			print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\tReturn the frequency estimates of given sequence(s) in multiple strain(s) *\n\t\t\t\t\t\t\n";
			&FreqEstimate("$startdirectory\/$outputdir","$pooleddata","$maxreads","$minreads","$minpop");
			die "\n";
		}
		else{
			die "* Must detect the presence/absence of given sequences before! So run T-lex first or give the good project name\n";
		}
	} else {
		print "$startdirectory\/Tfreqs_output\/Tresults";
		if (-e "$startdirectory\/Tfreqs_output\/Tresults"){
			print "\n\n **************************************************************************************************************************************\n";
			print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\tReturn the frequency estimates of given sequence(s) in multiple strain(s) *\n\t\t\t\t\t\t\n";
			&FreqEstimate("$startdirectory\/Tfreqs_output","$pooleddata","$maxreads","$minreads","$minpop");
			die "\n";
		}
		else{
			die "* Must detect the presence/absence of given sequences before! So run T-lex first or give the good project name\n";
		}
	}
	}
    if ($align){
	unless ($TE_list && $strains){
	    die "* Must specify the list of TEs, the sequencing data directory and the project name!\n";
	}
	else{
	    if($tsd){
		print "\n\n **************************************************************************************************************************************\n";
		print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\t\t\t TSD detection (requires the alignments) *\n\t\t\t\t\t\t\n";
		&MultiAlign("$outputdir","$TE_list","$strains","$flank","0","1");
	    }
	    else{
		print "\n\n **************************************************************************************************************************************\n";
		print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\t\t\t Multiple Alignment report only *\n\t\t\t\t\t\t\n";
	       
		unless ($shrimponly) {
		    &MultiAlign("$outputdir","$TE_list","$strains","$flank","1","0");
		}
		unless ($bwaonly) {
		    &MultiAlign("$outputdir","$TE_list","$strains","$flank","0","0");
		}
		else{
		    &MultiAlign("$outputdir","$TE_list","$strains","$flank","2","0");
		}
	    }
	    die "\n";
	}
    }
 
    unless ($TE_list || $TE_map || $refgenome || $strains) {
	&help();
    }
  
    unless ($TE_list) {
	die "* the file of the TE list is missing !\n";
    }
    unless ($TE_map) {
	die "* the TE annotation file is missing!\n";
    }
    unless ($refgenome) {
	die "* the file with the reference genomic sequences is missing!\n";
    }
    unless ($strains) {
	die "* the directory of the sequencing data is missing!\n";
    }

    if ($shrimponly && $bwaonly) {
	die "* Cannot specify [shrimponly] and [bwaonly] at the same time!\n";
    }
    unless ($processes) {
	$processes = 1;
    }
	
    my $original_TE_list=$TE_list;
    system("mkdir $outputdir");
    open (PARAM, ">$outputdir\/Tparam");
    
    print "\n\n **************************************************************************************************************************************\n";
    if ($bwaonly) {
	print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t  Detect the presence of given sequence(s) in strain(s) using BWA program *\n";
	print PARAM "T-lex release 1 Detect the presence of given sequence(s) in strain(s) using MAQ program\n";
    }
    elsif ($shrimponly) {
	print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t  Detection of the absence of given sequence(s) in strain(s) using SHRIMP program *\n"; 
	print PARAM "T-lex release 1 Detection of the absence of given sequence(s) in strain(s) using SHRIMP program\n"; 
    }
    else {
	print "\n\t\t\t\t\t\t\t * T-lex release 2 \n\t\t\t\tReport the presence/absence of given sequence(s) in strain(s) * \n\t\t\t\t\t\t and return their frequency\n";   
	print PARAM "T-lex release 2 Report the presence/absence of given sequence(s) in strain(s) and return their frequency\n";   
    }
    my $start_time = localtime();
    print "\t\t\t\t\t\t    * $start_time * \n";
    print "\n **************************************************************************************************************************************\n";

    print "* Output directory: $startdirectory/$outputdir/\n";
  
    print PARAM "List of TEs: $TE_list\nTE annotations: $TE_map\nReference genome sequence: $refgenome\nSolexa data: $strains\n";
    print PARAM "Species: $species\n";

    if ($noFilterTE){
	print PARAM "Filter out the problematic TE insertions from the TE dataset: no\n\n";
    }
    else{
	print PARAM "Filter TEs: yes\n";
	print PARAM "==>Min, flanking coverage: $minflankcov\n";
	print PARAM "==>Max. read length: $maxReadLength bp \n\n";
     } 
    print PARAM "TE presence detection:\n\==>Length of the TE: $junction bp\n==>Length of the TE sequence: $buffer bp\n==>Min match length require in the TE sequence: $limp bp\n==>Min sequence identity inside the TE sequence: $id \% \n==>Min quality score: $minqual\n\n";
    print PARAM "TE absence detection:\n==>Length of the flanking sequence: $flank bp\n==>Min length of match on the two sides of the TE: $var bp\n==>Min. non-repeated mapped region on each side of the sequence:$lima bp\n";
    
    if ($PE =~/yes/){
	print PARAM "==>Paired-end mapping\n\n";
    }
    else{
	print PARAM "==>Single-end mapping\n\n";
    }
    close PARAM;

	
	#########################################################################################################################################################
	############################################################## PREPARE THE INPUT DATA ################################################################### 
	#########################################################################################################################################################
	
	print "\n\n\t\t\t\t ******************** Prepare the input data ********************\n\n";
	
	print "\nSimplify the fasta file of the reference sequences ...\n";
	&simplify_fasta("$refgenome","$outputdir\/Tgenome\.fasta");
	
	my $flank_check=$flank;
	my $flank_coords_check;
	
	print "\n\n\t\t\t\t ******************** Tjunction analysis ********************\n\n";
	system("mkdir $outputdir\/Tanalysis");
	
	# Identification of TE insertions nested or flanked by repeats
	print "Identification of TE insertions nested or flanked by repeats....\n";
	$flank_coords_check = &extract_flanks("$TE_list","$TE_map","$flank_check","0","no");
	open (OUT, ">$outputdir\/Tflank_checking\_$flank_check\.map");
	print OUT $flank_coords_check;
	close OUT;
	&maptodb("$outputdir\/Tflank_checking\_$flank_check\.map","$outputdir\/Tgenome\.fasta","$outputdir\/Tflank_checking\_$flank_check\.fasta");	
	system("RepeatMasker -species $species $outputdir\/Tflank_checking\_$flank_check\.fasta");
	system("mv $outputdir\/Tflank_checking\_$flank_check\.fasta\.out $outputdir\/Tflank_checking\_$flank_check\.fasta\.masked  $outputdir\/Tanalysis");
	$flank_coords_check = &extract_flanks("$TE_list","$TE_map","$flank_check","0","yes");
	open (OUT1, ">$outputdir\/Tpoly\_$flank_check\.map");
	print OUT1 $flank_coords_check;
	close OUT1;
	&maptodb("$outputdir\/Tpoly\_$flank_check\.map","$outputdir\/Tgenome\.fasta","$outputdir\/Tpoly\_$flank_check\.fasta");
	system("cp $outputdir\/Tpoly\_$flank_check\.map $outputdir\/Tpoly\_$flank_check\.fasta  $outputdir\/Tanalysis");
	
        # Identification of TE insertions misannotated because of a longer Poly A/T tail
	print "\nIdentification of TE insertions misannotated because of a longer Poly A/T tail....\n";	
	&polyAT_tail_detection("$outputdir\/Tanalysis", "$TE_map", "$TE_list","$outputdir\/Tgenome\.fasta", "$flank_check");
	

	# Identification of TE insertions part of segmental duplications
	print "\nIdentification of TE insertions part of segmental duplications....\n";
	&segmental_dup_detection("$outputdir\/Tanalysis", "$TE_map", "$TE_list", "$outputdir\/Tflank_checking\_$flank_check\.fasta", "$outputdir\/Tgenome\.fasta","$maxReadLength","$id");
	system("mv $outputdir\/Tflank_checking\_$flank_check\.fasta\.* $outputdir\/Tanalysis");
	
	if ($noFilterTE) {
	    print "\t\t\t\t\t\t\t*******************TE FILTER step is bypassed*********************\n"; 
	}
	else{
	    my $start_time_TEfilter = localtime();
	    print "\t\t\t\t\t\t\t*******************FILTER TEs starts at $start_time_TEfilter*********************\n"; 
	    chdir("$outputdir\/Tanalysis");
		my $TE_list_dir = dirname($TE_list);

		if (index($TE_list_dir, "/") != -1) {
			open (IN, "$TE_list") or die $!;
			$TE_list =~ /.*\/(.*)$/;
            $new_TE_list = $1
		} else {
			open (IN, "$startdirectory\/$TE_list") or die $!;
			$new_TE_list = $TE_list;
			$TE_list_dir = $startdirectory;
		}
	    open (IN2, "Tflank_checking\_$flank_check\.fasta\.out") or die $!;
            # $TE_list =~ /.*\/(.*)$/;
            # $new_TE_list = $1;
	    open (OUT, ">$TE_list_dir\/$new_TE_list\_filtered");
	    open (OUT2, ">$TE_list_dir\/$new_TE_list\_FRD");
	    my $n1=0;
	    my $n2=0;
	    
	    while(<IN>){
		$n2=0;
		my $TE=$_;
		chomp $TE;
		my @positions_RIGHT;
		my @positions_LEFT;
		my $totcov_LEFT=0;
		my $totcov_RIGHT=0;
		my $density_LEFT=0;
		my $density_RIGHT=0;
		my $tot=0;
		
		open (IN2, "Tflank_checking\_$flank\.fasta\.out") or die $!;
		while(<IN2>){
		    if ($n2>2){
			my ($query_name, $query_start, $query_end, $left ) = (split)[4,5,6,7];
			my @name = split(/_/,$query_name);
			
			if ($query_name =~ /^$TE/){
			    if ($query_end < $query_start){
				$query_end = $query_start;
				$query_start =  $query_end;
			    }
			    if ($name[1] =~ /LEFT/){
				push @positions_LEFT, [$query_start, $query_end];
			    }
			    if ($name[1] =~ /RIGHT/){
				push @positions_RIGHT, [$query_start, $query_end];
			    }
			}
		    }
		    $n2++;
		}
		close(IN2);
		$n1++;

		if ($#positions_LEFT != -1) {
		    @positions_LEFT = sort {$a->[0] <=> $b->[0]} @positions_LEFT;
		    my $prec_start = $positions_LEFT[0][0];
		    my $prec_end = $positions_LEFT[0][1];
		    for my $i ( 1 .. $#positions_LEFT ) {
			if ($prec_end > $positions_LEFT[$i][0]){ 
			    if($prec_end < $positions_LEFT[$i][1]){
				print "Upating of the data for the TE $TE ...\n";	
				$positions_LEFT[$i-1][1] = $positions_LEFT[$i][1];
				$prec_start = $positions_LEFT[$i-1][0];
				$prec_end = $positions_LEFT[$i-1][1];
				splice(@positions_LEFT,$i,1);
			    }
			}			
		    }
		    
		    if ($#positions_LEFT>0){
			for my $j ( 0 .. $#positions_LEFT ) {
			    my $repeat_length = $positions_LEFT[$j][1] - $positions_LEFT[$j][0];
			    $repeat_length = $repeat_length + 1;
			    $totcov_LEFT= $totcov_LEFT + $repeat_length ;
			}
		    }
		    else{
			my $repeat_length = $positions_LEFT[0][1] - $positions_LEFT[0][0];
			$repeat_length = $repeat_length + 1;
			$totcov_LEFT= $totcov_LEFT + $repeat_length ;
		    }
		}
		
		if ($#positions_RIGHT != -1) {
		    @positions_RIGHT = sort {$a->[0] <=> $b->[0]} @positions_RIGHT;
		    my $prec_start = $positions_RIGHT[0][0];
		    my $prec_end = $positions_RIGHT[0][1];
		    for my $i ( 1 .. $#positions_RIGHT ) {
			if ($prec_end > $positions_RIGHT[$i][0]){ 
			    if($prec_end < $positions_RIGHT[$i][1]){
				print "Ovlerlap: upating of the data ...\n";	
				$positions_RIGHT[$i-1][1] = $positions_RIGHT[$i][1];
				print "\tAFTER UPDATE: precedent = $positions_RIGHT[$i-1][0]\-$positions_RIGHT[$i-1][1]\n";
				$prec_start = $positions_RIGHT[$i-1][0];
				$prec_end = $positions_RIGHT[$i-1][1];
				splice(@positions_RIGHT,$i,1);
			    }
			}			
		    } 		    
		    if ($#positions_RIGHT>0){
			for my $j ( 0 .. $#positions_RIGHT ) {
			    my $repeat_length = $positions_RIGHT[$j][1] - $positions_RIGHT[$j][0];
			    $repeat_length = $repeat_length + 1;
			    $totcov_RIGHT= $totcov_RIGHT + $repeat_length ;
			}
		    }
		    else{
			my $repeat_length = $positions_RIGHT[0][1] - $positions_RIGHT[0][0];
			$repeat_length = $repeat_length + 1;
			$totcov_RIGHT= $totcov_RIGHT + $repeat_length ;
		    }
		}
		
		$density_LEFT=sprintf("%.2f",($totcov_LEFT/($flank)));
		$density_RIGHT=sprintf("%.2f",($totcov_RIGHT/($flank)));
		if ($density_LEFT < $minflankcov && $density_RIGHT < $minflankcov){
		    print OUT "$TE\n";
		}
		print OUT2 "$TE\t$totcov_LEFT\t$totcov_RIGHT\t$density_LEFT\t$density_RIGHT\n";
		
	    }
	    close OUT;
	    close OUT2;
	    
	    $TE_list="$TE_list_dir\/$new_TE_list\_filtered";
	    print "* The new TE list is stored in the file : $TE_list\n";
       
	    chdir("..\/..\/");
	}
	
	unless ($shrimponly) {	
	    my $start_time_BWA = localtime();
	    print "\n\n\t\t\t\t   ******************** Launch Presence detection start at $start_time_BWA ********************\n\n"; 
	    print "Parameters for the detection of the *PRESENCE* of the given sequence(s):\n";
	    print "   - Length of the flanking region = $junction\n";
	    print "   - Length of the internal region = $buffer\n";
	    
	    print "Convert the TE coordinates $TE_list ...\n";
	    my $flank_coords = &extract_flanks("$TE_list","$TE_map","$junction","$buffer","no");	    
	    open (OUT, ">$outputdir\/Tcoords_junction\_$junction\_$buffer\.map");
	    print OUT $flank_coords;
	    close OUT;
	    
	    print "Extract the TE junctions ....\n";
	    &maptodb("$outputdir\/Tcoords_junction\_$junction\_$buffer\.map","$outputdir\/Tgenome\.fasta","$outputdir\/Tjunction\_$junction\_$buffer\.fasta");
	    system("mkdir $outputdir\/Tpresence");
		if ($binreads) {
			$reads = 1;
		}
	    &PresenceDetectionMultipleStrains("$strains","$reads","$outputdir\/Tpresence","$maxReadLength","$processes","$junction","$buffer","$limp","$id","$minqual","$TE_map","$refgenome");  
	    system("cp $outputdir\/Tpresence\/results $outputdir\/Tresults_presence");
	    system("rm $outputdir\/Tpresence\/results");
	    
	    my $end_time_BWA = localtime();
	    print "\n\n\t\t\t\t   ******************** Presence detection end at $end_time_BWA ************\n\n"; 
	}
	
	unless ($bwaonly){ 
	    my $start_time_SHR = localtime();
	    print "\n\n\t\t\t\t   ******************** Launch Absence detection start at $start_time_SHR********************\n\n";
	    print "Parameters for the detection of the *ABSENCE* of the given sequence(s):\n";
	    print "   - Length of the flanking region = $flank\n";
	    print "   - Length of the internal region = 0\n";
	    
	    print "Convert the TE coordinates ...\n";		
	    my $flank_coords = &extract_flanks("$TE_list","$TE_map","$flank","0","no");
	    open (OUT, ">$outputdir\/Tcoords_flank\_$flank\.map");
	    print OUT $flank_coords;
	    close OUT;
	    
	    print "Extract the TE flanked regions ....\n";
	    &maptodb("$outputdir\/Tcoords_flank\_$flank\.map","$outputdir\/Tgenome\.fasta","$outputdir\/Tflank\_$flank\.fasta");

	    open (IN, "$outputdir\/Tflank\_$flank.fasta");
	    open (OUT, ">$outputdir\/Tconcat\_$flank\.fasta");
	    while (<IN>) {
		unless ($_ =~ /_RIGHT$/) {
		    unless ($_ =~ /_LEFT$/) {
			print OUT $_;
		    }
		    if ($_ =~ /(.*)_LEFT$/) {
			print OUT "$1\n";
		    }
		}
	    }
	    
	    close IN;
	    close OUT;

	    system("mkdir $outputdir\/Tabsence");
	    &AbsenceDetectionMultipleStrains("$strains","$outputdir\/Tabsence","$TE_list","$outputdir\/Tconcat\_$flank\.fasta","$reads","$lima","$var","$flank","$PE","$id");
	    system("cp $outputdir\/Tabsence\/results $outputdir\/Tresults_absence");
	    system("rm $outputdir\/Tabsence\/results");
	    
	    my $end_time_SHRIMP = localtime();
	    print "\n\n\t\t\t\t   ******************** Absence detection end at $end_time_SHRIMP ************\n\n";
	}	

   
    if ($noFilterTE){
        $newTElist = 1;
        $TE_list=$original_TE_list;
        print "TE list NO cleaned\n";
    }
    else{
        print "TE list cleaned\n";
    }

    #############
    my $start_TSD = localtime();
    print  "* TSD detection process start at $start_TSD\n";
    if($tsd){
        &MultiAlign("$outputdir","$TE_list","$strains","$flank","0","1");
    }
    my $end_TSD = localtime();
    print  "* TSD detection end at $end_TSD\n"; 	
    
    #############
    my $start_align = localtime();
    print  "* Multiple Alignment process start at $start_align\n";
    unless ($shrimponly) {
        &MultiAlign("$outputdir","$TE_list","$strains","$flank","1","0");
    }
    unless ($bwaonly) {
        &MultiAlign("$outputdir","$TE_list","$strains","$flank","0","0");
    }
    else{
        &MultiAlign("$outputdir","$TE_list","$strains","$flank","2","0");
    } 	
    my $end_align = localtime();
    print  "* Multiple Alignment end at $end_align\n"; 	
    
    unless ($shrimponly || $bwaonly){
        &FinalResults("$outputdir","$TE_list", "$TE_map","$flank");
        &FreqEstimate("$outputdir","$pooleddata","$maxreads","$minreads","$minpop");
    }

    ####################
    unless ($noclean) {
        print "cleaning......";
        system("rm -rf $outputdir\/Tconcat_*.fasta $outputdir\/Tcoords* $outputdir\/Tflank* $outputdir\/Tgenome* Tjunction* $outputdir\/tmp*");
        system("rm -rf $outputdir\/Tpresence $outputdir\/Tabsence");
    }
    my $end_time_Tlex = localtime();
    print "\n\n\t\t\t\t   ******************** T-lex finished successfully at $end_time_Tlex ************\n\n";
    print "\n\n\t\t\t\t   ************************ Have a nice day! ************************n\n";
	
}


#########################################################################################################################################################
########################################################################### HELP ######################################################################## 
########################################################################################################################################################
	
sub help {
	print "T-lex\ release 3\n\n";
	
    print "T-lex3 works in modules (Presence, Absence, Combine, TSD detection, Frequency estimation) and each module requires a T-lex run\n\n";
	print "\t- Usage:\n";
	print "\t  ======\n\n"; 
	print"\t  tlex-open-v3.pl [ options ] [ -T TE_list ] [ -M TE_annotations ] [ -G reference_genome ] [ -R Next-Generation Sequencing (NGS) data ]\n\n";
	print "\t\-T \t\tfile with the list of the transposable element (TE) identifiers\n";
	print "\t\-M \t\ttabulated file with the TE annotations (line format : \'TE name location start postion end position\')\n";
	print "\t\-G \t\tfile of the reference genomic sequences in FASTA\n";
	print "\t\-R \t\tdirectory of the NGS data in FASTQ\n\n\n\n";

	print "\t- Options:\n";
	print "\t  ========\n\n";
	print "\t\t\-A \t\tmax. read length in the data set ( default : 100 bp )\n";
	print "\t\t\-O \t\tproject name\( by default the output directory will be \"tlex_output\"\)\n";
	print "\t\t\-noclean \tkeep the intermediate files\n";
	print "\t\t\-h or -help \tdisplay this help\n\n";

	print "\t\t*For the analysis of the flanking sequences of each annotated TE: \n";
	print "\t\t\-noFilterTE \tdo not filter TEs\n";
	print "\t\t\-s \t\tname the species studied ( default: drosophila, cf. RepeatMasker program)\n";
	print "\t\t\-d \t\tmin. repeat density at the flanking regions of the TEs ( default : 0.5, that corresponds to a repeat density of 50 %)\n\n";
	print "\t\t\-id \t\tmin. sequence mapping identity required with the TE sequence in % ( default: 95 )\n"; 
	
	print "\t\t*For the TE presence detection: \n";
	print "\t\t\-q \t\tlaunch only the presence detection approach\n";
	print "\t\t\-j \t\tlength of the junction sequences to extract in bp ( default: 1000 )\n";
	print "\t\t\-b \t\tlength of the internal region of the TE in bp ( default: 60 )\n";
	print "\t\t\-limp \t\tmin. length of match required with the TE sequence in bp ( default: 15 )\n";  
	print "\t\t\-minQ \t\tmin. quality Phred score for the read assembly ( default: 30 )\n"; 
	print "\t\t\-processes \tnumber of processes (only used for the NGS data reformatting; default: 1)\n\n";
    print "\t\t\-pairends \tpaired-end mapping ('yes' or 'no'; default: 'no'; that option requires NGS data such as <strain name>_reads<1 or 2>.fastq)\n\n";
	
	print "\t\t*For the TE absence detection: \n"; 
	print "\t\t\-p \t\tlaunch only the absence detection approach\n";
	print "\t\t\-f \t\tlength of the flanking sequences to extract and concatenate in bp ( default: 125 )\n";
	print "\t\t\-v \t\tmin. read length spanning the two TE sides in bp ( default: 20 )\n";
	print "\t\t\-lima \t\tmin. non-repeated region on each side of the sequence in bp ( default: 15 )\n";
	print "\t\t\-pairends \tpaired-end mapping ('yes' or 'no'; default: 'no'; that option requires NGS data such as <strain name>_reads<1 or 2>.fastq)\n\n";
	
	print "\t\t*For manually combine the results:\n";
    print "\t\t\-combRes \tcombine the presence/absence results from one dataset\n";
	print "\t\t\-combData \tcombine the presence/absence results from several datasets stored in a different tlex output directories\n";
	print "\t\t\-combAll \tcombine the frequency estimates, the analysis of the TE flanking regions and  the TSD detection\n";
    print "\t\t*Estimate the TE frequency: \n";
    print "\t\t\-freq \t\treturn the TE frequency based on the given strains\n";
    print "\t\t\-minP \t\tminimum number of TE data based on given strains ( default: 1 )\n";
	print "\t\t\-pooled \treturn the TE frequency based on pooled data (To use with the option -freq)\n";
	print "\t\t\-minR \t\tminimum number of reads to calculate frequency ( default: 3 )\n";
    print "\t\t\-maxR \t\tmaximum number of reads to calculate frequency ( default: 90 )\n";
    print "\t\t*For TSD and multialignment detection:\n";
    print "\t\t-align \t\treturn the multiple alignments\n"; 
	print "\t\t\-tsd \t\treturn the Target Site Duplication (TSD) for each given TE insertion detected as absent (use with -align & -p)\n";	
	
	die "\n";
 }	


#########################################################################################################################################################
################################################################ MULTIPLE ALIGNMENTS ####################################################################
#########################################################################################################################################################

sub MultiAlign(){

    my $outputdir = $_[0]; 
    my $TE_list = $_[1]; 
    my $strains = $_[2];
    my $flank = $_[3];  
    my $type = $_[4];  
    my $TSD = $_[5]; 

    print "TYPE: $type\n";
    print "Output directory: $outputdir\n\n";
    chdir("$outputdir"); 
    print "mkdir Talign\n";
    system ("mkdir Talign");
    
    if( $type == 1 || $type == 2){
	print "PRESENCE ALIGNMENT\n";
	if (-d "Talign\/presence_detection\/"){
	    print "Talign\/presence_detection directory created !\n";
	}
	else{
	    print "presence_detection directory does not exist !\n";
	    system("mkdir Talign\/presence_detection");
	}
	system("pwd");
	open (INT, "$TE_list") or die $!; 
	while (<INT>) {
	    my $TE = $_;
	    chomp $TE;
	    opendir (DIR, "$strains");
	    my $nline=0;
	    while (defined( my $strain_name = readdir (DIR))) {
		if ($strain_name !~ /\./) {
		    if ($nline == 0){
			$nline++;
			if(-e "Tpresence\/${strain_name}\/detection\/${TE}\_LEFT\.ref"){
			    system ("cat Tpresence\/${strain_name}\/detection\/${TE}\_LEFT\.ref  > Talign\/presence_detection\/${TE}\_LEFT\.contig_ref");
			}
			if(-e "Tpresence\/${strain_name}\/detection\/${TE}\_RIGHT\.ref"){
			    system ("cat Tpresence\/${strain_name}\/detection\/${TE}\_RIGHT\.ref > Talign\/presence_detection\/${TE}\_RIGHT\.contig_ref");
			}
		    }
		    if ($nline == 1){
			if (-e "Talign\/presence_detection\/${TE}\_LEFT\.contig_ref"){
			    system ("cat Talign\/presence_detection\/${TE}\_LEFT\.contig_ref Tpresence\/\*\/detection\/${TE}\_LEFT\.contig > Talign\/presence_detection\/${TE}\_LEFT\.contig_all");
			}
			if (-e "Talign\/presence_detection\/${TE}\_RIGHT\.contig_ref"){
			    system ("cat Talign\/presence_detection\/${TE}\_RIGHT\.contig_ref Tpresence\/\*\/detection\/${TE}\_RIGHT\.contig > Talign\/presence_detection\/${TE}\_RIGHT\.contig_all");	    
			}
			last;
		    }
		}
	    }
	}
	system ("rm Talign\/presence_detection\/*.contig_ref");
	close IN ;
	closedir DIR;
    }
    
    if( $type == 0 || $type == 2){
	print "ABSENCE ALIGNMENT\n";
       
	my %Selected=();
	my %SelectedM=();
	my %dejavu=();
	
	if (-d "Talign\/"){
	    print "Talign\/ directory exists !\n";
	}
	else{
	    system("mkdir Talign\/");
	}
	
	if (-d "Talign\/absence_detection\/"){
	    print "Talign\/absence_detection directory exists !\n";
	}
	
	else{
	    print "absence_detection directory does not exist !\n";
	    system("mkdir Talign\/absence_detection");
		    
	    opendir (DIR, "$strains");
	    while (defined( my $strain_name = readdir (DIR))) {
		if ($strain_name !~ /\./) {
		    if (-e "Tabsence\/${strain_name}/detection/results_gmapper_overjunction_nodup_format_nogap.fa.masked"){
			    system ("cp Tabsence\/${strain_name}\/mapping/results_gmapper2merge Talign\/absence_detection\/Talign_${strain_name}.fasta");
                system ("cp Tabsence\/${strain_name}\/detection\/results_gmapper_overjunction Talign\/absence_detection\/Talign_${strain_name}.merged");
			    system ("cp Tabsence\/${strain_name}\/detection\/results_gmapper_overjunction_nodup_format2select_format Talign\/absence_detection\/TselectedRead_${strain_name}.list");
		     }
		     else{
			 print "no file called Tabsence\/${strain_name}/detection/*_nogap.fa.masked";
			 next;
		     }
		 }
	     }
	     close DIR ;
	    
	    opendir (DIR, "$strains"); # Removed ../
	    while (defined( my $strain_name = readdir (DIR))) {
		if ($strain_name !~ /\./) {
		    print "\n\n ### $strain_name : \n";
		    my $read_headers;
		    my $TE_headers;
		    my $readseq;
		    my $refseq;
		    my $nickname;
		    my @fullname;
		    my @shortname;
		    my $pe;
		    my $me;
		    my $r1;
		    my $r2;
		    my $preTE="";
		    my $prereadt="";
		    
		    open (SE, "Talign\/absence_detection\/TselectedRead_${strain_name}.list");
		    while (<SE>) {
			my ($TE, $read)=(split)[0,1];
			my $ntype = " ";
			my @shortr=split("\/",$read);
			
			if($shortr[1]){
			    $ntype = $shortr[1];
			}
			else{
			    $ntype = "merge";
			}
			
			$Selected{$TE}{$read} = $ntype;
		    }
		    close SE;


		    open (SEF, "Talign\/absence_detection\/Talign_${strain_name}.merged");
		    while (<SEF>) {
			my ($read, $TE, $ref, $seq)=(split)[0,1,9,10];
			my $ntype = " ";
			my $nickm;
			my @fullm=split("\:",$read);
			if ( scalar(@fullm)>1 ) {
			    shift(@fullm); 
			    $nickm=join(":",@fullm);
			}
			
			my @shortm=split("\/",$nickm);
			
			if($shortm[1]){
			    $ntype = $shortm[1];
			}
			else{
			    $ntype = "merge";
			}

			if( $ntype =~ /merge/){
			    $SelectedM{$TE}{$nickm} = [$ref,$seq];
			}
		    }
		    close SEF;
		    
		    my $preTE_headers = "";
		    my $prenewname = "";
		    my $preme = "";
		    my $newname = "";

		    open (FA, "Talign\/absence_detection\/Talign_${strain_name}.fasta"); 
		    while(<FA>){
			my @l = split(); 
			my $lg = scalar(@l);  
		
			if ( $_ =~ /^\>/ && $lg >1 ) {
			    my @line = split(); 
			
			    $read_headers=$line[0];
			    $TE_headers=$line[1];
			    chomp $read_headers;
			    
			    
			    @fullname=split("\:",$read_headers);
			    if ( scalar(@fullname)>1 ) {
				shift(@fullname); 
				$nickname=join(":",@fullname);
			    }
			    elsif  ( scalar(@fullname)==1 ) {
				@fullname=split("\>",$read_headers);
				shift(@fullname); 
				$nickname=join("",@fullname);
			    }
			    else{
				print "error read name";
				die $!;
			    }
			    
			    @shortname=split("\/",$nickname);
			    $pe = $shortname[0];
			    $me = $shortname[1];
			}
			
			if ($_ =~ /^G/) {
			    $refseq= (split)[2];
			    chomp $refseq;			   
			}
			if ($_ =~ /^R/) {
			    $readseq= (split)[2];
			    chomp $readseq;
     		 
			    $newname =  $pe."/merge";
			    my $output_ref;
			    my $output_read;

			    if (exists $Selected{$TE_headers}{$nickname}){
				my $output.=">$TE_headers\n$refseq\n$read_headers\n$readseq\n";
				open (OUT, ">tmp_all");
				print OUT $output;
				close OUT;
				
				$output_ref.=">${TE_headers}_$nickname\n$refseq\n";
				open (OUTREF, ">tmp_ref");
				print OUTREF $output_ref;
				close OUTREF;	
				
				$output_read.=">${nickname}\n$readseq\n";
				open (OUTREAD, ">tmp_reads");
				print OUTREAD $output_read;
				close OUTREAD;
				
				unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}.fasta" ) {
				    system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}.fasta");
				}
				system(" cat tmp_all >> Talign\/absence_detection\/${TE_headers}_${strain_name}.fasta");
				
				unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}_ref.fasta" ) {
				    system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}_ref.fasta");
				}
				system(" cat tmp_ref >> Talign\/absence_detection\/${TE_headers}_${strain_name}_ref.fasta");	
				unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}_reads.fasta" ) {
				    system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}_reads.fasta");
				}
				system(" cat tmp_reads >> Talign\/absence_detection\/${TE_headers}_${strain_name}_reads.fasta");
			    }

			    elsif (exists $Selected{$TE_headers}{$newname}){
				if (($TE_headers eq $preTE_headers) && ($newname eq $prenewname) && ( ($me == 1 && $preme == 2) || ($me == 2 && $preme == 1)) ){
				    print "";
				}
				else{ 
				    
				    my @seqs = $SelectedM{$TE_headers}{$newname};
				    $refseq = $seqs[0][0];
				    $readseq = $seqs [0][1];
				    
				    my $output.=">${TE_headers}_$newname\n$refseq\n>${newname}\n$readseq\n";
				    open (OUT, ">tmp_merge");
				    print OUT $output;
				    close OUT;
				     
				    $output_ref.=">${TE_headers}_$newname\n$refseq\n";
				    open (OUTREF, ">tmp_ref_merge");
				    print OUTREF $output_ref;
				    close OUTREF;
				   
				    $output_read.=">${newname}\n$readseq\n";
				    open (OUTREAD, ">tmp_reads_merge");
				    print OUTREAD $output_read;
				    close OUTREAD;
				    
				    unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}.fasta" ) {
					system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}.fasta");
				    }
				    system(" cat tmp_merge >> Talign\/absence_detection\/${TE_headers}_${strain_name}.fasta");

				    unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}_ref.fasta" ) {
					system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}_ref.fasta");
				    }
				    system(" cat tmp_ref_merge >> Talign\/absence_detection\/${TE_headers}_${strain_name}_ref.fasta");

				    unless (-e "Talign\/absence_detection\/${TE_headers}_${strain_name}_reads.fasta" ) {
					system(" touch Talign/absence_detection\/${TE_headers}_${strain_name}_reads.fasta");
				    }
				    system(" cat tmp_reads_merge >> Talign\/absence_detection\/${TE_headers}_${strain_name}_reads.fasta");
				}
			    }
			    $preTE_headers = $TE_headers;
			    $prenewname = $newname;
			    $preme = $me;
			}		
		    }
		    close FA;
		}
	    }
	    close DIR;	
	}
	
	
	my %SelectedTSD=();
	if($TSD==1){ 
	    opendir (DIR, "$strains");
	    while (defined( my $strain_name = readdir (DIR))) {
		if ($strain_name !~ /\./) {
		    open (SE, "Talign\/absence_detection\/TselectedRead_${strain_name}.list") or die $!;		
		    while (<SE>) {
			my ($TE, $read)=(split)[0,1];
			my $readt = " ";
			my $type = " ";
			my @shortr=split("\/",$read);
			$readt = $shortr[0];
			
			if($shortr[1]){
			    $type = $shortr[1];
			}
			else{
			    $type = "merge";
			}
			
			$Selected{$TE}{$readt}{$type} = "selected";
		    }
		    close SE;
		}
	    }
	    
	    print "TARGET SITE DETECTION\n\n";
	    open (TTSD, ">Tannot_TSD");
	    print TTSD "TE\tputative_target_site\tTSD\tTSDdistance\tTSDresult\treference_sequence\tcontig_sequence\n"; 
	    for my $TE_headers ( keys %Selected ) {
		print "\n####  $TE_headers\n";
		system("cat Talign\/absence_detection\/${TE_headers}_*_reads.fasta > Talign\/absence_detection\/${TE_headers}_reads.fasta");
		system("cat Talign\/absence_detection\/${TE_headers}_*_ref.fasta > Talign\/absence_detection\/${TE_headers}_ref.fasta");
		my $nbreads=&fastacount("Talign\/absence_detection\/${TE_headers}_reads.fasta");
		if ($nbreads <= 3) {
		    print " *** no enough reads to build a contig ==>  $nbreads reads  ***\n";
		    # use the single reads to detect the TSD
		    
		    # 1/ split the fasta files : reads + ref
		    system("split -d -l 2 Talign\/absence_detection\/${TE_headers}_reads.fasta Talign\/absence_detection\/${TE_headers}_reads.fasta");
		    system("split -d -l 2 Talign\/absence_detection\/${TE_headers}_ref.fasta Talign\/absence_detection\/${TE_headers}_ref.fasta");
		    # 2/ Blat each 
		    if ($nbreads > 0) {
			system("blat Talign\/absence_detection\/${TE_headers}_reads.fasta00 Talign\/absence_detection\/${TE_headers}_ref.fasta00 Talign\/absence_detection\/${TE_headers}_reads.fasta00_blat -out=pslx &> stdoutput");
			system("cat Talign\/absence_detection\/${TE_headers}_reads.fasta00_blat >> tmp_reads_blat");
			if ($nbreads > 1) {
			    system("blat Talign\/absence_detection\/${TE_headers}_reads.fasta01 Talign\/absence_detection\/${TE_headers}_ref.fasta01 Talign\/absence_detection\/${TE_headers}_reads.fasta01_blat -out=pslx &> stdoutput");
			    system("cat Talign\/absence_detection\/${TE_headers}_reads.fasta01_blat >> tmp_reads_blat");
			    if ($nbreads > 2) {
				system("blat Talign\/absence_detection\/${TE_headers}_reads.fasta02 Talign\/absence_detection\/${TE_headers}_ref.fasta02 Talign\/absence_detection\/${TE_headers}_reads.fasta02_blat -out=pslx &> stdoutput");
				system("cat Talign\/absence_detection\/${TE_headers}_reads.fasta02_blat >> tmp_reads_blat");
			    }
			}
		    }
		    else{
			
		    }
		    my $l;
		    open (BLATIN, "tmp_reads_blat") or die $!; 
		    open (BLATOUT, ">Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_blat") or die $!; 
		    while(<BLATIN>){
			chomp $_;
			if (($_ =~ /^$/ )){
			    $l = "empty";
			}
			elsif (($_ =~ /^psLayout/ )){
			    $l = "psl";
			}
			elsif (($_ =~ /^\-/ )){
			    $l = "dashes";
			}
			elsif (($_ =~ /match/ )){
			    $l = "headers";
			}
			else{
			    print BLATOUT $_;
			}
		    }

		    close BLATIN;
		    close BLATOUT;
		    system("cp Talign\/absence_detection\/${TE_headers}_reads.fasta Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs");
		    system("cp Talign\/absence_detection\/${TE_headers}_ref.fasta Talign\/absence_detection\/${TE_headers}_ref.fasta.contigs");
		}

		else{ 
		    if (-s "Talign\/absence_detection\/${TE_headers}_reads.fasta"){
			system("sed s/'\-'/''/g Talign\/absence_detection\/${TE_headers}_reads.fasta > Talign\/absence_detection\/${TE_headers}_reads_nogap.fasta");
			system("phrap Talign\/absence_detection\/${TE_headers}_reads.fasta &> stdoutput ");
		    }
		    if (-s "Talign\/absence_detection\/${TE_headers}_ref.fasta"){
			system("sed s/'\-'/''/g Talign\/absence_detection\/${TE_headers}_ref.fasta > Talign\/absence_detection\/${TE_headers}_ref_nogap.fasta");
			system("phrap Talign\/absence_detection\/${TE_headers}_ref.fasta &> stdoutput");
		    }
		}
 
		if ((-s "Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs") && (-s "Talign\/absence_detection\/${TE_headers}_ref.fasta.contigs")){
		    my $nbcontigs=&fastacount("Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs");
		    system("blat Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs Talign\/absence_detection\/${TE_headers}_ref.fasta.contigs Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_blat -out=pslx &> stdoutput");
		}
		
		if (-s "Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_blat"){
		    my ($match, $mismatch, $Ns, $strand);
		    my $refseq="NA";
		    my $readseq="NA";
		    my $l;
		    
		    
		    open (BLATT, "Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_blat") or die $!;
		    while(<BLATT>){
			chomp $_;
			if (($_ =~ /^$/ )){
			    $l = "empty";
			}
			elsif (($_ =~ /^psLayout/ )){
			    $l = "psl";
			}
			elsif (($_ =~ /^\-/ )){
			    $l = "dashes";
			}
			elsif (($_ =~ /match/ )){
			    $l = "headers";
			}
			else{
			    my ($match, $mismatch, $Ns, $strand, $refseq, $readseq) = (split)[0,1,3,8,21,22]; 
			    my $new_readseq=(split("\,",$readseq))[0];
			    my $id = $mismatch/$match;    
			    my @gaps;
			    my $gap="";
			    my $start=-1;
			    my $end=-1;
			    my $count=1;
			    my $prebase="x";
			    my @seqsplit= split("",$new_readseq);

			    for my $base (@seqsplit) { 
				if ( $prebase =~ /[atcg]/ && $base =~ /n/){
				    $start=$count;
				}
				elsif ( $prebase =~ /n/ && $base =~ /[atcg]/){
				    $end=($count-1);
				    if ($start != -1){
					if ($Ns >= 3){
					    print "Gap coordinates: [ $start - $end ] / Length=$Ns\t";
					    my $values=$start." ".$end;
					    push @gaps, $values;
					}
				    }
				}
				$prebase=$base;
				$count=$count+1;
			    }
			    
			    
			    if(scalar(@gaps) == 0){
				print " *** No_Gap_detected *** \n"; 
				print TTSD "${TE_headers}\tNA\tNA\tNA\tnoGap\t$refseq\t$readseq\n"; 
			    }
			    
			    for my $g (@gaps) {
				my @val=split(" ",$g);
				my $GS = $val[0];
				my $GE = $val[1];
				my $GL = ($GE-$GS)+1;
				$gap=substr($refseq, ($GS-1), $GL);
				if( $gap ne " "){
				    my $nbmism = $Ns-int($Ns*(0.80));
				    print "\nMotif=$gap \t Max number of mismatches=$nbmism\t";
				    if (-s "Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs"){
					system("fastagrep -ga $GL $nbmism 0 -needle $gap -haystack Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs | sed s/'\:'/''/ > Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_fastagrep");
				    }
				  
				 		    
				    open (TSD, "Talign\/absence_detection\/${TE_headers}_reads.fasta.contigs_fastagrep") or die;
				    my @line;
				    my @coord;
				    my $NS;
				    my $NE;
				    my @motif;
				    my $copy;
				    my $dist;
				    my $pre_dist=($flank/2);
				    my @cand_arrays;
				    my @cand;
				    my $Ndist;
				    my $Nstart;
				    my $Nend;
				    my $Ncopy;
				    
				    while(<TSD>){
					chomp;
					@line = split();
				
					@coord = split(/-/,$line[2]);
					$NE = $coord[1];
					$NS = $NE-$GL; 
					@motif = split(/-/,$line[3]);
					$copy=$motif[1];
					
					if( $GS >= $NE) {
					    $dist = ($GS-$NE);
					    
					}
					else{
					    $dist = ($NS-$GE);
					    
					}
					
					if( $dist <= $pre_dist && $dist > -2){
					    $pre_dist = $dist;
					    push(@cand_arrays, [$dist, $NS, $NE, $copy]);
					}
				    }
				    
				    if ($#cand_arrays != -1) {
					@cand = sort {$a->[0] <=> $b->[0]} @cand_arrays;
					$Ndist=($cand[0][0]+1); ### T-lex v2.1
					$Nstart=$cand[0][1];
					$Nend=$cand[0][2];
					$Ncopy=$cand[0][3];
					print "==> TSD_welldetected\n";
					print TTSD "${TE_headers}\t$gap\t$Ncopy\t$Ndist\tTSD\t$refseq\t$readseq\n";
					$SelectedTSD{$TE_headers}{$gap}=$Ncopy."/".$readseq; 
				    }
				    else{
					print " ==> TSD_diverge/distant/TE_misannotation\n";
					print TTSD "${TE_headers}\t$gap\tNA\tNA\tnoTSD\t$refseq\t$readseq\n";
					
				    }
				} 
			    } 
			}
		    }
		    close BLATT;
		}
		else{
		    print TTSD "${TE_headers}\tNoblatFile\n";
		} 
	    }
	    close TTSD;
	    system ("sort Tannot_TSD > Tannot_TSDdetection");
	    system ("rm Tannot_TSD");
	}
    }
    
    chdir("..\/");
}

   
#########################################################################################################################################################
###################################################### COMBINE RESULTS FROM DIFFERENT STRAINS ########################################################### 
#########################################################################################################################################################


sub combine_results {
    
    my $nb=0;
    my $path = $_[0]; 
    print "path: $path\n";
    system ("mkdir $path\/Tfreqs_output");
    opendir(DIR, "$path");
    my @repertory = grep(/tlex\_/,readdir(DIR));
    foreach my $item (@repertory) {
        print "$item";
    }

    close DIR;
    
    open(OUT, ">$path\/Tfreqs_output\/Tresults");
    
  FILE: foreach (@repertory) {
      if( -d $_ ){
	  my $nl=0;
	  open(FILE,"${_}/Tresults" ) || ((warn "Can't open file ${_}/Tresults\n"), next FILE);
	  while (<FILE>) {
	      if ($nb>0){
            if($nl>0){
                print OUT $_;
            }
		  $nl++;
	      }
	      else{
		    print OUT $_; 
	      }
	  }
	  close(FILE);
      }
      $nb++;
  }
    
    close OUT;
    
}

 
#########################################################################################################################################################
##################################################################### BEGIN SUBROUTINES ################################################################# 
#########################################################################################################################################################
 

sub extract_flanks { 

    my $TE_list = $_[0];
    my $TE_map=$_[1];
    my $junction=$_[2];
    my $buffer=$_[3];
    my $strand=$_[4];

    open (IN, $TE_list) or die $!;
    my @TEs;
    while (<IN>) {
	chomp;
	push (@TEs,$_);
    }
    my $regex = "";
    foreach (@TEs) {
	$regex .= "\^$_|";
    }
    $regex = substr $regex, 0, -1;
    my $compiledregex = qr/$regex/;
    my $output="";
    close IN;

    open (IN, $TE_map);
    while (<IN>) { 
	
	my $left_flank_end;
	my $right_flank_start;
	my $left_flank_start;
	my $right_flank_end;
	my @lines=split();

	if ($_ =~ $compiledregex) {
	    if ($strand =~ /no/){ 
		$left_flank_start=$lines[2]-$junction-1;
		$right_flank_end=$lines[3]+$junction+1;
		if ($buffer != 0){
		    $left_flank_end=$lines[2]+$buffer-1;
		    $right_flank_start=$lines[3]-$buffer+1;
		}
		else{
		    $left_flank_end=$lines[2]-1;
		    $right_flank_start=$lines[3]+1;
		}
		$output.="$lines[0]".'_LEFT'."\t$lines[1]\t$left_flank_start\t$left_flank_end\n";
		$output.="$lines[0]".'_RIGHT'."\t$lines[1]\t$right_flank_start\t$right_flank_end\n";
	    }
	    else{ 
		$left_flank_start=$lines[2]-$junction-1;
		$left_flank_end=$lines[2];
		$right_flank_start=$lines[3]+1;
		$right_flank_end=$lines[3]+$junction+1;
		$left_flank_end=$lines[2];
		$right_flank_start=$lines[3];
		my $orient=$lines[4];

		$output.="$lines[0]"."_".$orient."_LEFT\t$lines[1]\t$left_flank_start\t$left_flank_end\n";
		$output.="$lines[0]"."_".$orient."_RIGHT\t$lines[1]\t$right_flank_start\t$right_flank_end\n";	
	    }
	}
    }
    return $output;
    close IN;
}

####################################################

sub simplify_fasta {

    my $filename=$_[0];
    my $outfile=$_[1];
    my $output;
 
    open (IN, $filename) or die $!;
    open (OUT, ">$outfile");
    while (<IN>) {
	if ($_ =~ /^>/) {
	    $output = substr $_, 0, index($_, ' ');
	    print OUT "$output\n";
	}
	else {
	    print OUT $_;    
	}
    }   
    close IN;
}

####################################################

sub maptodb { 

    my $flankmap=$_[0];
    my $genome=$_[1];
    my $output=$_[2];
    my $chromosome;
    my $counter=0;
    my %sequences;
   
    open (IN, $genome);
    while (<IN>) {
	if ($_ =~ /^>(.*)/) {
	    $chromosome=$1;
	    chomp $chromosome;
	}
	else {chomp $_;
	      $sequences{$chromosome}.=$_;}
    }

    open (IN2, $flankmap) or die $!;
    open (OUT, ">$output");
    while (<IN2>) {
	my @linesplit = split(/\t/, $_);
	my $flankname=$linesplit[0];
	my $flankchromosome=$linesplit[1];
	my $startposition=($linesplit[2]-1);
	my $junction=($linesplit[3]-$linesplit[2])+1;
	my $flanksequence=substr($sequences{$flankchromosome},$startposition,$junction);
	print OUT ">$flankname\n$flanksequence\n";
    }
    close IN;
    close IN2;
    close OUT;
}

####################################################

sub fastq2fasta {
    print "Convert to fasta and simplify the header...\n";
    my $input = $_[0];
    my $output = $_[1];
    my $split = $_[2];
    
    open (IN, $input) or die "Couldn't find input file!\n";
    open (OUT, ">$output\_1\.fa") or die "Can't get a write handle!\n";
 
    my $linecount=1;
    my $sequencecount=0;
    my $filecount=1;
    my $header;
    
  MAINLOOP: while (<IN>) {
      if ($sequencecount >= $split) {
	  $filecount++;
	  $sequencecount=0;
	  close OUT;
	  
	  if ( ($output =~ /1$/) || ($output =~ /2$/) ){ # <strain>_XXX1 or 2
	      my $outputshort = substr($output, 0, -1); # <strain>_XXX
	  }

	  open (OUT, ">$output\_$filecount\.fa") or die "Can't get a write handle!\n"; # <strain>_XXX<1 or 2>_<count>.fa
      }
	  if ($linecount==1) {
	      s/^.//;
	      $header = substr $_, 0, index($_, ' ');
	      print OUT ">$header\n";
	      $linecount++;
	      next MAINLOOP;
	  }
	  if ($linecount==2) {
	      print OUT $_;
	      $linecount++;
	      next MAINLOOP;
	  }
	  if ($linecount==3) {
	      $linecount++;
	      next MAINLOOP;
	  }
	  if ($linecount==4) {
	      $linecount=1;
	      $sequencecount++;
	      next MAINLOOP;
	  }
  }
    close IN;
}

####################################################

sub fastasplit {

    my $input = $_[0];
    my $output = $_[1];
    my $split = $_[2];
    
    open (IN, $input) or die "Couldn't find input file!\n";
   
    my $linecount=1;
    my $sequencecount=0;
    my $filecount=1;
    
    while (<IN>) {
	if ($sequencecount >= $split) {
	    $filecount++;
	    $sequencecount=0;
	    open (OUT, ">$output\_$filecount\.fa"); # <strain>_XXX<1 or 2>_<count>.fastq
	}
	if ($_ !~ /^>/) {
	    $sequencecount++;
	}
	print OUT $_;
	close OUT;
    }
    close IN;
    close OUTPE;
    close OUTS;
}


#################################################### HERE

sub list_readfiles {
   
    my $dir = $_[0];  
    my $out = $_[1];
    my $type = $_[2];

    print "DIR: $dir\n";
    if( $type =~/yes/){
	print "Identify and create the lists of paired-end files...\n";
	#system ("rm -f $dir\/list_${out}_PE");
	open (OUTPE, ">>$dir\/list_${out}_PE");
    }
    else{
	print "Identify and create the lists of single read files...\n";
	#system ("rm -f $dir\/list_${out}_S");
	open (OUTS, ">>$dir\/list_${out}_S");
    }

    opendir (DIR, "$dir") or die "Failed to open $dir!\n";
 
    my $pe;
    my $nb;

    while (defined( my $fasta = readdir (DIR))) {
	if ($fasta =~ /\.fa$/) {

	    my @readA = split(/\.fa/,$fasta);
	    my $reads = $readA[0];

	    # paired-end
	    if( $type =~/yes/){
		if ( ($reads =~ /1_/) || ($reads =~ /2_/) ){ 
		    my @peA = split(/\_/,$reads);
		    $pe = $peA[0]."_".substr($peA[1], 0, -1);
		    $nb = $peA[2]; 
		    print OUTPE "$pe\t$nb\n";
		}
		else{ 
		    die "Warning: Check the format of the sequecing data (Paired-end or single)\n";
		}
	    }
	    
	    # single reads  
	    else{
		my @sA = split(/\_/,$reads);
		my $c=1;
		my $s=$sA[0];
		foreach (@sA){
		    if( $c > 1 && $c < scalar(@sA)){
			$s=$s."_".$_;
		    }
		    $c++;
		}
		$nb = $sA[scalar(@sA)-1];
		print OUTS "$s\t$nb\n";
	    }
	}
    } 
    if( $type =~/yes/){
	close OUTPE;
    }
    else{
	close OUTS;
    }
}


####################################################

sub fastacount {

    my $input = $_[0];  
    open (IN, $input) or die "Couldn't find input file!\n";
   
    my $sequencecount=0;
    
    while (<IN>) {
	if ($_ =~ /^>/) {
	    $sequencecount++;
	}
    }
    return $sequencecount;
    close IN;
}

####################################################

sub polyAT_tail_detection {
    
    my $outpath = $_[0];
    my $TE_annot = $_[1];
    my $TE_list = $_[2];
    my $refgenome = $_[3];
    my $flank = $_[4];
    
    if ($TE_annot and $refgenome){
	my @linesplit;
	my @seqsplit;
	my @poly;
	my %sequences;
	my $TEid;
	my $seq;
	my $tract;
	my $orient;
	my $side;
	my $fseq;

	open (IN1, "$outpath\/Tpoly\_$flank\.fasta");
	while (<IN1>) { 
	    if ($_ =~ /^>(.*)/) {
		$TEid=$1;
		chomp $TEid;
	    }
	    else {
		$seq=$_;
		chomp $seq;
		$sequences{$TEid}=$seq;
	    }
	}
	
	open (IN2, "$outpath\/Tpoly\_$flank\.map");
	open (OUT2, ">$outpath\/Tpoly\_$flank\.fasta\.polyAT");

	while (<IN2>) {
	    chomp $_;
	    my @lines=split("\_",$_);
	    $TEid=$lines[0];
	    $orient=$lines[1];
	    my @pos=split("\t",$lines[2]);
	    $side=$pos[0];
	    my $chr=$pos[1];
	    my $start=$pos[2];
	    my $end=$pos[3];
	    
	    if ($side =~ /RIGHT/){
		$fseq=$sequences{"${TEid}\_${orient}\_RIGHT"};
		if ($fseq){
		    if ( ($orient =~/\+/) && ($fseq =~/^AAAA/) ){
			my $seqR=$sequences{"${TEid}\_${orient}\_RIGHT"};
			my @baseR=split("",$seqR);
			my $preR="A";
			my $R;
			my $length=0;
			
			foreach (@baseR){
			    $R=$_;
			    if ($preR =~ /A/ && $R !~ /A/){
				print OUT2 "$TEid\tpoly(A)$length\t$orient\tRIGHT\t$seqR\n";
				last;
			    }
			    $preR=$R;
			    $length=$length+1;
			}
		    }
		}
	    }
	    if ($side =~ /LEFT/){
		$fseq=$sequences{"${TEid}\_${orient}\_LEFT"};
		if ($fseq){
		    if ( ($orient =~ /\-/) && ($fseq =~/TTTT$/)){
			my $seqL=$sequences{"${TEid}\_${orient}\_LEFT"};
			my @baseL=split("",$seqL);
			my $preL="A";
			my $L;
			my $length=0;
			
			foreach (reverse(@baseL)){
			    $L=$_;
			    if ($preL =~ /T/ && $L !~ /T/){
				print OUT2 "$TEid\tpoly(A)$length\t$orient\tLEFT\t$seqL\n";
				last;
			    }
			    $preL=$L;
			    $length=$length+1;
			}
		    }
		}
	    }
	}
	close IN1;
	close IN2;
	close OUT2;
    }
}

####################################################

sub segmental_dup_detection {
    
    my $outpath = $_[0];
    my $TE_annot = $_[1];
    my $TE_list = $_[2];
    my $TEflank_seq = $_[3];
    my $refgenome = $_[4];
    my $read_length = $_[5];
    my $min_id = $_[6];

    system("pwd");
    print "blat $TEflank_seq $refgenome ${TEflank_seq}.blast9 -out=blast9\n";
    system("blat $TEflank_seq $refgenome ${TEflank_seq}.blast9 -out=blast9");
    
    my $min_length=($read_length)/2;
    my $max_length=$read_length;

    open (TE, "$TE_annot") or die "ERROR, no input file\n";
    my ($TE,$chr,$start,$end);
    my %TEtable;
    while (<TE>) {
	($TE,$chr,$start,$end)= split();    
	$TEtable{ $TE }=[$chr,($start-100),($end+100)];
    }
    close TE;
    
    open (BLAT, "${TEflank_seq}.blast9") or die "ERROR, no input file\n";
    my ($Q,$S,$idt,$al,$mi,$go,$Qs,$Qe,$Ss,$Se,$ev,$score);
    my %TEhit; 

    #init   
    open (IN, "$TE_list");
    while (<IN>) {
	chomp $_;
	$TEhit{$_}="no_sd";
    }
    close IN;
    
    my ($c1,$s1,$e1);
    my @ar;
    while (<BLAT>) {
	($Q,$S,$idt,$al,$mi,$go,$Qs,$Qe,$Ss,$Se,$ev,$score)= split();    
	if ( exists $TEtable{ $S }) {
	    @ar=$TEtable{ $S };
	    $c1=$ar[0][0];
	    $s1=$ar[0][1];
	    $e1=$ar[0][2];
	    
	    if ( ($c1 eq $Q) && ( ( ( ($Qs >= $s1) && ($Qe <= $e1) ) || ( ($s1 >= $Qs) && ($s1 <= $Qe) ) ) || ( ( ($Qs >= $s1) && ($Qs <= $e1) ) || ( ($e1 >= $Qs) && ($e1 <= $Qe) ) ) ) ) { # if original location 
		print ".";
	    }
	    else{
		
		if (($al >= $min_length) && ($idt >= $min_id)){ # check the id and match length
		    if (( $al >= (0.95*$max_length) )){ # check if one or two TE sides are part of the duplication 
			#print " =====> Two TE sides part of one duplication (no TE insertion): $S $c1 $s1 $e1 /  $_\n";
			$TEhit{ $S }="sd_noTE";
		    }
		    
		    else{
			#check which TE side is part of a duplication
			if ( max($Ss,$Se) <= (1.05*($max_length/2)) ){
			    #print " ==> Left TE side part of a duplication: $S $c1 $s1 $e1 /  $_ \n"; 
			    if ( exists $TEhit{ $S }){
				if( $TEhit{ $S } eq "sd_rigth"){
				    $TEhit{ $S }="sd";  
				}
			    }
			    else{
				$TEhit{ $S }="sd_left";
			    }
			}
			elsif ( min($Ss,$Se) >= (0.95*($max_length/2)) ){
			    #print " ==> Rigth TE side part of a duplication: $S $c1 $s1 $e1 /  $_ \n"; 
			    if ( exists $TEhit{ $S }){
				if( $TEhit{ $S } eq "sd_left"){
				    $TEhit{ $S }="sd";  
				}
			    }
			    else{
				$TEhit{ $S }="sd_rigth";
			    }
			}
		    }
		}
	    }
	}
    }
    close BLAT;
    print "\n";

    # only TEs overlaping a duplication
    open (OUT, ">${TEflank_seq}.blast9_sd");
    foreach my $TE (sort(keys %TEhit)) {
	print OUT $TE,"\t", $TEhit{ $TE }, "\n";
    }   
    close OUT;
}


#########################################################################################################################################################
################################################################ TE DECTECTION APPROACHES  ############################################################## 
#########################################################################################################################################################


#############################################################
#                                                           #
# Detection of the presence of given sequences using BWA    #
#                                                           #
#############################################################



sub PresenceDetectionMultipleStrains {

    my $startdirectory=`pwd`;
    chomp $startdirectory;
    my $inputdir = $_[0]; 
   
    if ($inputdir =~ /\/$/) {
	chop $inputdir;
    }
    my $binreads = $_[1]; 
    #my $binref = $_[2]; 
    my $outputdir = $_[2]; 
    my $maxReadLength = $_[3];
    my $processes = $_[4];
    my $junction = $_[5];
    my $buffer = $_[6];
    my $lim = $_[7]; 
    my $id = $_[8]; 
    my $qual=$_[9];
	my $TE_map=$_[10];
	my $refgenome=$_[11];
    
    opendir (DIR4, "$inputdir");
    print "\nStrain data: $inputdir\n";

    while (defined( my $dir = readdir (DIR4))) {
     
	if ($dir !~ /\./) {
	    print "\nStrain: $dir\n";

	    my @files = <$inputdir/$dir/*>;
	    my $filelist = join (' ', @files);
	    if ($filelist =~ /\.fastq/) {
		print "$outputdir\/$dir\n"; 
		system ("mkdir $outputdir\/$dir") ; 
        if ($binreads == 0) {
		    &generatefiles("$inputdir","$outputdir","$maxReadLength","$dir","$refgenome","$TE_map","$qual","$processes","$buffer");
        }
		&scorecounter("$dir","$outputdir\/$dir\/detection","$maxReadLength","$junction","$buffer","$lim","$id","$qual","$maxReadLength");
	    }	
	}
	else{
	    print "$dir\n";
	}
    }
    open (HDR,"\>file_header");
	# MODIFIED 17/10/17: added the filtered column for reads
    print HDR "strain\tTE\tleft_proxi_average_read_coverage\tleft_distal_average_read_coverage\tleft_match_length\tleft_match_id\tleft_nb_reads\tleft_filtered_nb_reads\tright_proxi_average_read_coverage\tright_distal_average_read_coverage\tright_match_length\tright_match_id\tright_nb_reads\tright_filtered_nb_reads\tresult\n";
    close HDR;

    system ("cat file_header $outputdir\/*\/detection\/results > $outputdir\/results");
    system ("rm -rf file_header");
    closedir DIR4;
    close IN;
}

# Added in release 3: BWA
sub generatefiles {
    my $inputdir = $_[0];
    my $outputdir=$_[1];
    my $maxReadLength=$_[2];
    my $strain=$_[3];
    my $refgenome=$_[4];
    my $TE_map=$_[5];
    my $qual=$_[6];
    my $processes=$_[7];
    my $buffer=$_[8];
	my $strainname = $strain;
    my $nucl = '';
	my @nucl = '';
	my $i = 0;
	my $num_fastq=0;

	print "\nMapping with bwa-mem";
    # Loop for detecting single ends or paired ends
    print("$inputdir\/$strain");
	opendir(fastqDIR, "$inputdir\/$strain") or die $!;
	while (my $fastq = readdir(fastqDIR)) {
		$num_fastq++;
		print "$fastq\n";
	}
    print("$num_fastq");
    # Check if the index is already done
    if (-e "$refgenome.bwt") {
        print("Bwa index already done!\n");
    } else {
        print("Making bwa index\n");
        system("bwa index $refgenome");
    }
    if (-e "$refgenome.fai") {
        print("Fasta index already done!\n");
    } else {
        print("Making FASTA index\n");
        system("samtools faidx $refgenome");
    }
	if ($num_fastq == 3) {
        print("Single ends\n");
		system("bwa mem -t $processes $refgenome $inputdir\/$strain\/$strainname\.fastq > $outputdir\/$strain\/$strainname\_bwa.sam");
	} else {
        print("Paired ends\n");
        system("bwa mem -t $processes $refgenome $inputdir\/$strain\/$strainname\_1.fastq $inputdir\/$strain\/$strainname\_2.fastq > $outputdir\/$strain\/$strainname\_bwa.sam");
	}
	print "\nCreating security bam file\n";
	system("samtools view -Sb $outputdir\/$strain\/$strainname\_bwa.sam > $outputdir\/$strain\/$strainname\_bwa_copy.bam");
	print "\nFilter SAM file for read identity\n";
    #Filtering SAM file by numnber of matches and missmatches in the CIGAR
	print "\nFiltering SAM file by numnber of matches and missmatches in the CIGAR\n";
	open (SAM,"$outputdir\/$strain\/$strainname\_bwa.sam");
	open my $outsam, '>', "$outputdir\/$strain\/$strainname\_bwa_identity.sam" or die "Can't write new file: $!";
	my $identity = $maxReadLength*0.95;
	while (<SAM>) {
		chomp $_;
		my @lines = split (/\n/, $_);
		foreach my $line (@lines) {
			if ($line =~ /^@/) {
				print $outsam $line;
				print $outsam "\n";
			} else {
				my $cigar = (split /\t/, $line)[5];
				my @mscore = ($cigar =~ m/(\d*)M/g);
				my $totalm = 0;
				foreach my $element (@mscore) {
					$totalm = $totalm + $element;
				}
				if ($totalm >= $identity) {
					print $outsam $line;
					print $outsam "\n";
				}
			}
		}
	}
	close $outsam;
	print "\nTransform sam file to bam file\n";
   	system("samtools view -Sb $outputdir\/$strain\/$strainname\_bwa_identity.sam > $outputdir\/$strain\/$strainname\_bwa.bam");
    # Checkpoint
	# my $cmd = "-s $outputdir\/$strain\/$strainname\_bwa.sam";
    #  if    (system($cmd) != 0)
    if (-s "$outputdir\/$strain\/$strainname\_bwa.sam")
	{
    	print "\nDelete SAM file\n";
		system("rm -f $outputdir\/$strain\/$strainname\_bwa.sam");
    }
    else { print "Something is going wrong!\n";}
	# my $cmd_identity = "-s $outputdir\/$strain\/$strainname\_bwa_identity.sam";
    #  if    (system($cmd_identity) != 0)
	if (-s "$outputdir\/$strain\/$strainname\_bwa_identity.sam")
    {
    	print "\nDelete SAM file\n";
		system("rm -f $outputdir\/$strain\/$strainname\_bwa_identity.sam");
    }
    else { print "Something is going wrong!\n";} 
    print "\nSort BAM file\n";
    system("samtools sort $outputdir\/$strain\/$strainname\_bwa.bam > $outputdir\/$strain\/$strainname\_bwa_sorted.bam");
    # Filter bam by mapping quality 30
    system("samtools view -bq $qual $outputdir\/$strain\/$strainname\_bwa_sorted.bam > $outputdir\/$strain\/$strainname\_bwa_MQ_identity.bam"); 
    # Filter bam by readlength
    system("samtools view -h $outputdir\/$strain\/$strainname\_bwa_MQ_identity.bam | awk 'length(\$10) > $maxReadLength-1 || \$1 ~ /^@/' | samtools view -bS - > $outputdir\/$strain\/$strainname\_bwa_MQ.bam");
    # Make the pileup for all genome; at the same time as the previous step it can be generated the fastq file
    system("samtools mpileup -O -aa -f $refgenome $outputdir\/$strain\/$strainname\_bwa_MQ.bam > $outputdir\/$strain\/$strainname\_bwa_MQ_pileups & samtools mpileup -aa -uf $refgenome $outputdir\/$strain\/$strainname\_bwa_MQ.bam | bcftools call -c | vcfutils.pl vcf2fq > $outputdir\/$strain\/$strainname\_cns_MQ.fastq");
    
    if (! -e "$outputdir\/$strain\/chromosomes_up") {
 		mkdir("$outputdir\/$strain\/chromosomes_up");
	}
	if (! -e "$outputdir\/$strain\/ref_fasta") {
 		mkdir("$outputdir\/$strain\/ref_fasta");
	}
    # Create individuals pileups (required TE_pileup.py)
    my $wd = "$outputdir\/$strain";
    my @ref; # Annotation file
    open(my $data, '<', $TE_map) or die "File can't be open: '$TE_map' $!\n";
    while (my $line = <$data>) {
        chomp $line;
        my @fields = split "\t" , $line;
        push @ref, \@fields;
    }
    close($data);
    # Generating Pileup data
    # We are going to read the pileup file line by line to separate each chromosome entries in a different file.
	system("mkdir -p $wd\/chromosomes"); # Create folder in case it was not created previously
    print "[PILEUP] Removing previous separated chromosomas...\n";
    system("rm -rf $wd\/chromosomes/*"); # Remove previous pileup by chomosomes

    print "[PILEUP] Separating pileup by chromosoma...\n";
    system("awk -F\$'\\t' '{filename = \"$wd\/chromosomes\/\" \$1; if (length(\$1) <=2) print > filename }' $outputdir\/$strain\/$strainname\_bwa_MQ_pileups");

    # NOTE: Chromosmoes should be sorted by positon
    system("mkdir -p $wd\/pileups"); # Create folder in case it was not created previously
    print "[PILEUP] Removing previous separated pileups...\n";
    system("rm -rf $wd\/pileups/*"); # Remove previous pileup by chomosomes
    
    # Separate in different lists by chromosome
    my @chromosomes;
    foreach my $line (@ref)
    {
        my $chr = $line->[1];
        if (!grep( /^$chr$/, @chromosomes ) ) {
            push @chromosomes, $chr;
        }
    }
    my $choromosomes_print = join(", ", @chromosomes);
    print "[PILEUP] Chromosomes used: $choromosomes_print \n";
    print "[PILEUP] Generating breakpoints for pileups...\n";
    my %global_breakpoints;
    foreach my $line (@ref)
    {
        my $chromosome = $line->[1];
        my $coord_left = $line->[2] - 1000;
        my $coord_right = $line->[3] - $buffer;
        my $pileup_id = $line->[0];
        my $path_right = "$wd\/pileups\/$pileup_id\_RIGHT.pileup";
        my $path_left = "$wd\/pileups\/$pileup_id\_LEFT.pileup";

        if (!exists($global_breakpoints{$chromosome}{$coord_right})) {
            $global_breakpoints{$chromosome}{$coord_right} = ();
        }

        if (!exists($global_breakpoints{$chromosome}{$coord_left})) {
            $global_breakpoints{$chromosome}{$coord_left} = ();
        }

        push @{$global_breakpoints{$chromosome}{$coord_left}}, $path_left;
        push @{$global_breakpoints{$chromosome}{$coord_right}}, $path_right;
    }

    foreach my $chr (sort keys %global_breakpoints) {
    my $no_line = 0;
    my %breakpoints = %{$global_breakpoints{$chr}};
    my @opened_files = ();
    print "[PILEUP] Generating pileups of chromosome $chr...\n";

    my $file_to_open = "$wd\/chromosomes\/$chr";
    open(my $file, '<', $file_to_open) or die "File can't be open: '$file_to_open' $!\n";

    while (<$file>) {
        $no_line = $no_line + 1;

        if (exists($breakpoints{$no_line})) {
        my @path = @{$breakpoints{$no_line}};
        my @lines_arr;
        my @data = (
            $no_line,
            [@lines_arr],
            [@path]
        );

        push(@opened_files, [@data]);
        }

        foreach my $opened_file (@opened_files) {
        if (defined($opened_file)) {
            push(@{$opened_file->[1]}, $_);
        }
        }

        my @new_opened_files = ();

        foreach my $opened_file (@opened_files) {
        if (defined($opened_file)) {
            if ($no_line == ($opened_file->[0] + (1000 + $buffer))) {
            if (scalar(@{$opened_file->[1]}) != (1001 + $buffer)) {
                my $filename = $opened_file->[2];
                my $size = scalar(@{$opened_file->[1]});
                die("Number of lines for file $filename is not 1061, it's $size");
            }

            foreach my $path (@{$opened_file->[2]}) {
                open my $write_file, '>', $path or die "Cannot open $path: $!";
                foreach (@{$opened_file->[1]}) {
                print $write_file "$_"; # Print each entry in our array to the file
                }
                close($write_file);
            }
            undef $opened_file;
            }
        }
        }
        }
    close($file);
    }

    print "[PILEUP] Generating vertical pileups... \n";
    foreach my $line (@ref) {
        my $te = $line->[0];
        system("cut -f 1,2,3 $wd\/pileups\/$te\_LEFT.pileup > $wd\/pileups\/$te\_LEFT_3.pileup");
        system("cut -f 1,2,3 $wd\/pileups\/$te\_RIGHT.pileup > $wd\/pileups\/$te\_RIGHT_3.pileup");
    }
    # Create the consensus sequence
	open (FQ, "<$outputdir\/$strain\/$strainname\_cns_MQ.fastq") or die "cannot open file";
	while(<FQ>){
		$_ =~ s/\n//g;
		if ($_ =~/^\@(\w{1,2})$/){
			open (CHROM, ">$outputdir\/$strain\/ref_fasta/$1") or die "cannot create";
			open (MAY, ">$outputdir\/$strain\/chromosomes_up/$1") or die "cannot create";
			$i = 0;
		}elsif ($_ =~ /^\+/ || $_ =~ /^\@/){
			$i++;
		}elsif ($i==0){
			@nucl = split ("", $_);
			foreach $nucl (@nucl){
				if ($nucl =~/[ATCGRYMKSWHBVDatcgrymkswhbvd]/){
					my $norm = uc $nucl;
					print MAY $norm;
					my $min = lc $nucl;
					print CHROM $min;
				}elsif ($nucl =~/[Nn]/){
					my $normn = uc $nucl;
					print MAY $normn;
					my $may = uc $nucl;
					print CHROM $may;
				} 
			}
		}
	}
	# Create view files
	print("python scripts\/view_2.py $TE_map $strainname $outputdir\/$strain");

	# Mimic a view file from Tlex using pileup and reference fasta outputs
	my @ref;
	open(my $data, '<', $TE_map) or die "File can't be open: '$TE_map' $!\n";

	my $line;

	while ( $line = <$data>) {
	chomp $line;

	my @fields = split "\t" , $line;
		push @ref, \@fields;
	}

	close($data);

	my $wd = "$outputdir\/$strain";

	my $readlength =  $maxReadLength;
	my $readlength_m = $readlength+$buffer;

	system("rm -rf $wd\/fasta/*");
	system("mkdir -p $wd\/fasta");

	system("rm -rf $wd\/view_files/*");
	system("mkdir -p $wd\/view_files");


	foreach my $line (@ref) {
		my $te = $line->[0];
		my $chrom = $line->[1];
		my $left = $line->[2];
		my $right = $line->[3];
		my @left_boundaries;
		my @right_boundaries;

		print "[FASTA] Procesing $te ...\n";
        print("$te");
        print("$chrom");
        print("$left");

		push @left_boundaries, $left-1000;
		push @left_boundaries, $left+$buffer;
		push @right_boundaries, $right-$buffer;
		push @right_boundaries, $right+1000;

		system("cut -c $left_boundaries[0]-$left_boundaries[1] $wd\/ref_fasta\/$chrom > $wd\/fasta\/$te\_LEFT_.contig_UpperN");
		system("echo \">$te\_LEFT.contig\" > $wd\/fasta\/$te\_LEFT.header");
		system("cat $wd\/fasta\/$te\_LEFT.header $wd\/fasta\/$te\_LEFT_.contig_UpperN > $wd\/fasta\/$te\_LEFT.contig_UpperN");
		system("cut -c $right_boundaries[0]-$right_boundaries[1] $wd\/ref_fasta\/$chrom > $wd\/fasta\/$te\_RIGHT_.contig_UpperN");
		system("echo \">$te\_RIGHT.contig\" > $wd\/fasta\/$te\_RIGHT.header");
		system("cat $wd\/fasta\/$te\_RIGHT.header $wd\/fasta\/$te\_RIGHT_.contig_UpperN > $wd\/fasta\/$te\_RIGHT.contig_UpperN");
		system("cut -c $left_boundaries[0]-$left_boundaries[1] $wd\/chromosomes_up\/$chrom > $wd\/fasta\/$te\_LEFT_up.fasta");
		system("cut -c $right_boundaries[0]-$right_boundaries[1] $wd\/chromosomes_up\/$chrom > $wd\/fasta\/$te\_RIGHT_up.fasta");

		# Parsing fasta consensus sequence:
		system("fold -w 1 $wd\/fasta\/$te\_LEFT_up.fasta > $wd\/fasta\/$te\_LEFT_c.fasta");
		system("fold -w 1 $wd\/fasta\/$te\_RIGHT_up.fasta > $wd\/fasta\/$te\_RIGHT_c.fasta");

		# Mimic MAQ's view file by pasting three first columns from pileup with fasta consensus sequence as a column
		system("paste $wd\/pileups\/$te\_LEFT_3.pileup $wd\/fasta\/$te\_LEFT_c.fasta > $wd\/view_files\/$te\_LEFT_.view");
		system("paste $wd\/pileups\/$te\_RIGHT_3.pileup $wd\/fasta\/$te\_RIGHT_c.fasta > $wd\/view_files\/$te\_RIGHT_.view");

		# 
		system("head -n $buffer $wd\/view_files\/$te\_RIGHT_.view > $wd\/view_files\/$te\_RIGHT.view");
		system("tail -n $buffer $wd\/view_files\/$te\_LEFT_.view > $wd\/view_files\/$te\_LEFT.view");
		system("head -n $readlength_m $wd\/view_files\/$te\_RIGHT_.view | tail -n $readlength  > $wd\/view_files\/$te\_RIGHT_OUT_TE.view");
		system("tail -n $readlength_m $wd\/view_files\/$te\_LEFT_.view | head -n $readlength  > $wd\/view_files\/$te\_LEFT_OUT_TE.view");
	}
	# Delete big pileup file
	#my $cmd_pileup = "-s $outputdir\/$strain\/$strainname\_bwa_MQ_pileups";
    #if    (system($cmd_pileup) != 0)
    if (-s "$outputdir\/$strain\/$strainname\_bwa_MQ_pileups")
	{
    	print "\nDelete SAM file";
		system("rm -f $outputdir\/$strain\/$strainname\_bwa_MQ_pileups");
    }
    else { print "Something is going wrong!";} 
    system ("mkdir $outputdir\/$strain\/detection") ;
    # Save all files in the detection folder
	system("mv $outputdir\/$strain\/pileups\/*T.pileup $outputdir\/$strain\/detection");
    system("mv $outputdir\/$strain\/view_files\/*T.view $outputdir\/$strain\/detection");
    system("mv $outputdir\/$strain\/view_files\/*E.view $outputdir\/$strain\/detection");
    system("mv $outputdir\/$strain\/fasta\/*T.contig_UpperN $outputdir\/$strain\/detection");

}

sub scorecounter {

    my $startdirectory = `pwd`;
    chomp $startdirectory;
    my $strain = $_[0]; 
    my $aligned = $_[1]; 
    my $maxReadLength = $_[2];
    my $junction = $_[3];
    my $buffer = $_[4]; # 60
    my $lim = $_[5]; # Same as Limp = 15
    my $id = $_[6];
    my $qual=$_[7];
    my $readlength=$_[8];
    my @array;
    my $TEname;
    my $countL = 0;
    my $countR = 0; 
    my $matchL= 0;
    my $matchR= 0;
    my $covL = 0;
    my $covR = 0;
    my $contigfile = '';
    my $detection;
    my $ref;
    my $cns;
    my $mismatchL;
    my $mismatchR;
    my $totmatchL;
    my $totmatchR;
    my $idL;
	my $idL_N; # Check if the last 5 nt on the left side of the TE are N
    my $idR;
	my $idR_N; # Check if the last 5 nt on the right side of the TE are N
    my $infoL;
    my $infoR;  
    my $OUTcountL = 0;
    my $OUTcountR = 0; 
    my $OUTmatchL= 0;
    my $OUTmatchR= 0;
    my $OUTmismatchL;
    my $OUTmismatchR;
    my $OUTtotmatchL;
    my $OUTtotmatchR;
    my $OUTidL;
	my $OUTidL_30; # Check the identity of the 30 nt outside the TE on the left side
	my $OUTidL_N; # Check if the last 5 nt on the left side of the junction are N
	my $view_L; # Counter for the N counting
    my $OUTidR;
	my $OUTidR_30; # Check the identity of the 30 nt outside the TE on the right side
	my $OUTidR_N; # Check if the last 5 nt on the right side of the junction are N
	my $view_R; # Counter for the N counting
    my $OUTinfoL;
    my $OUTinfoR; 
    my $matchA;
    my $matchT;
    my $preref;
    my $check_left;
    my $check_right;
    my $polyTractATL="NA";
    my $polyTractATR="NA";
    
    open (OUT,"\>$aligned\/results");
    print "Detect presence of each TE in all the strains ...\n";
    opendir (DIR, "$aligned\/") or die "Invalid directory $aligned\/ specified!\n";

    while (my $file = (readdir (DIR))) { 
 	if ($file =~ /(.*)_LEFT\.contig_UpperN/) {
	    $TEname=$1;
	    my $Nbuffer=($buffer+$maxReadLength);
	  
	    chdir("$aligned");

	    open (INMR,"$TEname\_RIGHT.view");
	    open (INML,"$TEname\_LEFT.view");
	    open (OUTMR,"$TEname\_RIGHT_OUT_TE.view");
	    open (OUTML,"$TEname\_LEFT_OUT_TE.view");
	    open (COVR,"$TEname\_RIGHT.pileup");
	    open (COVL,"$TEname\_LEFT.pileup");
	    open (OUTCOV,">$TEname.tercoverage");
	    
	    my ($threshold,$mint,$maxt,$minb,$maxb);
	    if ( $maxReadLength < (2 * $buffer) ) { # In the case the read length is less than 120 bp
		$threshold=( $maxReadLength / 2);
		$mint=($threshold - ($buffer/2));
		$maxt=($threshold + ($buffer/2));
		$minb=$mint;
		$maxb=$maxt;
		print "###\n\n$TEname\n THRESHOLD = $threshold [ $mint - $maxt ]\n";
	    }
	    else{
		$threshold=$maxReadLength-$buffer; 
		$mint=($threshold - ($buffer/2));
		$maxt=($threshold + ($buffer/2)); 
		$minb=($buffer - ($buffer/2));
		$maxb=($buffer + ($buffer/2));
		print "THRESHOLD = $threshold [ $minb - $maxb ] AND [ $mint - $maxt ]\n";
	    }

	    my $nL = 0;
	    my $covJunctionL = 0;
	    my $covL = 0;
	    my $readposL;
	    my $npLcount=0;
		my $i=0;
	    while (<COVL>) {
			$i++;
		my @lines = split (/\t/, $_);
		my $localcov = $lines[3];
		chomp $localcov;

		if ( $i == $junction+1 ){ 
		    $readposL = $lines[6];
		    chomp $readposL;
		    if ( $readposL ){
			my @read=(split /,/,$readposL);
			foreach (@read){
			    my $oneread=$_;
				# This loop is in charge of selecting which reads has enough bp inside and outside the junction to be considered as valid.
				# It takes into account the maximum read length. 
			    if( ( ( $mint <= $oneread ) && ( $maxt >= $oneread ) ) || ( ( $minb <= $oneread ) && ( $maxb >= $oneread ) ) ){
				print"==> HERE: $lines[0] \t $lines[1] \t $localcov \t $oneread\n";
				$npLcount++;
			    }
			}
		    }
		}
		
        if ($nL >= $junction - $maxReadLength + 1){
		    $covJunctionL = $covJunctionL + $localcov; 
		}
		else{ 
		    $covL = $covL + $localcov; 
		}
		$nL++;
	    }
	    
	    my $nR = 0;
	    my $covJunctionR = 0;
	    my $covR = 0;
	    my $readposR;
	    my $npRcount=0;
		my $r=0;
	    while (<COVR>) {
			$r++;
		my @lines = split (/\t/, $_);
		my $localcov = $lines[3];
		chomp $localcov;
		
		if ( $r == $buffer+1 ){ 
		    $readposR = $lines[6];
		    chomp $readposR;
		    if ( $readposR ){
			my @read=(split /,/,$readposR);
			foreach (@read){
			    my $oneread=$_;
			    if( ( ( $mint <= $oneread ) && ( $maxt >= $oneread ) ) || ( ( $minb <= $oneread ) && ( $maxb >= $oneread ) ) ){
				print"==> HERE: $lines[0] \t $lines[1] \t $localcov \t $oneread\n";
				$npRcount++;
			    }
			}
		    }
		}
		
		if ($nR <= $maxReadLength+$buffer){
		    $covJunctionR = $covJunctionR + $localcov; 
		}
		else{ 
		    $covR = $covR + $localcov; 
		}
		$nR++;
	    }
	    
	    my $avgJunctionL = sprintf("%.2f",($covJunctionL / ($maxReadLength+$buffer)));
	    my $avgJunctionR = sprintf("%.2f",($covJunctionR / ($maxReadLength+$buffer)));
	    my $avgL = sprintf("%.2f",($covL / $junction));
	    my $avgR = sprintf("%.2f",($covR / $junction));
	   
	    print OUTCOV "$TEname\t$avgJunctionL\t$avgJunctionR\t$avgL\t$avgR\t$npLcount\t$npRcount\n";
	    close OUTCOV;
	    
	    # First the left side #####################
	    my $leftside = '';	    
	    {
		local $/=undef;
		open (IN, "$file") or die "file error on file $file!\n";
		binmode IN;
		$leftside = <IN>;
		close IN;
	    }
	    @array = (split />/,$leftside);
	    shift @array;
	    foreach my $flank (@array) {
		if ($flank =~ /^(.*?)\n/) {
		    $contigfile=$1;
		}


		chomp $flank;
		my $mbuffer=($buffer+($lim*2))*-1;
		my $INflank=substr $flank, $mbuffer, $buffer+($lim*2); 
		if (length($INflank) == 90) {
			$countL = ((uc $INflank) =~ tr/N//); # Number of N in the left side of the TE
			$matchL = ($mbuffer+$countL)*-1;
		} else {
			$matchL = 0;
		}
		$mismatchL = 0;
		$totmatchL = 0;
		$idL = 0;
		$idL_N = 0;
		$view_L = 0;
		while (<INML>) {
		    my @line = split();
		    $ref = $line[2];
		    $cns = $line[3];
		    
		    if ( ($ref !~ /N/) && ($cns !~ /N/) ){
			if ( $ref ne $cns){
			    $mismatchL++; # Number of mismatches in the identity
			}
			$totmatchL++;
		    }
            $view_L++;
			if ( $view_L <= 5 && $ref ne $cns ) {
				$idL_N++;
			}
		}
		if ($totmatchL > 0){
		    $idL = sprintf("%.2f",(($totmatchL-$mismatchL)/$totmatchL)*100);
		}
		$infoL = "$idL%($mismatchL/$totmatchL)";
		
		
		my $OUTflank=substr $flank, (($buffer+$maxReadLength)*-1), $maxReadLength;
		$OUTmatchL = $maxReadLength;
		$OUTmismatchL = 0;
		$OUTtotmatchL = 0;
		$OUTidL = 0;
		$OUTidL_30 = 0;
		$OUTidL_N = 0;
		while (<OUTML>) {
		    my @line = split();
		    $ref = $line[2];
		    $cns = $line[3];
		    if ( $ref ne $cns){
			$OUTmismatchL++;
		    }
		    $OUTtotmatchL++;
			if ( $OUTtotmatchL >= $readlength-4 && $ref ne $cns ) {
				$OUTidL_N++;
			}
			if ( $OUTtotmatchL >= $readlength-29 && $ref ne $cns ) {
				$OUTidL_30++;
			}
		}
		if ($OUTtotmatchL > 0){
		    $OUTidL = sprintf("%.2f",(($OUTtotmatchL-$OUTmismatchL)/$OUTtotmatchL)*100);
		}
		$OUTinfoL = "$OUTidL%($OUTmismatchL/$OUTtotmatchL)";
		print "$TEname / LEFT: inside the TE seq = $infoL / outside the TE seq = $OUTinfoL / $OUTidL_30 /$OUTidL_N /$idL_N \n";
	    }
	    
	    # Now the right side #####################	
	    $file =~ s/_LEFT/_RIGHT/;
	    my $rightside = '';
	    {
		local $/=undef;
		open (IN, "$file") or die "file error on file $file!\n";
		binmode IN;
		$rightside = <IN>;
		close IN;
	    }
	    @array = (split />/,$rightside);
	    shift @array;
	    foreach my $flank (@array) {
		if ($flank =~ /^(.*?)\n/) {
		    $contigfile=$1;
		}
		chomp $flank;
		my $mbuffer=($buffer+($lim*2));
		my $INflank=substr $flank,(length $contigfile)+1,$mbuffer;
		$countR = ((uc $INflank) =~ tr/N//);
		$matchR = ($mbuffer-$countR);	
		$mismatchR = 0;
		$totmatchR = 0; 
		$matchA = 0;
		$matchT = 0;
		$idR = 0;
		$idR_N = 0;
		$view_R = 0;
		while (<INMR>) {
		    my @line = split();
		    $ref = $line[2];
		    $cns = $line[3];
		    if ( ($ref !~ /N/) && ($cns !~ /N/) ){
			if ( $ref ne $cns){
			    $mismatchR++;
			}
			$totmatchR++;
		    }
			$view_R++;
			if ( $view_R >= 56 && $ref ne $cns ) {
				$idR_N++;
			}
		}
		if ($totmatchR > 0){
		    $idR = sprintf("%.2f",(($totmatchR-$mismatchR)/$totmatchR)*100);
		}
		$infoR = "$idR%($mismatchR/$totmatchR)";
                
		my $OUTflank=substr $flank, $buffer, $maxReadLength;
		$OUTmatchR = $maxReadLength;
		$OUTmismatchR = 0;
		$OUTtotmatchR = 0;
		$OUTidR = 0;
		$OUTidR_30 = 0;
		$OUTidR_N = 0;
		while (<OUTMR>) {
		    my @line = split();
		    $ref = $line[2];
		    $cns = $line[3];
		    if ( $ref ne $cns){
			$OUTmismatchR++;
		    }
		    $OUTtotmatchR++;
			if ( $OUTtotmatchR <= 5 && $ref ne $cns ) {
				$OUTidR_N++;
			}
			if ( $OUTtotmatchR <= 30 && $ref ne $cns ) {
				$OUTidR_30++;
			}
		}
		if ($OUTtotmatchR > 0){
		    $OUTidR = sprintf("%.2f",(($OUTtotmatchR-$OUTmismatchR)/$OUTtotmatchR)*100);
		}
		$OUTinfoR = "$OUTidR%($OUTmismatchR/$OUTtotmatchR)";
		print "$TEname / RIGTH: inside the TE seq = $infoR / outside the TE seq = $OUTinfoR / $OUTidR_30 /$OUTidR_N /$idR_N \n";
	    }
	    
	    my $pmL = $matchL-($lim*2);
	    my $pmR = $matchR-($lim*2);
		my $npLcount_filtered = '';
		my $npRcount_filtered = '';

		# MODIFICATION 17/10/2017
		# This loop considers the number of N inside the 60 TE bp and 30 bp outside the junction.
		# There should be a minimum of 45 non-N bp, so it ensures a minimum of 15 bp inside the TE (lim parameter).
	    if ( (($pmL > $lim && $idL >= $id && $npLcount != 0) || ($pmR > $lim && $idR >= $id && $npRcount != 0))) {
		$detection="present";
		$npLcount_filtered = $npLcount;
		$npRcount_filtered = $npRcount;
	    }
	    else{
            # MODIFICATION DECEMBER 2017
			if ($idL == 0 && $pmL >= -10 && $OUTidL_30 <= 11 || $idR == 0 && $pmR >= -10 && $OUTidR_30 <= 11 ) {
				$detection="absent";
				$npLcount_filtered = "0";
				$npRcount_filtered = "0";
				}
            # MODIFICATION 22/12/2017
			elsif ( $idL_N + $OUTidL_N < 10 && $pmL > -5 ||  $idR_N + $OUTidR_N < 10 && $pmR > -5) {
				$detection="absent";
				$npLcount_filtered = "0";
				$npRcount_filtered = "0";
			}
            ## ORIGINAL ABSENT CONDITION
			# elsif ($pmL > -5 || $pmR > -5) {
			# 	$detection="absent";
			# 	$npLcount_filtered = "0";
			# 	$npRcount_filtered = "0";
			# }
			else {
				$detection="no_data";
				$npLcount_filtered = "0";
				$npRcount_filtered = "0";
			}
	    }
	    my $read_number_presence=(($npRcount+$npLcount)/2);

	    print "$strain\t$TEname\t$avgJunctionL\t$avgL\t$pmL\t$infoL\t$npLcount\t$npLcount_filtered\t$avgJunctionR\t$avgR\t$pmR\t$infoR\t$npRcount\t$npRcount_filtered\t$detection\n";
	    print OUT "$strain\t$TEname\t$avgJunctionL\t$avgL\t$pmL\t$infoL\t$npLcount\t$npLcount_filtered\t$avgJunctionR\t$avgR\t$pmR\t$infoR\t$npRcount\t$npRcount_filtered\t$detection\n";
	}
    }
    close OUT;
    chdir("$startdirectory");
}


#############################################################
#                                                           #
# Detection of the absence of given sequences using SHRIMP  #
#                                                           #
#############################################################


sub AbsenceDetectionMultipleStrains{
    my $startdirectory=`pwd`; 
    chomp $startdirectory;
    my $strains=$_[0];
    my $output=$_[1]; 
    my $TE_list=$_[2];
    my $binref=$_[3];
    my $binreads=$_[4];
    my $lim=$_[5];
    my $var=$_[6];
    my $flank=$_[7];
    my $PE=$_[8];
    my $minid=$_[9];
    my @folderlist;
    opendir (DIR, $strains);
    STRAINSLOOP: while (defined (my $strainsdir = readdir (DIR))) {
      if ($strainsdir !~ /\./) {
	  push (@folderlist, $strainsdir);
      }
  }  
    foreach my $direc (@folderlist) { 
	print "$direc\n";
	&LaunchSHRIMP("$strains\/$direc","$output\/$direc","$TE_list","$binref","$binreads","$lim","$var","$flank","$PE","$minid");
    }
    chdir ("$startdirectory");
    system("echo Strain\tTEname\tDetection\tSelectedReadNb \> $output\/results");
    system("cat $output\/*\/detection\/results \>> $output\/results");
} 

sub LaunchSHRIMP {
    my $startdirectory=`pwd`;
    print "$startdirectory";
    chomp $startdirectory;
    my $input=$_[0]; 
    my $output=$_[1]; 
    my $TE_list=$_[2];
    my $binref=$_[3]; 
    my $binreads=$_[4];
    my $lim=$_[5];
    my $var=$_[6];
    my $flank=$_[7];
    my $PE=$_[8];
    my $minid=$_[9];
    print "$minid\n\n\n\n\n";
    my @list=split(/\//,$input);
    my $len=@list;
    my $strain=$list[$len-1];

    my %scoretable;
    my %pgenometable;
    open (IN, "$TE_list") or die "Invalid TE_list specified!\n";
    while (<IN>) {
	chomp $_;
	$scoretable{$_}=0;
	$pgenometable{$_}=0.00;
    }
    close IN;
   
    system ("mkdir $output");
    system ("mkdir $output\/mapping");
    
    if ($binreads == 0) {
	system("mkdir $output\/fasta") ;
	opendir (DIR, $input) or die "Failed to open $input!\n";
	print "Convert the read file from fastq to fasta format ...\n";
	while (defined (my $dir = readdir(DIR))) { 
	    if ($dir =~ /(.*)\.(fastq|fq)/) {
		print " - $dir\n";
		&fastq2fasta("$input\/$dir","$output\/fasta\/$1","500000");
	    }
	}
    }
    else {
	print "Read files already converted in fasta files !!!\n";
    }
    my $end_formatting = localtime();
    print "* End of the format step (SHRIMP): $end_formatting\n";
    my $refseq="$binref";

    my %hrefseq=();
    my $hseq ="";
    my $linecount=0;
    my (@hTE, $hTEn, $fa);
    open(REF, $refseq);
    while (<REF>) {
      if ($_ =~ /^>/){
	  chomp;
	  @hTE = split(/>/,$_);
	  $hTEn = $hTE[1];
	  $linecount=1;
      }
      else{   
	  if ($linecount==1) {
	      $fa= $_;
	      chomp $fa;
	      $linecount=2;
	      next;
	      
	  }
	  if ($linecount==2) {
	      $fa= $fa."".$_;
	      chomp $fa;
	      $hrefseq{$hTEn}=$fa;
	      $linecount=0;
	  }
      }
    }
    close REF;

    ##############################################################################################
    print "Mapping step ...\n";
  
    my @readR; 
    my %arpe;  
    my ($editstring, $peid, $readid, $TEid, $strand, $qstart, $qend , $seq);
    my ($editstring_next, $peid_next, $readid_next, $TEid_next, $strand_next, $qstart_next, $qend_next, $seq_next);
    my @peida;
    my $peidname;
    my $merge="no";
    my @O=split(/\//,$output);
    my $outn=$O[2];
    print "$outn\n\n";
    &list_readfiles("$output\/fasta\/","$outn","$PE");
    
    # Paired-end reads
    if( $PE =~/yes/){
	system("cp $startdirectory\/$output\/fasta\/list_${outn}_PE list_${outn}_PE");
	print " \n#### PAIRED-END READS ####\n";
    	open (PE, "list_${outn}_PE" ) or die "Failed to open file list_${outn}_PE\n";
	while (<PE>){
	    my ($fn,$fc) = split();
	    system ("gmapper-ls -p opp-in -P -R -i -20 -g -40 -q -40 -e -1 -f -1 -1 $startdirectory\/$output\/fasta\/${fn}1_${fc}\.fa -2 $startdirectory\/$output\/fasta\/${fn}2_${fc}\.fa $refseq \> $startdirectory\/$output\/fasta\/${fn}_${fc}_fa_gmapper");
	}
	close PE;
	
	system ("mv $output\/fasta\/*fa_gmapper $output\/mapping\/"); 
	system ("cat $output\/mapping\/*fa_gmapper \> $output\/mapping\/results_gmapper2merge"); 
	print "PE Mapping done\n";

	open (INPE, "$output\/mapping\/results_gmapper2merge") or die "Failed to open file $output\/mapping\/results_gmapper2merge\n";
	while (<INPE>) {   
	    if ($_ =~ /^>/){
		my @tmp=split();
		my $ll = scalar(@tmp);
		if ( $ll > 1){
		    ($editstring, $readid, $TEid, $strand, $qstart, $qend, $seq) = (split)[9,0,1,2,3,4,10];
		    @peida = split (/\//,$readid);
		    $peidname = $peida[0];
		    push(@readR, [$peidname, $readid, $TEid, $strand, $qstart, $qend, $editstring, $seq]);  
		} 
	    }
	} 

	print "# SORT...\n";
	my @read_info = sort {$a->[0] cmp $b->[0]} @readR;
	$peid = $read_info[0][0]; 
	$TEid = $read_info[0][2];
	$readid = $read_info[0][1];
	$strand = $read_info[0][3]; 
	$qstart = $read_info[0][4];
	$qend = $read_info[0][5];
	$editstring = $read_info[0][6];
	$seq = $read_info[0][7];
 
	$peid_next = $read_info[1][0];
	$TEid_next = $read_info[1][2];
	$readid_next = $read_info[1][1];
	$strand_next = $read_info[1][3]; 
	$qstart_next = $read_info[1][4];
	$qend_next = $read_info[1][5]; 
	$editstring_next = $read_info[1][6];
	$seq_next = $read_info[1][7];

     
	print "# MERGE...\n";
	for my $i ( 1 .. $#read_info ) {
	    if ($i < $#read_info){
		my $next = $i+1;
		$peid = $read_info[$i][0]; 
		my $speid=substr($peid,1, length($peid));
		$TEid = $read_info[$i][2];
		$peid_next = $read_info[$next][0];
		$TEid_next = $read_info[$next][2];
		
		if ( ( $peid =~ $peid_next ) && ( $TEid =~ $TEid_next) ){
		    $readid = $read_info[$i][1];
		    $readid_next = $read_info[$next][1];
		    if ( (( $readid =~ /\/1$/) && ($readid_next =~ /\/2$/)) || (($readid =~ /\/2$/) && ($readid_next =~ /\/1$/)) ){
			$strand = $read_info[$i][3]; 
			$qstart = $read_info[$i][4];
			$qend = $read_info[$i][5];
			$seq = $read_info[$i][7];
			$strand_next = $read_info[$next][3]; 
			$qstart_next = $read_info[$next][4];
			$qend_next = $read_info[$next][5]; 
			$seq_next = $read_info[$next][7];
			open (OUTR1,'>tmp1.fa'); 
			print OUTR1 "$readid\_$TEid\n$seq\n";
			close OUTR1;
			open (OUTR2,'>tmp2.fa');
			print OUTR2 "$readid_next\_$TEid_next\n$seq_next\n";
			close OUTR2;
			my @coord = ($qstart,$qstart_next,$qend,$qend_next);
			my @coord_sort = sort {$a <=> $b} @coord;
			my $mseq;
			my $over = $coord_sort[2]-$coord_sort[1];
			
			if ($over >= $var){
			    if (((($qend >= $qstart_next) && ($qend <= $qend_next)) || ( ($qstart_next >= $qstart) && ($qstart_next <= $qend))) || ((($qstart >= $qstart_next) && ($qstart <= $qend_next)) || ( ($qend_next >= $qstart) && ($qend_next <= $qend)))){
				$merge="yes";
			      	
				if ( ( $strand_next =~ /\+/ ) && ( $strand =~ /\-/) ){	    
				    delete $read_info[$i]; 
				    #print "Overlap: [$qstart - $qend] [$qstart_next - $qend_next]\n";
				    $read_info[$next][4] = $coord_sort[0];
				    $read_info[$next][5] = $coord_sort[3];
				    
				    system("fastx_reverse_complement -i tmp1.fa > tmp1.fa_rc");	
				    system("merger -asequence tmp1.fa_rc -bsequence tmp2.fa -outfile merger.align -outseq merger.fa_${TEid} &> stdoutput");
				    open (IN, "merger.fa_${TEid}");
				    while (<IN>){
					chomp;
					if ($_ !~ /^>/){
					    $mseq .="$_";
					}
				    }
				    $read_info[$next][7] = $mseq;
				    $read_info[$next][1] = $peid_next."/merge"; 
				}
				
				if ( ( $strand =~ /\+/ ) && ( $strand_next =~ /\-/) ){ 
				    delete $read_info[$i]; 
				    #print "Overlap: [$qstart - $qend] [$qstart_next - $qend_next]\n";
				    
				    $read_info[$next][4] = $coord_sort[0];
				    $read_info[$next][5] = $coord_sort[3]; 
				    
				    system("fastx_reverse_complement -i tmp2.fa > tmp2.fa_rc");
				    system("merger -asequence tmp1.fa -bsequence tmp2.fa_rc -outfile merger.align -outseq merger.fa_${TEid} &> stdoutput");
				    open (IN, "merger.fa_${TEid}");
				    while (<IN>){
					chomp;
					if ($_ !~ /^>/){
					    $mseq .="$_";
					}
				    }
				    $read_info[$next][7] = $mseq;
				    $read_info[$next][1] = $peid_next."/merge"; 
				}
			    
				
				# Extract the ref seq
				my $TEfa = $TEid.".fa";
				open(OUTREF, "> $TEfa");
				print OUTREF ">".$TEid."\n".$hrefseq{$TEid}."\n";
				close OUTREF;
				system("blat merger.fa_${TEid} $TEfa  merger.fa_${TEid}blat -out=axt &> stdoutput");
				if( -s "merger.fa_${TEid}blat"){
				    open(IB, "merger.fa_${TEid}blat");
				    my @lines = <IB>;
				    my @l = split(/ /,$lines[0]);
				    $read_info[$next][4] = $l[5];
				    $read_info[$next][5] = $l[6]; 
				    chomp $lines[2];
				    $read_info[$next][6] = $lines[2];
				    chomp $lines[1];
				    $read_info[$next][7] = $lines[1];
				}
			    }
			}
			else{ # no overlap but some of them can be linked 
			    $read_info[$i][1] = $peid;
			    $read_info[$next][1] = $peid_next;
			}
		    }
		}
	    }
	}

	system("rm -f tmp1* tmp2* *merger.");
	open (OUTM, ">results_gmapper");
	for my $i ( 1 .. $#read_info ) { 
	    if($read_info[$i][1]){
		$TEid = $read_info[$i][2];
		$readid = $read_info[$i][1];
		$strand = $read_info[$i][3]; 
		$qstart = $read_info[$i][4];
		$qend = $read_info[$i][5];
		$editstring = $read_info[$i][6];
		$seq = $read_info[$i][7];
		print OUTM "$readid\t$TEid\t$strand\t$qstart\t$qend\t1\tX\tX\t-1\t$editstring\t$seq\n";
	    }
	}
	close OUTM;
	system("mv results_gmapper $output\/mapping\/results_gmapper");
    }
    
    else{
	# Single reads
	print " \n#### SINGLE READS ####\n";
	system("pwd");
	system("cp $startdirectory\/$output\/fasta\/list_${outn}_S list_${outn}_S"); 
	open (S, "list_${outn}_S" ) or die "Failed to open file list_${outn}_S\n";
	### HERE ###
	while (<S>){
	    print $_;
	    my ($fn,$fc) = split();
	    system ("gmapper-ls -P -R -i -20 -g -40 -q -40 -e -1 -f -1 $startdirectory\/$output\/fasta\/${fn}_${fc}\.fa $refseq \> $startdirectory\/$output\/fasta\/${fn}_${fc}_fa_gmapper");
	}
	close S;
	system("rm list_${outn}_S");
	
	system ("mv $output\/fasta\/*fa_gmapper $output\/mapping\/");
	system ("cat $output\/mapping\/*fa_gmapper \> $output\/mapping\/results_gmapper");	    
	
	print "Mapping done\n";
    }
    
    ###########################################################
    # Last step of the Absence detection:                     #
    # ----------------------------------                      #
    # - Pull out all the reads matching across the boundaries #
    # - Count them and output them in aligned fasta format    #
    ###########################################################
    
    print "detection of reads/PEs overlapping the two TE sides...\n";
    system("mkdir $output\/detection");
    system ("cp $output\/mapping\/results_gmapper $output\/detection\/results_gmapper");	   
    chdir ("$output\/detection");
  
   
    open (IN, 'results_gmapper');
    open (OUT, '>results_gmapper_overjunction');
    my $minlim=$flank-$var; 
    my $maxlim=$flank+$var;
    print "VAR = $var\n";
    while (<IN>) {
	if ($_ =~ /^>/) {
	    my @line = split (/\t/, $_);
	    if ( $line[0] =~ /merge$/) {
		if (($line[3] < $minlim && $line[4] > $maxlim) || ($line[3] > $maxlim && $line[4] < $minlim)) {
		    my $editstring = $line[9];
		    if( $editstring !~ m/(\d+)/g){
			print OUT $_;
		    }
		    else{
			print "This is an edit string $_";# test avec une fake seq et un grand gap!!!!
		    }
		}
	    }
	    else{ 
		if (($line[3] < $minlim && $line[4] > $maxlim) || ($line[3] > $maxlim && $line[4] < $minlim)) {
		    my $editstring = $line[9];
		    my @a1 = [];
		    my @a2 = [];
		    my @a3 = [];
		    my @a4 = [];
		    my @id3 = [];
		    my $sum=0;
		    my $id=0;
		    my $idL=0;
		    my $idR=0;
		    my $mlength=0;
		    my $termatchL=0;
		    my $termatchR=0;

		    #print "FULL: $editstring\n";	    
		    @a1 = split("---",$editstring);
		    #print scalar(@a1);
		    #print "\n";
		    if(scalar(@a1)==1){
			$mlength=0;
			@a2 = ($a1[0] =~ m/(\d+)/g);
			@a4 = ($a1[0] =~ m/[ACGT]/g);
			$sum=0;
			foreach (@a2) {
			    $sum= $_ + $sum;
			}
			$mlength = scalar(@a4);
			$termatchL = $sum;
			$termatchR = $termatchL;
			$idL = ($sum/($mlength+$sum))*100;
			$idR = $idL;
			print "Editstring: $editstring \t $sum \t $mlength \t $idL\n";
		    }
		    else{
			foreach (@a1) {
			    if($_){
				#print "editstring: $_\n";
				$mlength=0;
				@a2 = ($_ =~ m/(\d+)/g); # nb of matches 
				@a4 = ($_ =~ m/[ACGT]/g);# nb of mismatches
				$sum=0;
				foreach (@a2) {
				    $sum= $_ + $sum;
				}
				$mlength = scalar(@a4);
				if ($mlength ==0){
				    if ($sum ==0){
					$id = 0;
				    }
				    else{
					$id = 100;
				    }
				}
				else{
				    $id = ($sum/($mlength+$sum))*100;
				}
				#print "Editstring: $editstring \t $sum \t $mlength \t $id\n";
				push(@a3, "$sum");
				push(@id3, "$id");
			    }
			}
			$termatchL = $a3[1];
			$termatchR = $a3[scalar(@a3)-1];
			$idL = $id3[1];
			$idR = $id3[scalar(@id3)-1];
		    }
		    
		    if (( $termatchL>=$var ) && ($termatchR>=$var) && ($idL>=$minid) && ($idR>=$minid)) { 
			#print "$editstring read length[$termatchL - $termatchR]  / read id [$idL-$idR]\n";
			print OUT $_;
		    }
		}
	    }
	}
    }
    close OUT;
    close IN;
    
    open (IN, 'results_gmapper_overjunction');
    my @readl;    
    while (<IN>) {
	my ($editstring, $readid, $TEid, $qstart, $qend , $readseq) = (split)[9,0,1,3,4,10];
	push(@readl, [$editstring, $readid, $TEid, $qstart, $qend, $readseq]);
    } 
    
    my @read_list = sort {$a->[0] cmp $b->[0]} @readl;
    
    my $prec_editstring = $read_list[0][0]; 
    my $prec_readid = $read_list[0][1]; 
    my $prec_TEid = $read_list[0][2];
    my $prec_qstart = $read_list[0][3];
    my $prec_qend = $read_list[0][4];
    my $i=1;
    
    while ($i<=$#read_list ) {
	if ($prec_editstring =~ /$read_list[$i][0]/ ){
	    if ($prec_TEid =~ /$read_list[$i][2]/){
		if ($prec_qstart == $read_list[$i][3] &&  $prec_qend == $read_list[$i][4]){ 
		    splice(@read_list,$i,1); 
		    $i--;
		}
	    }
	}
	$prec_editstring = $read_list[$i][0]; 
	$prec_readid = $read_list[$i][1];
	$prec_TEid = $read_list[$i][2];
	$prec_qstart = $read_list[$i][3];
	$prec_qend = $read_list[$i][4];
	$i++;
    }
    
    open (OUT, '>results_gmapper_overjunction_nodup');
    open (OUT2,'>results_gmapper_overjunction_nodup_format');
    open (OUT3,'>results_gmapper_overjunction_nodup_format.fa');
    for my $i ( 1 .. $#read_list ) {
	$read_list[$i][1] =~ s/^>//;
	print OUT "$read_list[$i][1]\t$read_list[$i][2]\n";
	$read_list[$i][1] =~ s/^.*?://;
	my $nstr = "$read_list[$i][2]\/$read_list[$i][1]";
	print OUT2 "$nstr\n";
	print OUT3 ">$nstr\n$read_list[$i][5]\n";
    }
    close OUT;
    close OUT2;
    close OUT3; 
    close IN;
    
    system("sort results_gmapper_overjunction_nodup_format | uniq > results_gmapper_overjunction_nodup_format_uniq");
    system("mv results_gmapper_overjunction_nodup_format_uniq results_gmapper_overjunction_nodup_format");

    open (OUT2, '>results');
    
    if( -z "results_gmapper_overjunction_nodup_format.fa"){
	print "file empty ************************\n";
    }
    else{
	print "#### RM ####\n"; 
	system("sed s/'\-'/''/g results_gmapper_overjunction_nodup_format.fa > results_gmapper_overjunction_nodup_format_nogap.fa");
	# Mask the single repeats of the reads
	system("RepeatMasker -noint results_gmapper_overjunction_nodup_format_nogap.fa"); 
	
	# Select only the reads with non-masked terminal regions (of more than 5 bp)
	open (IN,'results_gmapper_overjunction_nodup_format_nogap.fa.out');
	open (OUT,'>results_gmapper_overjunction_nodup_format2drop');
	my $n = 0 ;
	while (<IN>) { 
	    my $line = $_ ;
	    if ($n > 2 ) {
		$line =~ s/\(//g;
		$line =~ s/\)//g;
		my @lines = split (/ +/, $line);
		my $id = $lines[5];
		my $start = $lines[6];
		my $end = $lines[7];
		my $max = ($lines[7]+$lines[8]) - $lim;
		if ($start < $lim || $end > $max) { 
		    print OUT "$id\n";
		   
		}
	    }
	    $n++;
	}
	close OUT;
	close IN;
	
	# Diff of the all list of the reads in the cross file and the reads to eliminate
	system("sort results_gmapper_overjunction_nodup_format > results_gmapper_overjunction_nodup_format_sort");
	system("sort results_gmapper_overjunction_nodup_format2drop > results_gmapper_overjunction_nodup_format2drop_sort");
	system("diff results_gmapper_overjunction_nodup_format_sort results_gmapper_overjunction_nodup_format2drop_sort > results_gmapper_overjunction_nodup_format2select");
	
	open (IN, 'results_gmapper_overjunction_nodup_format2select');
	open (OUT, '>results_gmapper_overjunction_nodup_format2select_format');
	#open (IN, 'results_gmapper_overjunction_nodup_format2select_short');
	#open (OUT, '>results_gmapper_overjunction_nodup_format2select_format_short');
	# to remove ### open (IN2,'results_gmapper_overjunction_nodup_format_nogap.fa.masked');
	
	while (<IN>) {
	    if ($_ =~ /^</) {
		my $lines = $_;
		my @line = split (/< /,$lines);
		my $linf = $line[1]; 
		my @line1 = split (/\//, $linf);
		my $TEname = $line1[0];
		my $Readname = $line1[1];
		print " $lines -> $linf -> $TEname / $Readname \n";
		my $Merge = "";
                if ($Readname =~ /merge/){
		    $Merge = "read_merged";
		}
		chomp $TEname;
		chomp $Readname;
		print OUT "$TEname\t$Readname\t$Merge\n";
		$scoretable{$TEname}++;
	    }
	}
    }
    print Dumper %scoretable;
    while ( my ($key, $value) = each(%scoretable) ) {
	if ($value == 0){
	    system("grep $key results_gmapper > tmp_$key "); ### slow step ###
	    
	    if (-z "tmp_$key"){
		print OUT2 "$strain\t$key\tno_data\t$value\n";
	    }
	    else{
		print OUT2 "$strain\t$key\tpresent\t$value\n";
	    }
	}
	else{
	    print OUT2 "$strain\t$key\tabsent\t$value\n";
	}
    }
    close OUT;
    close IN;
    close OUT2;
    #close IN2;
    chdir("$startdirectory");
}



#########################################################################################################################################################
##################################################### Combine Presence and Absence detection results ####################################################
#########################################################################################################################################################

sub FinalResults {   
    my $output = $_[0];
    my $TE_list = $_[1]; 
    my $TE_annot= $_[2];
    my $flank= $_[3];  

    print "$output";
    chdir("$output");
    
    open (OUT,">Tresults");
    
	# MODIFIED 17/10/17: including filtered reads after absence detection; following modifications are due to the same issue
    print OUT "strain\tTE\tabsence_detection\tpresence_detection\tcombination\tread_number_absence\tleft_match_length\tleft_match_id\tpolyAT_left\tleft_coverage\tleft_repeat\tleft_read_number_presence\tleft_filtered_read_number_presence\tright_match_length\tright_match_id\tpolyAT_right\tright_coverage\tright_repeat\tright_read_number_presence\tright_filtered_read_number_presence\n";	
    

    open (RES1, 'Tresults_absence');
    my %ABS = ();
    my %NBR = ();
    my %Pgenome = ();
    my %Pchance = ();
    my $na = 0;
    my ($straina, $TEa, $resa, $NR);
    
    while (<RES1>) {
	if ($na > 0 ){
	    ($straina, $TEa, $resa, $NR)=split(/\t/,$_);
	    $ABS{ $straina }{ $TEa } = $resa;
	    $NBR{ $straina }{ $TEa } = $NR;
	}
	$na++;
    }
    close RES1;

    open (RES2, 'Tresults_presence'); 
    my %PRE = ();
    my %CR = ();
    my %CL = ();
    my %MR = ();
    my %ML = ();  
    my %IR = ();
    my %IL = (); 
    my %CNR = ();
    my %CNL = ();
    my %RR;
    my %RL; 
    my %TAL;
    my %TAR;
	my %CNRF = ();
	my %CNLF = ();
 
    # MODIFIED 17/10/17
    my ($strainp, $TEp, $coverageL, $covflankL, $matchL, $idL, $coverageR, $covflankR, $matchR, $idR, $nbreadL, $nbreadLF, $nbreadR, $nbreadRF, $resp);
    while (<RES2>) {
	chomp;

	if( $_ !~ /^strain/){
		# MODIFIED 17/10/17
	    ($strainp, $TEp, $coverageL, $covflankL, $matchL, $idL, $nbreadL, $nbreadLF, $coverageR, $covflankR, $matchR, $idR, $nbreadR, $nbreadRF, $resp)=split(/\t/,$_);
	    $PRE{ $strainp }{ $TEp } = $resp;
	    $CR{ $strainp }{ $TEp } = $coverageR; 
	    $CL{ $strainp }{ $TEp } = $coverageL;
	    $MR{ $strainp }{ $TEp } = $matchR; 
	    $ML{ $strainp }{ $TEp } = $matchL; 
	    $IR{ $strainp }{ $TEp } = $idR; 
	    $IL{ $strainp }{ $TEp } = $idL; 
	    $CNR{ $strainp }{ $TEp } = $nbreadR; 
	    $CNL{ $strainp }{ $TEp } = $nbreadL;
	    $CNRF{ $strainp }{ $TEp } = $nbreadRF; 
	    $CNLF{ $strainp }{ $TEp } = $nbreadLF;

	}
    }
    close RES2;

    open (IN, "../$TE_list");
    my @TEs;
    while (<IN>) {
	chomp;
	push (@TEs,$_);
    }
    close IN;
    
    print "INIT\n";
    open (IN, "../$TE_list");
    while (<IN>) {
	print $_;
	chomp $_;
	$RR{$_}="no_repeat";
	$RL{$_}="no_repeat";
	$TAL{"left"}{$_}="no_polyAT";
	$TAR{"left"}{$_}="no_polyAT";
	$TAL{"right"}{$_}="no_polyAT";
	$TAR{"right"}{$_}="no_polyAT";
    }
    close IN;

    my $nl=0;
    open (RES3, "Tanalysis\/Tflank_checking_${flank}.fasta.out"); 
    while (<RES3>) {
	if ($nl>2){
	    my ($TE, $repeat) = (split)[4,9];
	    if ($TE =~ /RIGHT$/) { 
		my @SplitTEname = split (/\_/, $TE);
		my $nTE = $SplitTEname[0];
		$RR{ $nTE } = $repeat;
	    }
	    if ($TE =~ /LEFT$/) { 
		my @SplitTEname = split (/\_/, $TE);
		my $nTE = $SplitTEname[0];
		$RL{ $nTE } = $repeat;
	    } 
	}
	$nl++;
    }
    close RES3;

    open (RES4, "Tanalysis\/Tpoly_${flank}.fasta.polyAT"); 
    while (<RES4>) {
	my ($TE, $side, $tail) = (split)[0,3,1];
	if( $side =~ /LEFT$/){
	    $TAL{"left"}{ $TE } = $tail;
	}
	if( $side =~ /RIGHT$/){
	    $TAR{"right"}{ $TE } = $tail;
	}
    }    
    close RES4;

    ### HERE ### add duplication info open (OUT, ">$outpath\/${TE_annot}_sd");
    #open (RES5, "Tanalysis\/${TE_annot}_sd"); #    my $TE_annot = $_[1];
    #while (<RES5>) {
    #}
    ############################################ HERE ########################



    my $conclusion;
    for my $S ( sort keys %PRE ) {
        for my $TE ( keys %{$PRE{ $S }} ) {
	    my $resab = $ABS{ $S }{ $TE };
	    my $respr = $PRE{ $S }{ $TE };
	    my $nbreada = $NBR{ $S }{ $TE };
	    chomp $nbreada;
	    my $covR = $CR{ $S }{ $TE };
	    my $covL = $CL{ $S }{ $TE };
	    my $repeatR = $RR{ $TE };
	    my $repeatL = $RL{ $TE };
	    $matchR = $MR{ $S }{ $TE }; 
	    $matchL = $ML{ $S }{ $TE }; 
	    my $nbreadR = $CNR{ $S }{ $TE }; 
	    my $nbreadL = $CNL{ $S }{ $TE };
	    my $nbreadRF = $CNRF{ $S }{ $TE }; # MODIFIED 17/10/17
	    my $nbreadLF = $CNLF{ $S }{ $TE }; # MODIFIED 17/10/17
	    my $infoR = $IR{ $S }{ $TE }; 
	    my $infoL = $IL{ $S }{ $TE }; 
	    my $PolyTAL = $TAL{ "left" }{ $TE }; 
	    my $PolyTAR = $TAR{ "right" }{ $TE }; 
	    my @ar = split(/%/,$infoR); 
	    $idR = $ar[0];
	    my $mR = $ar[1];
	    @ar = split(/%/,$infoL);
	    $idL = $ar[0];
	    my $mL = $ar[1];
	    
	    
            ################################### 
	    # The different classes:          #
	    # ---------------------           #
            #  absent                         #
            #  absent/polymorphic             #
            #  present                        #
	    #  present/polymorphic            #
	    #  polymorphic                    #
            #  no_data                        #
	    ###################################
	   
       ### Modification in release 3: elimination of absent/poylmorphic and present/polymorphic
	    if ($respr eq 'present'){ 
		if ($resab eq 'present'){ 
		    $conclusion = "present"; 
		}
		elsif($resab eq 'no_data'){
		    $conclusion = "no_data"; 
		}
		else{
		    $conclusion = "polymorphic";
		}
	    }
	    elsif ($respr eq 'absent'){ 
		if ($resab eq 'absent'){
		    $conclusion = "absent";
		}
		else{
		    $conclusion = "no_data";
		}
	    }
	    elsif ($respr eq 'no_data'){
		if ($resab eq 'absent'){
		    $conclusion = "no_data"; 
		}
		else{
		    $conclusion = "no_data"; 
		}
	    }
		# MODIFIED 17/10/17
	    print OUT "$S\t$TE\t$resab\t$respr\t$conclusion\t$nbreada\t$matchL\t$idL$mL\t$PolyTAL\t$covL\t$repeatL\t$nbreadL\t$nbreadLF\t$matchR\t$idR$mR\t$PolyTAR\t$covR\t$repeatR\t$nbreadR\t$nbreadRF\n"; 
	}
    }
    close RES1;
    close RES2; 
    close RES3;
    close RES4; 
    #close RES5;
    close OUT;
}

#########################################################################################################################################################
################################################################## TE FREQUENCY ESTIMATES ############################################################### 
#########################################################################################################################################################

sub FreqEstimate {

    my $output = $_[0];
    my $pooled = $_[1];
    my $maxreads = $_[2];
    my $minreads = $_[3];
    my $minpop = $_[4];
    chdir("$output");
    
    my $nb=0;
    my $strain_name;
    my $oldstrain_name="";
    open(IN, "Tresults");
    while(<IN>){
	$strain_name = (split)[0];
	print "$strain_name\n";

    if($strain_name !~ $oldstrain_name){
	    $nb++;
	}
	$oldstrain_name=$strain_name;
    }
    my $totstrain=$nb;
    close IN;
   

    if ($pooled){

    # Release 3: modification in the pooled frequency estimation

	print "frequency estimates using pooled data...\n";  
	open (OUTP,">Tfreq_pooled"); 
	print OUTP "TE\tnumber_reads_absence\tnumber_reads_presenceL\tnumber_reads_presenceR\tTEfrequency\n";

	my %RA;
	my %RPL;
	my %RPR;
	my ($strain, $TE, $abs, $preL, $preR);
	my $n = 0;
	open (IN,"Tresults");
	while (<IN>) {
	    if ($n > 0){
		my @line=split(/\t/,$_);
		$TE=$line[1];
		$RA{$TE}=0;
		$RPL{$TE}=0;
		$RPR{$TE}=0;
	    }
	    $n++;
	}
    print "Numero total de TES $n";
	close IN;

	$n=0;
	open (IN,"Tresults");
	while (<IN>) {
	    if ($n >0){
        ($strain, $TE, $abs, $preL, $preR)=(split /\t/)[0,1,5,12,19];
		$RA{$TE}=$RA{$TE}+$abs;
		$RPL{$TE}=$RPL{$TE}+$preL;
		$RPR{$TE}=$RPR{$TE}+$preR;
	    }
	    $n++;
	}
	my $preMax;
	my $truefreq_sum;
	for my $te ( sort keys %RA ) {
	    if( ($RPR{$te}+$RPL{$te}+$RA{$te})>=0){
            if( ($RPR{$te}+$RPL{$te}) == 0 ){
                if( $RA{$te} > $minreads && $RA{$te} < $maxreads ){
                print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\t0\n";
                }
                else{
                print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\tno_data\n";
                }
            } elsif ($RPR{$te} > 0 && $RPL{$te} == 0) {
                if ($RA{$te}+$RPR{$te} > $minreads && $RA{$te}+$RPR{$te} < $maxreads) {
                    $truefreq_sum=(($RPR{$te})/($RPR{$te}+$RA{$te}))*100;
                    print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\t$truefreq_sum\n";
                } else {
                    print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\tno_data\n";
                }
            } elsif ($RPL{$te} > 0 && $RPR{$te} == 0) {
                if ($RA{$te}+$RPL{$te} > $minreads && $RA{$te}+$RPL{$te} < $maxreads) {
                    $truefreq_sum=((+$RPL{$te})/($RPL{$te}+$RA{$te}))*100;
                    print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\t$truefreq_sum\n";
                } else {
                    print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\tno_data\n";
                }
            } else{
                if( $RA{$te} == 0 && (($RPL{$te}+$RPR{$te})/2) > $minreads && (($RPL{$te}+$RPR{$te})/2) < $maxreads){
                    print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\t100\n";
                } else {
                    if (((($RPL{$te}+$RPR{$te})/2)+$RA{$te}) > $minreads && ((($RPL{$te}+$RPR{$te})/2)+$RA{$te}) < $maxreads) {	
                        $truefreq_sum=((($RPL{$te}+$RPR{$te})/2)/((($RPL{$te}+$RPR{$te})/2)+$RA{$te}))*100;
                        print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\t$truefreq_sum\n";
                    } else {
                        print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\tno_data\n";
                    }
                }
            }
	    }
	    
	    else{
		print OUTP "$te\t$RA{$te}\t$RPL{$te}\t$RPR{$te}\tno_data\n";	
	    }
	    
	}
	
	#######################
	
    }
    else{

    # Release 3: modification in the individual frequency estimation
    
	print "frequency estimates using single strains...\n";  
	open (IN,"Tresults");  
	open (OUT,">Tfreq"); 
	print OUT "TE\tpresence_results\tpolymorphic_results\tabsence_results\ttotal_results\tTEfrequency\n";
	
	my %PRES;
	my %POLY;
	my %APOLY;
	my %PPOLY;
	my %DATA;
	 my($strain, $TE, $abs, $pre, $conclusion);
	my $n = 0;
	while (<IN>) {
	    if ($n > 0){
		my @line=split(/\t/,$_);
		$TE=$line[1];
		$PRES{$TE}=0;
		$POLY{$TE}=0;
		$APOLY{$TE}=0;
		$PPOLY{$TE}=0;
		$DATA{$TE}=0;
	    }
	    $n++;
	}
	
	$n=0;
	close IN;
	open (IN,"Tresults");
	while (<IN>) {
	    if ($n >0){
		($strain, $TE, $abs, $pre, $conclusion)=(split)[0,1,2,3,4];
		if ($conclusion eq 'present'){
		    $PRES{ $TE }++;
		}
		elsif ($conclusion =~ /^polymorphic$/){
		    $POLY{ $TE }++;
		}
		if ($conclusion ne 'no_data'){
		    $DATA{ $TE }++;
		}
	    }
	    $n++;
	}
	
	while ( my ($te, $count) = each(%PRES) ) {
	    my $fq= 0;
	    my $data = $DATA{ $te };
	    my $poly = $POLY{ $te };
	    my $abs = 0;
	    my $tot = 0;
	    if ($data >= $minpop ) {
		if( $count > 0 || $poly > 0 ){
            $fq = sprintf("%.2f",((($count+($poly/2))/$data)*100));
		    $tot = $count+$poly;
		    $abs = $data-$tot;
		}
		else{
		    $tot = $count+$poly;
		    $fq = sprintf("%.2f",(($tot/$data) *100));
		    $abs = $data-$tot;
		}
	    }
	    else{
		$fq="no_data"
	    }
	    print OUT "$te\t$count\t$poly\t$abs\t$data\t$fq\n";
	}
    }
    close OUT;
    close IN;
    chdir("../");
}

#################################################################################################################################################################
################################################################## Combine Tfreq + Tresults + TSD ############################################################### 
#################################################################################################################################################################

sub CombineAll {

    my $output = $_[0];
    my $pooled = $_[1];
    chdir("$output");
    system('pwd');
    open (OUT,">Tfreq_combine");
    print OUT "TEname\tnb_evidences_absence\tnb_evidences_presence_L\tnb_evidences_presence_R\tTEfrequency_estimates\tTSDdetection\tTEflanking_L\tTEflanking_R\n";
       
    my %HR;
    my $n = 0;
    my ($TEr, $polyL, $repeatL, $polyR, $repeatR);
    open (R,"Tresults") or die "Failed to open Tresults file !\n";
    while (<R>) {
	if ($n >0){
	    ($TEr, $polyL, $repeatL, $polyR, $repeatR)=(split)[1,8,10,15,17]; 
	    $HR{ $TEr } = [$polyL, $repeatL, $polyR, $repeatR];
	}
	$n++;
    }
    close R;
    
    my %HT;
    my $jt="";
    my $TEt="";
    my $target="";
    my $copy="";
    my $res="";
    my $preTEt="";
    my $pret="";
    my $prec=""; 
    my $prer="";

    open (T,"Tannot_TSDdetection") or die "Failed to open Tannot_TSDdetection file !\n";
    while (<T>) {
	my ($TEt, $target, $copy,$res)=(split)[0,1,2,4];
	print "$TEt $target $copy $res\n";
	if ($preTEt eq $TEt){
	    if ($prer ne $res){
		if ( $res =~/noGap/ ) {
		    $jt=$jt."/".$res;
		}
		else{
		    if (length($target) >= 2) {
			$jt=$jt."/".$res."[".$target."-".$copy."]";
		    }
		    else{
			$jt=$jt."/noGap";
		    }
		}
		$HT{ $TEt } = $jt;
	    }
	    else{
		if ($pret ne $target) {
		    if ( $res =~/noGap/ ) { 
			$jt=$jt."/".$res; 
		    }
		    else{
			if (length($target) >= 2) {
			    $jt=$jt."/".$res."[".$target."-".$copy."]";
			}
			else{
			    $jt=$jt."/noGap";
			}
		    }
		    $HT{ $TEt } = $jt;
		}
	    }
	}
	else{ 
	    if ( $res =~/noGap/ ) {  
		$jt=$res; 
	    }
	    else{
		if (length($target) >= 2) {
		    $jt=$res."[".$target."-".$copy."]";
		}
	    }
	    $HT{ $TEt } = $jt;
	}
	$preTEt = $TEt;
	$prer = $res;
	$pret = $target;
	$prec = $copy;
    }
    close T;
    
    
    my ($TEf, $na, $npL, $npR, $freq);
    my $jf="";
    my @jr;
    my $jrj="";
    if ($pooled){
	print "frequency were estimated using POOLED data...\n";  
	open (F,"Tfreq_pooled") or die "Failed to open Tfreq_pooled file !\n";
	$n = 0;
	while (<F>) { 
	    if ($n >0){
		($TEf, $na, $npL, $npR, $freq)=(split)[0,1,2,3,4];
		
		if ( ( $HT{ $TEf } ) && ( $HR{ $TEf } ) ) {
		    $jf = $HT{ $TEf };
		    $jf =~ s/noGap\///g;    
		    $jf =~ s/\/noGap//g;   
		    $jf =~ s/\/noGap\///g;    
		    @jr = $HR{ $TEf };
		    $jrj=join('/',@{$jr[0]}); ### unjoined
		    print OUT "$TEf\t$na\t$npL\t$npR\t$freq\t$jf\t$jrj\n";
		}
		else{
		    print OUT "$TEf\t$na\t$npL\t$npR\t$freq\tNA\tNA\tNA\n";
		}
	    }
	    $n++;
	}
	close F;
	
    }
    else{ 
	print "For combining all results in pooled data you need to ad -pooled flag\n"; 
    }
    
    close OUT; 
    chdir("../");
}
