# **Welcome to _T-lex 2.3_ manual**

**Release 2.3:** Maria Bogaerts <maria.bogaerts@ibe.upf-csic.es>, Josefa Gonzalez <josefa.gonzalez@ibe.upf-csic.es>

[You can find us on Github!](https://github.com/GonzalezLab/T-lex)

T-lex is a computational pipeline that detects presence and/or absence of annotated individual transposable elements (TEs) using next-generation sequencing (NGS) data. 
### **UPDATES** 
**Version 2.3:**
- In the PRESENCE module, **MAQ** is changed by **BWA-MEM**
- Filtering steps are added due to the change of mapper and adds a more accurate result
- Tresults file now has two more columns providing the _"filtered number of reads"_ in the left and in the right flaking region of the TE for more accurate pool frequency estimation
- _Absent/polymorphic_ and _present/polymorphic_ calls are considered as _no_data_ 
- Minimum and maximum number of reads are required for calculating TE frequencies from pool data
- Minimum number of strains with calls for each TEs are required to calculate frequencies from individual strains

**Previous versions:**    
	    **T-lex** _Fiston-Lavier et al., 2015_ |
[Paper](https://academic.oup.com/nar/article/43/4/e22/241098) | 
[Repository](http://petrov.stanford.edu/cgi-bin/Tlex.html) |

### **PRE-REQUISITES**
1. Unix system with Perl version 5.10.0 or higher | [Perl](http://www.perl.org/get.html) |
2. RepeatMasker and libraries open-3.2.8 or higher (_Smit 1996-2010_) | [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) | [Libraries](http://www.girinst.org) |
3. SHRIMP2 - Short Read Mapping Package- Release 2.2.1/Oct. 31, 2011 (_David et al 2011_) | [SHRIMP2](http://compbio.cs.toronto.edu/shrimp) |
4. BLAT version 34 or greater (_Kent 2002_) | [BLAT](http://www.soe.ucsc.edu/~kent) |
5. Samtools version 1.4 or higher (_Li et al 2009_)
6. Bcftools version 1.3.1 or higher (_Li et al 2009_)
7. BWA version 0.7.15 or higher (_Li and Durbin 2009_)

**Only for the Target Site Duplication (TSD) detection:**

8. Phrap version 1.090518 ("Phil's Revised Assembly Program" (_Green, 1999_) | [Phrap]( http://www.phrap.org) |
      - For downloading this program, you may need to contact the original authors to get permission 
9. FastaGrep developed by the Department of Bioinformatics at the University of Tartu & Estonian Biocentre | [FastaGrep](http://bioinfo.ebc.ee) |


### **INPUT FILES** 
T-lex requires at least four input data: 
TE list             - the list of the TEs to be analyzed with one TE identifier per line (e.g for a TE in Drosophila 'FBti0019293'). 
TE annotations      - the annotations of these TEs. Tabulated file with five columns: 
                        TE name / location / start nucleotide position/ end nucleotide position / strand (+ or -),
                        e.g. 'FBti0019293 3R 405387 406627 - '. 
reference genome    - the reference genomic sequences where TEs have been annotated in fasta format (e.g. '>3R XXXXX'). 
NGS data            - the NGS data in official [fastq format](http://en.wikipedia.org/wiki/Fastq). 

NGS data has to be stored in one directory in which each subdirectory corresponds to a single strain NGS data such as: 

       [input strain directory]/
             [strain name 1]/
                   [strain name 1].fastq
             [strain name 2]/
                   [strain name 2].fastq

To handle paired-end reads, separate the reads from a same pair ([strain name]) in two fastq files ([strain name]_1.fastq and [strain name]_2.fastq). Each fastq file name should be renamed as 'XXX_X.fastq'. 

       [input strain directory]/
               [strain name]/
                   [strain name]_1.fastq
                   [strain name]_2.fastq

### **RUNNING T-LEX / USAGE**
T-lex2.3 is composed by five different modules. All the modules can be run in the same job, excepting calculating frequencies for individual strains and TSD detection (see OPTIONS). Each module can be also run independently in case only one part of the analysis is required.

General usage:

    tlex-open-v2.3.pl [ Options ] [ -T TE list ] [ -M TE annotations ] [ -G reference genome ] [ -R NGS data ]

First time running T-lex2.3 (for different options, see "OPTIONS/PARAMETERS"):

1. General run for one strain: TE-filtering, presence, absence, combine
        
        perl tlex-open-v2.3.pl -O projectname -A 95 -pairends yes -s 'drosophila' -T TElist_file.txt -M TEannotation_file.txt -G genome_file.fa -R path_to_input/strain
    >**Output:** path_to_folder/tlex_projectname/Tresults

 2. Frequency estimation in several individual strains
        
        perl tlex-open-v2.3.pl -combData
    >**Output:** path_to_folder/Tfreqs_output/Tresults  

        perl tlex-open-v2.3.pl -freq
    >**Output:** path_to_folder/Tfreqs_output/Tfreq

3. General run in a pool strain including frequency estimation

        perl tlex-open-v2.3.pl -O projectname -A 95 -pairends yes -s 'drosophila' -pooled -T TElist_file.txt -M TEannotation_file.txt -G genome_file.fa -R path_to_input/strain

    >**Output:** path_to_folder/tlex_projectname/Tresults
    >**Output:** path_to_folder/tlex_projectname/Tfreq

4. TSD analysis in any strain (requires ABSENT module alredy run)

        perl tlex-open-v2.3.pl -tsd -align -p -O projectname -A 95 -pairends yes -s 'drosophila' -T TElist_file.txt -M TEannotation_file.txt -G genome_file.fa -R path_to_input/strain

    >**Output:** path_to_folder/Tannot_TSDdetection       

If only one of the Presence or Absence modules is run, results are obtained in the following paths:

        path_to_folder/tlex_projectname/Tresults_present
        path_to_folder/tlex_projectname/Tresults_absent


After first time run we recommend the user:

* T-lex generates a list with filtered TEs, TElist_file.txt_new, which is located in the same folder as the original input TE list. We recommend to use this new file as input list parameter, and bypass the TE-analysis with the command -noFilterTE (see options).
* T-lex generates a new file according with the TSD detection. We recommend changing the annotation input file according with this information and use this new file as input annotation parameter. For following runs, the TSD module could be bypassed.

### **OPTIONS**

| Parameters | Type | Description |
|----|----|---- |
| **General parameters:** |
| -A | int | maximum read length in the data set in bp ( default = 100 ) |
| -O     |string     | project name
| -R     |string     | NGS data in FASTQ path
| -G     |string     | Path to reference genome
| -T     |string     | Path to TE list
| -M     |string     | Path to TE annotation file
| -noclean |         |   keep the intermediate files
| -binreads |        |   bypasses the generation of data in the absence and presence module
| -h or -help |      |   display this help
| **For the TE filtering step:** |
| -s  | string  | name of the species studied ( for RepeatMasker, e.g. 'drosophila' ) |
| -d  |  float  |   minimum repeat density at the flanking regions of the TEs ( default: 0.5 (50%) ) |
| -noFilterTE  | int  | do not filter TEs |
|  **Module specific parameters:** |
| _PRESENCE_ |
|-q |      | launch only the presence detection approach |
|-j | int  |length of the junction sequences to extract in bp ( default: 1000 ) |
|-b | int  | length of the internal region of the TE in bp ( default: 60 )|
|-limp      | int |  minimum match length required with the TE sequence in bp ( default: 15 )|
|-id        | int |  minimum sequence identity required with the TE sequence in % ( default: 95 )|
|-pairends  |       string   |   PE mapping at 'no' by default ('yes' or 'no')|
|-minQ      | int  | minimum quality Phrep score required to select a read ( default: 30 )|
|-processes |       int  | number of processes (used only for the NGS data reformatting; default: 1 )|
| _ABSENCE_ | 
| -p        |       | launch only the absence detection approach |
| -f        |  int  | length of the flanking sequences to extract and concatenate in bp ( default: 125 )|
| -v        |  int  | minimum read length spanning the two TE sides in bp ( default: 15 )|
| -pairends |        string  |    PE mapping at 'no' by default ('yes' or 'no')|
| -lima     |  int  | minimum non-repeated region on each side of the sequence in bp ( default: 10 )    |       
| -binref   |       | bypasses the generation of the new reference genome without insertions|
| _COMBINE_ |
| -combRes |    |  combine the presence/absence results from one strain ( REQUIRED: PRESENCE and ABSENCE modules ) (this parameter only needs to be specified if you have run the PRESENCE and ABSENCE modules independently and you want to combine those results) |
| -combData |  |  combine the presence/absence results from different strains ( REQUIRED: at least two individual strains Tresults data file )(this parameter can only be run once there are more than 1 individual strains already run with the normal pipeline)|
| -combAll |   |  combine the frequency estimates with the TE breakpoint analysis ( for pooled data only )|
| _TSD (requires the ABSENCE module '-p' )_ |
| -align  |    |  return the multiple alignments|
| -tsd  |     |   return the TSDs for the TE insertions detected as absent|
| _POOL FREQUENCIES_ |
| -freq   |         |         return the TE frequency based on the given strains |
| -pooled  |         |        return the TE frequency based on pooled data ( to use with the opti| -freq when running only this module )|
| -minR  |     int  |   minimum number of reads to calculate frequency ( default: 3 )|
| -maxR |      int   |  maximum number of reads to calculate frequency ( default: 90 )|
| _INDIVIDUAL FREQUENCIES (after using combData for combining all individual strain results)_ |
| -freq | | return the TE frequency based on the given strains|
| -minP | | minimum number of TE calls to estimate the TE frequency ( default: 1 )|

### **EXAMPLES++
A small dataset example is provided. This dataset correspond only to a few TEs in an individual Drosophila melanogaster strain. Run the following command lines in the "example" folder.

    perl tlex-open-v2.3.pl -O example -A 95 -pairends yes -noFilterTE -T TElist_example.txt -M TEannotation_example.txt -G genome_example.fa -R fastq_files/example

For TSD detection:

    perl tlex-open-v2.3.pl -tsd -align -p -O example -A 95 -pairends yes -noFilterTE -T TElist_example.txt -M TEannotation_example.txt -G genome_example.fa -R fastq_files/example

      
### **OUTPUTS**
T-lex produces several of output directories and files. The output is stored in a working directory named by default: 'tlex_output' or 'tlex_[project name]'. By default, only the final results (the 'Tresults' file) and the data necessary for the manual curation (the 'Tanalysis' and 'Talign' sub-directories) are returned.

- TE flanking sequence analysis:  
The 'Tanalysis' sub-directory includes two RepeatMasker output files: One contains the submitted sequences in which the repeats have been masked (_.masked_). The repeated regions are represented by 'N's'. The other file contains the table summarizing the detected repeats (_.out_). Tanalysis also includes the detection of poly-A tails by looking for strech of 'A' or 'T' in three prime of the TE flanking sequence. The files are organized such as:

       Tanalysis/
          Tflank_checking_[length of the flanking region analyzed].fasta.masked
          Tflank_checking_[length of the flanking region analyzed].fasta.out
          Tpoly_[length of the flanking region analyzed].fasta.polyAT

- TE detection results  
The 'Tresults' file contains the TE presence/absence detection results. Information about the read and repeat density at the flanking regions of each TE are also returned. This file is composed of 20 columns described below: 

    | Column  | Field  |  Description |
    |----|---|---|
    | 1  | strain  | name of the strain  |
    | 2  | TE  | name of the TE  |
    | 3  | TE absence detection  |  TE detection result from the absence detection approach |
    | 4  | TE presence detection  |  TE detection result from the presence detection approach |
    | 5  | TE detection conclusion  | final result combining the results from the two approaches  |
    | 6  | absence read number  |  number of reads spanning the two TE flanking sequences |
    | 7  | left match length  | length of the terminal match overlapping the TE sequence of the left contig. This value can be negative if the full TE sequence is missing |
    | 8  | left match id  |  equence identity of the terminal match overlapping the TE sequence of the left contig. The number of mismatches/matches is reported in parenthesis |
    | 9  | poly(A) detection on left TE side  | poly(A) detection result after analysis of the left TE side. If no poly(A) or poly(T) detected = 'no_polyAT'  |
    | 10 | left_coverage  | mean read coverage at the left flanking side of the TE  |
    | 11 | left_repeat  | 	name of the repeat located at the left side of the TE. If no repeat is detected = 'no repeat'  |
    | 12 | presence read number left  | number of reads spanning the TE left junction flanking sequences  |
    | 13 | presence read number left filtered  | number of reads spanning the TE left junction flanking sequences taking into account for pool frequency estimation  |
    | 14 | right match length  | length of the terminal match overlapping the TE sequence of the right contig. This value can be negative if the full TE sequence is missing  |
    | 15 | right match id  | sequence identity of the terminal match overlapping the TE sequence of the right contig. The number of mismatches/matches is reported in parenthesis  |
    | 16 | poly(A) detection on right TE side  | poly(A) detection result after analysis of the right TE side. If no poly(A) or poly(T) detected = 'no_polyAT'  |
    | 17 | right_coverage  | mean read coverage at the right flanking side of the TE  |
    | 18 | right_repeat  | name of the repeat located at the right side of the TE. If no repeat is detected = 'no repeat'  |
    | 19 | presence read number right  | number of reads spanning the TE right junction flanking sequences  |
    | 20 | presence read number right filtered  | number of reads spanning the TE right junction flanking sequences taking into account for pool frequency estimation  |


- TE frequency estimates
Using the combination of the results, T-lex estimates the frequency for each TE in the population dividing the number of strains for which the TE is present by the total number of strains for which T-lex returns data. The 'Tfreq' file will be located in a folder called "Tfreqs_output" and contains 6 columns:

    | Column  | Field  | Description  |
    |---|---|---|
    | 1  | TE  | TE name  |
    | 2  | presence results  | number of strains for which T-lex classified the TE as 'present'  |
    | 3  | polymorphic results  | number of strains for which T-lex classified the TE as 'polymorphic'  |
    | 4  | absence results  | number of strains for which T-lex classified the TE as 'absent'  |
    | 5  | total results  | number of strains for which T-lex returns a result  | 
    | 6   | TE frequency  | TE frequency estimate (i.e. adding the number of strains for which the TE is 'present' and one-half times the number of 'polymorphic' strains, and dividing by the total number of strains for which T-lex returns results)  |

Using the option '-pooled', T-lex can also estimate the frequency for each TE in the population using pooled NGS data [_Fiston-Lavier et al. 2011_](https://www.ncbi.nlm.nih.gov/pubmed/21172826). This estimate is based on the number of reads supporting the absence and the presence of the TE insertions reported by T-lex. The frequency estimates are stored in the file called 'Tfreq_pooled'. This file contains five columns such as:

| Column  | Field  | Description  |
|---|---|---|
| 1  | TE  | Name of the TE  |
| 2  | absence read number  | Number of reads providing evidence of absence  |
| 3  | presence read number/left  | Number of reads providing evidence of presence on the left TE side  |
| 4  | presence read number/right  | Number of reads providing evidence of presence on the right TE side  |
| 5  | TE frequency  | TE frequency estimate(i.e. Sum of the frequency estimates using each TE detection approach (see Fiston-Lavier et al. 2011)  |


### **FAQ**
>To be completed by the users

### **REFERENCES** 
- RepeatMasker. version Open 3.2.8 A.F.A. Smit, R. Hubley & P. Green RepeatMasker [here](http://www.repeatmasker.org/).
- Jurka,J., Kapitonov,V.V., Pavlicek,A., Klonowski,P., Kohany,O., Walichiewicz,J. (2005) Repbase update, a database of eukaryotic repetitive elements. Cytogenet Genome Res 110:462_467
- Li,H., Ruan,J., Durbin,R. (2008) Mapping short DNA sequencing reads and calling variants using mapping quality scores. Genome Res. 18:1851_1858.
- Rumble,S.M., Lacroute,P., Dalca,A.V., Fiume,M., Sidow,A., Brudno,M. (2009) SHRiMP: Accurate Mapping of Short Color-space Reads. PLoS Comput. Biol. 5(5): e1000386. doi:10.1371/journal.pcbi.1000386.
- David M, Dzamba M, Lister D, Ilie L, Brudno M. (2011) SHRiMP2: sensitive yet practical SHort Read Mapping. Bioinformatics. Apr 1;27(7):1011-2. Epub 2011 Jan 28.
- Kent W.J. (2002) BLAT - the BLAST-like alignment tool. Genome Res. 12:656-664 
- de la Bastide M, McCombie WR. (2007) Assembling genomic DNA sequences with PHRAP. Curr Protoc Bioinformatics. 2007 Mar;Chapter 11:Unit11.4.
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 1754-1760. [PMID: 19451168]
- Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- Fiston-Lavier AS, Carrigan M, Petrov DA and Gonzalez J. T-LEX: A program for fast and accurate assessment of transposable element presence using next-generation sequencing data. Nuc. Acids. Res. 2011 Mar 1;39(6):e36. Epub 2010 Dec 21
- Fiston-Lavier AS, Barron MG, Petrov DA, Gonzalez J. T-lex2: genotyping, frequency estimation and re-annotation of transposable elements using single or pooled next-generation sequencing data. Nucleic Acids Res. 2015;43(4):e22. pmid:25510498

