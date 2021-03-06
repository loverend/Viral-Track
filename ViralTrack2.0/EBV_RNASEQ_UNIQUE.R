## ---------------------------
## Script name: RNA_SEQ_UNIQUE PIPELINE: EBV_HUMAN GENOME MAPPING 
## Function: Map EBV and HUMAN genes in BULK RNAseq.
## Author: Lauren Overend with some functions borrowed from Pierre Bost.
## Wellcome Trust Centre for Human Genetics
##
## Date Created: 2021 - JAN 
##
## Email: lauren.overend@oriel.ox.ac.uk
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
## Reading in Commandline Arguments. 

#module purge
#module use -a /apps/eb/dev/ivybridge/modules/all
#module load Python/3.8.2-GCCcore-9.3.0 
#module load UMI-tools/1.0.1-foss-2020a-Python-3.8.2 
#module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
#module load SAMtools/1.10-GCC-9.3.0
#module load STAR/2.7.3a-GCC-9.3.0
#module load Subread/2.0.1-GCC-9.3.0

## ---------------------------
if(nzchar(system.file(package = "optparse"))==FALSE){
  stop("optparse not installed. Terminating. \n")
}
suppressMessages(library(optparse))
parser <- OptionParser()
option_list <- list( 
  make_option(c("-n", "--nThreadmap"), action="store", default=8, type="integer", help="runThreadN for Star Mapping. Note will also be used as threads for Feature Counts [default]"),
  make_option(c("-o", "--outputdir"), action="store", default='/gpfs2/well/immune-rep/users/kvi236/VIRUS/TESTING', type="character", help="Path to output directory"),
  make_option(c("-i", "--indexgenome"), action="store", type="character", default="/well/immune-rep/shared/10X_GENOMICS/EBV_ANNOTATION/GRCh38_EBV_BRNASEQ/", help="Path to EBV reference [default]"),
  make_option(c("-s", "--nThreadsort"), action="store", type="integer", default=1, help="outBAMsortingThreadN for STAR Mapping [default] - usually < runThreadN"),
  make_option(c("-m", "--minreads"), action="store", type="integer", default=1, help="Minimum number of reads per virus prior to use in QC analysis on[default]"),
  make_option(c("-t", "--thresholdmappedreads"), action="store", type="integer", default=50, help="Minimum number of reads per virus to pass filtering [default]"),
  make_option(c("-b", "--bins"), action="store", type="integer", default=50, help="outBAMsortingBinsN for STAR Mapping [default]"),
  make_option(c("-f", "--fastq"), action="store", type="character", default = "/gpfs2/well/immune-rep/users/kvi236/SEPSIS_RNASEQ/gains8033923/gains8033923.Unmapped.out.mate1.fa,/gpfs2/well/immune-rep/users/kvi236/SEPSIS_RNASEQ/gains8033923/gains8033923.Unmapped.out.mate2.fa", help="Path to input FASTQ file [default]"),
  make_option(c("-r", "--runname"), action="store", type="character", default="RNASEQ_EBV_Unique", help="Run Name [default]"),
  make_option(c("-a", "--auxfunctions"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/ViralTrackProgram/RPipeline/Viral-Track/AuxillaryFunctions/auxillary_viral_track_functions.R", help="Path to ViralTrack Auxillary Functions [default]"),
  make_option(c("-g", "--gtffile"), action="store", type="character", default="/well/immune-rep/shared/10X_GENOMICS/EBV_ANNOTATION/REFERENCE_FILES/genes_final.gtf", help="Path to GTF file [default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE))

# Failures here will cause pipeline to Fail!! 
if (is.null(opt$fastq)){
  stop("No FASTq file provided. Terminating. \n", call.=FALSE)
}
if (is.null(opt$outputdir)){
  stop("Output directory not present: must be specified. Terminating. \n", call.=FALSE)
}
if (opt$runname == "Viral_Track"){
  warning("Using Default Run-Name: Reconsider using Custom Run-Name. \n")
}
if (!dir.exists(opt$outputdir)) {
  warning("Output directory does not exist! Creating it! \n")
  dir.create(opt$outputdir)
}
if (!is.numeric(opt$nThreadmap) | opt$nThreadmap < 1 ) {
  stop("STAR Mapping Threads are Incorrect. Terminating. \n")
}
if (!is.numeric(opt$nThreadsort)) {
  stop("STAR BAM Sorting Threads are incorrect. Terminating. \n")
}
if (opt$nThreadsort > 1 ) {
  warning("STAR BAM Sorting Threads are higher than recommended (Default: 1). Consider adjusting paramaters.  \n")
}
if (!dir.exists(opt$indexgenome)) {
  stop("Index Genome Directory Does Not Exist. Terminating. \n")
}
if (opt$nThreadsort > opt$nThreadmap ) {
  warning("Bam Sorting Threads exceeds Bam Mapping Threads. May be incompatible with memory request. \n")
}
## Check all packages are installed. 
if(nzchar(system.file(package = "Biostrings"))==FALSE){
  stop("Biostrings not installed. Terminating. \n")
} 
if(nzchar(system.file(package = "ShortRead"))==FALSE){
  stop("ShortRead not installed. Terminating. \n")
}
if(nzchar(system.file(package = "doParallel"))==FALSE){
  stop("doParallel not installed. Terminating. \n")
}
if(nzchar(system.file(package = "GenomicAlignments"))==FALSE){
  stop("GenomicAlignments not installed. Terminating. \n")
}
if(nzchar(system.file(package = "ggplot2"))==FALSE){
  stop("ggplot2 not installed. Terminating. \n")
}
if(nzchar(system.file(package = "Matrix"))==FALSE){
  stop("Matrix not installed. Terminating. \n")
}
if(nzchar(system.file(package = "cowplot"))==FALSE){
  stop("cowplot not installed. Terminating. \n")
}
if(nzchar(system.file(package = "stringr"))==FALSE){
  stop("stringr not installed. Terminating. \n")
}
if(nzchar(system.file(package = "reshape2"))==FALSE){
  stop("reshape2 not installed. Terminating. \n")
}

##-------------------------------------------------------
## Load Packages: 
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))
suppressMessages(library(ShortRead))
suppressMessages(library(doParallel))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(reshape2))
##-------------------------------------------------------

## Additional Functions and Log file
## Notin Function:
`%notin%` <- Negate(`%in%`)

## Setting up log.file: 
if (length(unlist(str_split(opt$f, ","))) >= 2) {
	cat("More than one input file detected: R1 and R2 File input, generating file name. \n")
	x <- unlist(str_split(opt$f, ","))
	x <- x[1]
	name <- unlist(strsplit(x,"/",fixed = T))
	sample_name <- name[length(name)]
	sample_name = gsub('.fastq|.fa|.fq|.gz','',sample_name) 
} else {
	name <- unlist(strsplit(opt$f,"/",fixed = T))
	sample_name <- name[length(name)]
	sample_name = gsub('.fastq|.fa|.fq|.gz','',sample_name)
} 
log <-  paste0(opt$outputdir, "/RNA_EBV_UniqueMapping_", sample_name, ".log")

##-------------------------------------------------------

##  SETTING UP LOG HEADER 
cat("RNA_BULK_UniqueMapping by Lauren Overend. \n", file=log, append=TRUE)
cat(paste0("Run Name: ", opt$runname, "\n"), file=log, append=TRUE)
start_time_1 <- Sys.time()
cat(paste0("Start time: ", start_time_1, "\n"), file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
##-------------------------------------------------------

cat(paste0("Analysing input files. \n"), file=log, append=TRUE)

## Check Validiamety of Input FASTQ File
## Function for testing file input:
testfiles  <- function(optfastq){
List_target_path = c()
if (!is.null(optfastq)) {
  if(file.exists(optfastq)){
    cat("FASTQ File present. \n", file=log, append=TRUE) 
    if(any(grepl(".fa|.fq|.fasta", optfastq))==TRUE){
      cat("FASTQ File Type is Valid. \n", file=log, append=TRUE)
      List_target_path = optfastq
    } else {
      stop("Fastq File Provided is Not of Type '.fasta/.fq/.fa'. Terminating. \n", file=log, append=TRUE)
    }
  } else {
    stop("FASTQ Provided but Path is Invalid. Terminating. \n", file=log, append=TRUE)
  }
} else {
   stop("No FASTQ File Provided. Terminating. \n", file=log, append=TRUE)
}
}

## Test the file path and wether it is a fasta file. 
if (length(unlist(str_split(opt$f, ","))) >= 2) {
	cat("More than one input file detected: R1 and R2 File input. \n", file=log, append=TRUE)
	x <- unlist(str_split(opt$f, ","))
	for (i in x){
		testfiles(i)
	}
	cat("All files Valid. \n", file=log, append=TRUE)
} else {
	testfiles(opt$f)
} 

##-------------------------------------------------------
## Checking the parameters values
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("RUN PARAMETERS. \n", file=log, append=TRUE)
cat(paste0("Output directory: ", opt$outputdir, "\n"), file=log, append=TRUE)
cat(paste0("FASTQ File: ", opt$fastq, "\n"), file=log, append=TRUE)
cat(paste0("--nThreadmap: ", opt$nThreadmap, "\n"), file=log, append = TRUE)
cat(paste0("--outputdir: ", opt$outputdir, "\n"), file=log, append = TRUE)
cat(paste0("--indexgenome: ", opt$indexgenome, "\n"), file=log, append = TRUE)
cat(paste0("--thresholdmappedreads: ", opt$t, "\n"), file=log, append = TRUE)
cat(paste0("--nThreadsort: ", opt$nThreadsort, "\n"), file=log, append = TRUE)
cat(paste0("--thresholdmappedreads: ", opt$t, "\n"), file=log, append = TRUE)
cat(paste0("--bins: ", opt$bins, "\n"), file=log, append = TRUE)
cat(paste0("--fastq: ", opt$fastq, "\n"), file=log, append = TRUE)
cat(paste0("--runname: ", opt$runname, "\n"), file=log, append = TRUE)
cat(paste0("--gtffile: ", opt$gtffile, "\n"), file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

##-------------------------------------------------------
## Renaming Paramaters to original variable names as in ViralTrack1.0 and sourcing auxillary funtions
N_thread = opt$nThreadmap 
N_thread_sort = opt$nThreadsort
N_bins = opt$bins
Output_directory = paste0(opt$outputdir, "/RNA_EBV_UNIQUE_Mapping_Analysis")
dir.create(Output_directory)
Name_run = opt$runname    
Index_genome = opt$indexgenome 
path_to_gtf = opt$gtffile
source(opt$auxfunctions) 


## ------------------------------------------------------------------------------------
## Mapping: 
## Setting Up Directory: 
is_gz_file = any(grepl(pattern = ".gz", opt$f))
name_target = sample_name  #Cleaning the name to get the original Amplification batch number
temp_output_dir = paste0(Output_directory, "/", name_target)
dir.create(temp_output_dir)


## ------------------------------------------------------------------------------------
## MAPPING 
cat(paste0("Mapping: ",name_target,".fastq file \n"), file=log, append = TRUE)
start_time <- Sys.time()
cat(paste0("Start time: ", start_time, "\n"), file=log, append=TRUE)
name_prefix = paste0(temp_output_dir, "/", name_target)
  
#We construct a complex command 
STAR_mapping_command = paste("STAR --runThreadN ",N_thread," --outBAMsortingThreadN ",N_thread_sort," --outBAMsortingBinsN ", N_bins, " --genomeDir ",Index_genome," --readFilesIn ", opt$fastq," --outSAMattributes NH HI AS nM NM XS ",
                               "--outFileNamePrefix ",name_prefix," --outSAMtype BAM SortedByCoordinate --twopassMode Basic ",
                               "--outFilterMatchNmin 35 --outMultimapperOrder Random --runRNGseed 1 --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6 > ", name_prefix, "_STAR_MAPPING.log", sep="")

#If the file is in the format .gz then we need to add an additional paramter :
if (is_gz_file==TRUE) {
    STAR_mapping_command = paste(STAR_mapping_command, "--readFilesCommand zcat", sep=" ") #Allowing to read .gz files
    }
## Launching STAR COMMAND
system(STAR_mapping_command)
end_time <- Sys.time()
## DOCUMENTING 
cat(paste0("Is File Gzipped: ", is_gz_file, "\n"), file=log, append=TRUE)
cat(paste0("End time: ", end_time, "\n"), file=log, append=TRUE)
cat(paste0("Mapping: ",name_target,".fastq DONE! \n"), file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

#-----------------------------------------------------------------------------

cat("Commencing QC analysis: \n", file=log, append = TRUE)
start_time <- Sys.time()
cat(paste0("QC START TIME: ", start_time, "\n"), file=log, append=TRUE)

# 1. Indexing the STAR Aligned.sortedByCoord.out.bam file:
cat("\t 2. Indexing BAM: \n", file=log, append = TRUE)
temp_sorted_bam = paste0(temp_output_dir, "/", name_target, "Aligned.sortedByCoord.out.bam")#[1]
#To begin with : the ordered .BAM file need to indexed
SAMtools_indexing_command = paste("samtools index",temp_sorted_bam)
system(SAMtools_indexing_command)
cat(paste0("\t Indexing of the bam file for ",name_target," is done. \n"), file=log, append = TRUE)

# 2. Compute The number of Mapped Reads For Each Chromosome/virus
cat("\t 3. Computing The number of Mapped Reads For Each Chromosome/Virus. \n", file=log, append = TRUE)
temp_chromosome_count_path = paste(temp_output_dir,"/Count_chromosomes.txt",sep = "")
SAMtools_chromosome_count_command = paste("samtools idxstats",temp_sorted_bam,">",temp_chromosome_count_path)
system(SAMtools_chromosome_count_command)
cat(paste0("\t Indexing and Stat file for ",name_target,"is done. \n"), file=log, append = TRUE)


### ---------------------------------------------------------------------------------------------------------------
# Creating VIRAL BAM FILES 
# Loading Count_Chromosome.txt and Filtering 
temp_chromosome_count = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(temp_chromosome_count) = c("Chromosome_length","Mapped_reads","Unknown")

# Calculate number of reads mapped to EBV 
virus <- grep("EBV", rownames(temp_chromosome_count), value=TRUE)
temp_chromosome_count <- temp_chromosome_count[rownames(temp_chromosome_count) %in% virus, ]
no_viruses <- length(temp_chromosome_count[,1])
cat(paste0("\t 4. Identified ", temp_chromosome_count$Mapped_reads, "EBV Reads in Sample \n"), file=log, append = TRUE)


# Create a sub-directory to export the sam files corresponding to EBV reads
newdir <- paste(temp_output_dir,"/EBV_BAM_file",sep = "")
dir.create(newdir) 
# Extracting Viral Read BAM file per virus into new Directory: 
cat("\t 5. Extracting BAMs files to new directories. \n", file=log, append = TRUE)
# Identify Number of Parrallel threads to use for Bam sorting depending on the input number of threads and the number of identified viruses
if (length(temp_chromosome_count[,1]) <= N_thread) {
  if(length(temp_chromosome_count[,1])==0){
    N_threadR <- 1
  } else {
    N_threadR <- length(temp_chromosome_count[,1])
  } 
  } else {N_threadR <- N_thread 
}
cat("\t 5.a  Setting up Parrallel Environment in R using: ", N_threadR, " Threads. \n", file=log, append = TRUE)
cl =makeCluster(N_threadR)
registerDoParallel(cl)
cat("\t Done. \n", file=log, append = TRUE)
# Create One Bam File Per Virus: 
foreach(i=rownames(temp_chromosome_count)) %dopar% {
  temp_export_bam_command = paste("samtools view -b ",temp_sorted_bam," \'",i,"\'"," > \'",temp_output_dir,"/EBV_BAM_file/",i,".bam\'",sep = "")
  system(temp_export_bam_command)
}

cat("\t 5.b  Viral BAM Files Extracted. \n", file=log, append = TRUE)

### ---------------------------------------------------------------------------------------------------------------
# Creating HUMAN BAM FILES
temp_chromosome_count_human = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(temp_chromosome_count_human) = c("Chromosome_length","Mapped_reads","Unknown")
virus_sequence_names <- grep("EBV", rownames(temp_chromosome_count_human), value=TRUE)
human_chromosomes <- rownames(temp_chromosome_count_human)[rownames(temp_chromosome_count_human) %notin% virus_sequence_names]
human_chromosomes <- human_chromosomes[human_chromosomes!="*"]
temp_chromosome_count_human = temp_chromosome_count_human[rownames(temp_chromosome_count_human)%in% human_chromosomes ,]
dir.create(paste(temp_output_dir,"/HUMAN_BAM_files",sep = ""))
foreach(i=rownames(temp_chromosome_count_human)) %dopar% {
  temp_export_bam_command = paste("samtools view -b ",temp_sorted_bam," \'",i,"\'"," > \'",temp_output_dir,"/HUMAN_BAM_files/",i,".bam\'",sep = "")
  system(temp_export_bam_command)
}
cat("\t 5.c Human BAM Files Extracted. \n", file=log, append = TRUE)
stopImplicitCluster()
cat("\t 5.d Terminating Parrallel Environment. \n", file=log, append = TRUE)
### ---------------------------------------------------------------------------------------------------------------
### ---------------------------------------------------------------------------------------------------------------

cat("\t 7. Calculating QC Metrics. \n", file=log, append = TRUE)
## Generation of the QC report:
dir <- paste(temp_output_dir,"/EBV_BAM_file/",sep = "")
if(length(list.files(dir))==0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t NO VIRUSES IDENFIFIED WITH MIN READS MAPPED >=", Minimal_read_mapped, ". \n"), file=log, append=TRUE)
  cat("\t NO QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  }  

if(temp_chromosome_count[,2]>0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t ", length(list.files(dir))," EBV READS IDENFIFIED . \n"), file=log, append=TRUE)
  cat("\t QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  QC_result <- NULL
## Generation of the QC report
  for(i in 1:length(rownames(temp_chromosome_count))) {
    z <- rownames(temp_chromosome_count)[i]
    BAM_file_1=readGAlignments(paste(temp_output_dir, "/EBV_BAM_file/", z, ".bam", sep=""),param=ScanBamParam(what=scanBamWhat()))
    BAM_file_2= BAM_file_1[BAM_file_1@elementMetadata$flag==0 | BAM_file_1@elementMetadata$flag==16]
	BAM_file= BAM_file_2[BAM_file_2@elementMetadata$mapq==255]	## Identifies only primary alignments
	if(length(BAM_file) >0){
		#Let's check the diversity of the reads
		#Lets use just the reads that map uniquely to calculate the statistics!!!
		Viral_reads = unique(BAM_file@elementMetadata$seq)
		Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob=TRUE)
		Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]
		
		#calculate percentage of nucleotides in the viral reads for Read Entropy
		# If only one reads present returns numeric vector rather than dataframe:
		if (class(Viral_reads_contents)=="numeric") {
		  Viral_reads_contents = matrix(Viral_reads_contents_mean,ncol = 4)
		}
		
		Viral_reads_contents_mean =colMeans(Viral_reads_contents)
		Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)
		#Caclulate the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
		Covered_genome = coverage(BAM_file)[[z]]
		Covered_genome = as.numeric(Covered_genome)
		Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
		Covered_genome = rle(sign(Covered_genome))
		Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
		
		##The mean reads quality
		Reads_quality = as.character(BAM_file@elementMetadata$qual)
		Reads_quality = PhredQuality(Reads_quality)
		Reads_quality = as(Reads_quality,"IntegerList")
		Reads_quality = as.numeric(as.matrix(Reads_quality))
		Mean_read_quality = mean(Reads_quality)
		Sd_read_quality = sd(Reads_quality)
		
		##... the number of mapped reads and unique mapped reads
		N_unique_mapped_reads = sum(BAM_file_2@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
		N_mapped_reads = length(BAM_file_2)
		Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
		
		##DUSTy score identifies low-complexity sequences, in a manner inspired by the dust implementation in BLAST
		Mean_dust_score = NA
		Percent_high_quality_reads = NA
		if ("ShortRead"%in%installed.packages()){
		  DUST_score = dustyScore(BAM_file@elementMetadata$seq)
		  Mean_dust_score = mean(DUST_score)
		  Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
		}
	} else { 
		N_unique_mapped_reads = sum(BAM_file_2@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
		N_mapped_reads = length(BAM_file_2)
		Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
		Mean_read_quality <- "NA"
		Sd_read_quality <- "NA"
		A = "NA" 
		C = "NA" 
		G = "NA" 
		T = "NA"
		Read_entropy="NA"
		Spatial_distribution="NA" 
		Longest_contig="NA"
        Mean_dust_score="NA"
		Percent_high_quality_reads="NA"
		}
	   
    ##Summary Statistics Per Virus 
    QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
                Mean_read_quality,Sd_read_quality,
                Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
                Mean_dust_score,Percent_high_quality_reads)
    QC_result = rbind(QC_result, QC_temp)
  }
} else {
  QC_result = data.frame(N_mapped_reads = numeric(),N_unique_mapped_reads=numeric(),Percent_uniquely_mapped=numeric(),
              Mean_read_quality=numeric(),Sd_read_quality=numeric(),
              A = numeric(), C= numeric(), G = numeric(), T = numeric(), Read_entropy=numeric(),Spatial_distribution=numeric(),Longest_contig=numeric(),
              Mean_dust_score=numeric(),Percent_high_quality_reads=numeric())
}
cat("\t Calculating VIRAL QC Metrics.... DONE!. \n", file=log, append = TRUE)


## ------------------------------------------------------------------------------------
## Now we perform filtering based ont the Calculaed QC statsitics: 
## Editied by LEO -- otherwise fails if only one virus is detected as produces numeric vector rather than dataframe. 
if (class(QC_result)=="numeric"){
  cat('Only one virus detected - output is not a dataframe: converting \n', file=log, append = TRUE)
  QC_result <- as.data.frame(t(as.data.frame(QC_result)))
}
if (length(QC_result[,1])>0){
	colnames(QC_result) = c("N_reads","N_unique_reads","Percent_uniquely_mapped",
							"Mean_read_quality","Sd_read_quality",
							c("A","C","G","T"),"Sequence_entropy","Spatial_distribution","Longest_contig",
							"DUST_score","Percent_high_quality_reads")
	QC_result <- as.data.frame(QC_result)
	QC_result <- sapply( QC_result, as.numeric )
	QC_result <- as.data.frame(t(QC_result))
	rownames(QC_result) = rownames(temp_chromosome_count)
} 
QC_result <- as.data.frame(QC_result)
QC_result$genome <- rownames(QC_result)

### ------------------------------------------------------------------------------------
## Now we need to extract information on the mapping by itself to get information about the QC
## This Requires R auxillalry Functions: 
path_to_Log_file = paste0(temp_output_dir, "/", name_target, "Log.final.out")
Mapping_information = Extraction_Log_final(path_to_Log_file)
Mean_mapping_length = Mapping_information$Length_vector[1]
cat("\t 8. Filtering on QC Metrics. \n", file=log, append = TRUE)

## DO FILTERING 
detected_virus = rownames(QC_result[QC_result$N_reads > opt$t & QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])

## Add on metric True/False filtering to QC tabel which is used later in plotting: 
if (length(QC_result[,1])>0){
	QC_result[,"PassedFiltering"] <- NA
	for (j in 1:length(QC_result$genome)){
		i <- QC_result$genome[j]
		if (i %in% detected_virus){
			QC_result$PassedFiltering[j] <-"PASS"
		} else {
			QC_result$PassedFiltering[j] <- "FAIL"
		}
	}
}

## Returning QC Metrics for viruses that passed Filtering 
Filtered_QC=QC_result[detected_virus,]
unfilteredoutfile <- paste0(temp_output_dir, "/", name_target, "QC_EBV_ANALYSIS.txt")

cat("\t 9. Exporting QC Metrics for EBV.txt files. \n", file=log, append = TRUE)
##Exporting the tables of the QC analysis
write.table(file = unfilteredoutfile, x = QC_result, quote = F,sep = "\t")
cat("\t Exporting QC Metrics... DONE! \n", file=log, append = TRUE)


## ------------------------------------------------------------------------------------
# Note removed splice table info as this didnt seem to be used anywhere else?
##------------------------------------------------------------------------------------

## Calculating some additional Statistics needed for Plots: 
## Additional info : % of reads mapped to viral vs host
cat("\t 10. Calculating Broad Metrics on Mapping. \n", file=log, append = TRUE)
Read_count_temp = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(Read_count_temp) = c("Chromosome_length","Mapped_reads","Unknown")
Read_count_temp = Read_count_temp[Read_count_temp$Mapped_reads!=0,]
Read_count_temp$Chr <- rownames(Read_count_temp)

virus_sequence_names <- grep("EBV", rownames(Read_count_temp), value=TRUE)
human_chromosomes <- rownames(Read_count_temp)[rownames(Read_count_temp) %notin% virus_sequence_names]
human_chromosomes <- human_chromosomes[human_chromosomes!="*"]

## Need to rename human chr to append chr to make it easy to grep:
for (x in human_chromosomes){
  Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == x , paste0("chr", Chr), Chr))
}

## Set as Rownames
rownames(Read_count_temp) <- Read_count_temp$Chr
Read_count_temp$Chr <- NULL

## Calculate the percentage of human reads etc 
host_mapping_count = sum(Read_count_temp[grepl(pattern = "chr",rownames(Read_count_temp)),"Mapped_reads"])
viral_mapping_count = sum(Read_count_temp[grepl(pattern = "EBV",rownames(Read_count_temp)),"Mapped_reads"])
total_mapping = viral_mapping_count + host_mapping_count
Ratio_host_virus = matrix(data = c(host_mapping_count,viral_mapping_count)/total_mapping,ncol = 1)*100
percent_viral = viral_mapping_count / total_mapping * 100 

## Number of uniquely mapped reads and other reads for the filtered virus 
Mapping_selected_virus = QC_result
rownames(Mapping_selected_virus) <- rownames(QC_result)
if(length(Mapping_selected_virus[, 1])>0) {
	Mapping_selected_virus = Mapping_selected_virus[order(Mapping_selected_virus$N_unique_reads,decreasing = TRUE),]
}
cat("\t 10. Calculating Broad Metrics on Mapping...DONE \n", file=log, append = TRUE)
## ------------------------------------------------------------------------------------
## Plotting the Statisitics 

cat("\t 11. Plotting QC plots to pdf. \n", file=log, append = TRUE)
Mapping_selected_virus$Name_id <- rownames(Mapping_selected_virus)
Mapping_selected_virus$Name_sequence <- rownames(Mapping_selected_virus) 
Mapping_selected_virus$Complete_segment_name <- rownames(Mapping_selected_virus)

if(length(Mapping_selected_virus[, 1])==0) {
  Mapping_selected_virus <- data.frame(N_reads = numeric(), N_unique_reads = numeric(), Name_id = character(), Name_sequence = character(), Complete_segment_name=character())
}


if(length(Mapping_selected_virus[, 1])>0) {
	Mapping_selected_virus$fullname <- Mapping_selected_virus$Name_id
	Nucleotide_usage <- data.frame(Mapping_selected_virus[, c("A", "C", "T", "G", "fullname")])
	Nucleotide_usage <- melt(Nucleotide_usage, id.vars=c("fullname"), measure.vars=c("A", "C", "T", "G"), value.name="NTUsage")
}


pdf_name <-paste0(temp_output_dir, "/", name_target, "_Summary.pdf")
## Open PDF file: 
pdf(pdf_name, height = 20, width = 20)

# Plotting the proportion of uniquely mapped reas, unmapped etc...
Mapping_result <- data.frame(Mapping_information$Mapping_result)
colnames(Mapping_result) <- "Percentage"
Mapping_result[, "Mapping"] <- "Alignment"
Mapping_result[, "Alignment"] <- "NA"
Mapping_result$Alignment[1] <- "Uniquely Mapped"
Mapping_result$Alignment[2] <- "Mapped to multiple Loci"
Mapping_result$Alignment[3] <- "Unmapped"
mapping_plot <- ggplot(Mapping_result, aes(fill=Alignment, y=Percentage, x=Mapping)) + geom_bar(position="stack", stat="identity") + ylab("Percentage of Reads (%)") + xlab(" ") + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + ggtitle("Percentage of Reads Mapped")
  
# Size of mapping, insertion and deletion
Mapping_summary <- data.frame(Mapping_information$Length_vector)
colnames(Mapping_summary) <- "Average_Length"
Mapping_summary[, "Category"] <- "Category"
Mapping_summary[, "Type"] <- "NA"
Mapping_summary$Type[1] <- "Mapped Read"
Mapping_summary$Type[2] <- "Insertion"
Mapping_summary$Type[3] <- "Deletion"
mapping_summary <- ggplot(Mapping_summary, aes(x=Type, y=Average_Length)) + geom_bar(stat="identity") + ylab("Average Length in Nucleotides") + xlab(" ") + theme_classic() + ggtitle("Mapped Read Length") 
  
# Rates as percentages
Mapping_rate <- data.frame(Mapping_information$Rate_vector)
colnames(Mapping_rate) <- "Rate"
Mapping_rate[, "Category"] <- "Category"
Mapping_rate[, "Type"] <- "NA"
Mapping_rate$Type[1] <- "Mismatch"
Mapping_rate$Type[2] <- "Insertion"
Mapping_rate$Type[3] <- "Deletion"
Mapping_rate <- ggplot(Mapping_rate, aes(x=Type, y=Rate)) + geom_bar(stat="identity") + ylab("Rate per base(%)") + xlab(" ") + theme_classic() + ggtitle("Rate of Mismatch")
  
# Ratio 
Ratio_host_virus <- data.frame(Ratio_host_virus)
colnames(Ratio_host_virus) <- "Percentage"
Ratio_host_virus[, "Mapping"] <- "Alignment"
Ratio_host_virus[, "Class"] <- "NA"
Ratio_host_virus$Class[1] <- "Host"
Ratio_host_virus$Class[2] <- "Virus"
mapping_host_virus <- ggplot(Ratio_host_virus, aes(fill=Class, y=Percentage, x=Mapping)) + geom_bar(position="stack", stat="identity") + ylab("Mapped Reads (%)") + xlab(" ") + theme_classic() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + ggtitle("Origin of Mapped Reads")
 

# Making sure binary categories availible for pass/fail filtering.  
QC_result$PassedFiltering <- factor(QC_result$PassedFiltering)
if ("PASS" %notin% levels(QC_result$PassedFiltering)){
	levels(QC_result$PassedFiltering) <- c(levels(QC_result$PassedFiltering), "PASS")
	QC_result$PassedFiltering <- relevel(QC_result$PassedFiltering, "PASS")
}

if ("FAIL" %notin% levels(QC_result$PassedFiltering)){
	levels(QC_result$PassedFiltering) <- c(levels(QC_result$PassedFiltering), "FAIL")
	QC_result$PassedFiltering <- relevel(QC_result$PassedFiltering, "PASS")
}


# Plot Number of Reads > 50 (filtering threshold)
if(length(rownames(QC_result))>0){
  p1 <- ggplot(QC_result, aes(shape=PassedFiltering, color=genome, x=N_reads, y=(Spatial_distribution*100))) + geom_point() + scale_color_discrete(drop=FALSE) + scale_shape(drop=FALSE) + labs(color="Genome", shape="Passed Filtering") + theme_classic() + xlab("Number of Mapped Reads") + ylab("% Mapped genome") + geom_vline(xintercept=(opt$t), linetype="dashed", color="grey")+ geom_hline(yintercept=5, linetype="dashed", color="grey") + ylim(0, 100) + ggtitle("Viral Summary Unique Reads: Multimapping Read Count")
  p2 <- ggplot(QC_result, aes(shape=PassedFiltering, color=genome, x=N_unique_reads, y=(Spatial_distribution*100))) + geom_point() + scale_color_discrete(drop=FALSE) + scale_shape(drop=FALSE) + labs(color="Genome", shape="Passed Filtering") + theme_classic() + xlab("Number of Uniquely Mapped Reads") + ylab("% Mapped genome") + geom_hline(yintercept=5, linetype="dashed", color="grey") + ylim(0, 100) + ggtitle("Viral Summary Unique Reads: Unique Read Count")
  p3 <- ggplot(QC_result, aes(shape=PassedFiltering, color=genome, x=Sequence_entropy, y=(Spatial_distribution*100))) + geom_point() + scale_color_discrete(drop=FALSE) + scale_shape(drop=FALSE) + labs(color="Genome", shape="Passed Filtering") + theme_classic() + xlab("Sequence Entropy") + ylab("% Mapped genome") + geom_hline(yintercept=5, linetype="dashed", color="grey") + ylim(0, 100) + geom_vline(xintercept=1.2, linetype="dashed", color="grey") + ggtitle("Viral Summary Unique Reads: Unique Read Sequence Complexity")
  p4 <- ggplot(QC_result, aes(shape=PassedFiltering, color=genome, x=Longest_contig, y=DUST_score)) + geom_point() + scale_color_discrete(drop=FALSE) + scale_shape(drop=FALSE) + labs(color="Genome", shape="Passed Filtering") + theme_classic() + xlab("Longest Contig (nt)") + ylab("DUST Score") + geom_vline(xintercept=(3*Mean_mapping_length), linetype="dashed", color="grey") + ggtitle("Viral Summary Unique Reads: DUST Score")
  p5 <- ggplot(QC_result, aes(shape=PassedFiltering, color=genome, x=Mean_read_quality, y=Sd_read_quality)) + geom_point() + scale_color_discrete(drop=FALSE) + scale_shape(drop=FALSE) + labs(color="Genome", shape="Passed Filtering") + theme_classic() + xlab("Mean Read Quality") + ylab("SD Read Quality") + ggtitle("Viral Summary Unique Reads: Read Quality")
  plot(plot_grid(mapping_plot, mapping_host_virus, mapping_summary, Mapping_rate, p1, p2, p3, p4, p5, ncol=3, labels = "AUTO"))
} 
if(length(rownames(QC_result))==0){
	plot(plot_grid(mapping_plot, mapping_host_virus, mapping_summary, Mapping_rate,  ncol=3, labels = "AUTO"))
} 
#Number of reads for each filtered virus
if (length(Mapping_selected_virus[, 1])>0) {
 t1 <- ggplot(Mapping_selected_virus, aes(x = Complete_segment_name, y = N_reads, fill=Name_id)) + geom_bar(stat="identity") + theme_classic() + xlab("Virus") + ylab("Number of Mapped Reads") + coord_flip() + labs(fill="NC Identifier")
 t2 <- ggplot(Mapping_selected_virus, aes(x = Complete_segment_name, y = N_unique_reads, fill=Name_id)) + geom_bar(stat="identity") + theme_classic() + xlab("Virus") + ylab("Number of Uniquely Mapped Reads") + coord_flip() + labs(fill="NC Identifier")
 t4 <- ggplot(Nucleotide_usage, aes(fill=variable, y=NTUsage, x=fullname)) + geom_bar(position="fill", stat="identity") + theme_classic() + xlab("Virus") + ylab("Nucleotide Usage") + labs(fill="Nucleotide") + coord_flip()
 plot(plot_grid(t1, t2, t4, nrow=2, ncol=2, labels = "AUTO"))
} 
dev.off() 


cat("\t 11. Plotting QC statistics.. Done. \n", file=log, append = TRUE)
cat("\n ----------------------------------------------\n", file=log, append=TRUE)
cat("QC PDF and Summary Statistics Done and Saved.  \n", file=log, append = TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)


## -----------------------------------------------------------------------------------


## END of Log Script: 
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("ViralTrack_Scanning.2.0 COMPLETE \n", file=log, append=TRUE)
end_time <- Sys.time()
cat(paste0("End time: ", end_time, "\n"), file=log, append=TRUE)
Duration <- start_time_1 - end_time
cat(paste0("End time: ", Duration, "\n"), file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)


## ------------------------------------------------------------------------------------
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("ViralTrack_FeatureCounting  \n", file=log, append=TRUE)
start_time <- Sys.time()
cat("Start Time: ", start_time, "\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

## -----------------------------------------------------------------------------------
cat("Aggregating All Detected Viral and Human BAM files into Reads to Demultipled. \n", file=log, append=TRUE)

## First aggregating all the reads from the detected viruses:
List_bam_files =c()
## Identify which viral bam files we want to demultiplex (those for filtered viruses)
Identified_viral_fragments <- detected_virus  
## Assuming viral reads are present: 
for (segment_temp in QC_result$genome) {
  List_bam_files = c(List_bam_files, paste(temp_output_dir,"/EBV_BAM_file/",segment_temp,".bam",sep = ""))
}
## Also we want to demultiplex human Reads
path_to_human <- paste0(temp_output_dir, "/HUMAN_BAM_files")
List_bam_files <- c(List_bam_files, list.files(path_to_human, recursive = TRUE, full.name=TRUE))
List_bam_files = paste("\'",List_bam_files,"\'",sep = "")
  
## Merge these bam files into a Read to demultiplex file:   
command_merge = base::paste("samtools merge ", temp_output_dir,"/Reads_to_demultiplex.bam -f ",paste(List_bam_files,collapse = " "),sep="")
system(command_merge)

cat("Aggregation Complete. \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("Performing Feature Counts on Reads To Demultiplex.  \n", file=log, append=TRUE)
start_time_f <- Sys.time()
cat(paste0("-T: Threads for featureCounts: ", N_thread, " \n"), file=log, append=TRUE)
cat("Start Time Feature Counts: ", start_time_f, "\n", file=log, append=TRUE)
## Assigning reads to transcripts using Rsubread Featurecounts
featurecommand = paste("featureCounts -T ", N_thread, " -t exon -R BAM -g gene_id -a ", path_to_gtf, " -o ", temp_output_dir, "/counts.txt ", temp_output_dir, "/Reads_to_demultiplex.bam 2> ", name_prefix, "_FEATURE_COUNTS.log", sep="")
system(featurecommand)
end_time_f <- Sys.time()
cat("End Time Feature Counts: ", end_time_f, "\n", file=log, append=TRUE)
cat("Feature Counts Complete.  \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

#--------------------------------------------------------------------------------
## Read in counts.txt file (output from feature counts) 

RNA_counts <- read.delim(paste0(temp_output_dir, "/counts.txt"), header = TRUE, sep = "\t", skip="1")
colnames(RNA_counts)[7] <- "count"

## -------------------------------------------------------------------
## Summarised the information in a counts file 

cat("\n 9. Exporting RNA_SEQ COUNTS to files. \n", file=log, append = TRUE)
##Exporting the tables of the QC analysis

EBV_COUNTS_SUMMARY <- paste0(temp_output_dir, "/EBV_COUNTS_SUMMARY_", sample_name, ".txt")

write.table(file = EBV_COUNTS_SUMMARY, x = RNA_counts, quote = F,sep = "\t")
cat("\n Exporting RNA_SEQ Counts summary... DONE! \n", file=log, append = TRUE)

## -----------------------------------------------------------------------------------


## Now we just move files of interest into the Report directory: 
cat("Creating Report Directory. \n", file=log, append=TRUE)
Report_dir <- paste0(temp_output_dir, "/Reports")
dir.create(Report_dir)
cat("Moving Useful Files into Report Directory. \n", file=log, append=TRUE)
move <- paste0("mv ", EBV_COUNTS_SUMMARY, " ",  Report_dir )
system(move)
move <- paste0("mv ", pdf_name, " ", Report_dir )
system(move)
move <- paste0("mv ", unfilteredoutfile, " ", Report_dir )


system(move)
move <- paste0("mv ", temp_output_dir, "/Reads_to_demultiplex.bam ", Report_dir )
system(move)


features <- paste0(temp_output_dir, "/", sample_name, "_FEATURE_COUNTS.log")
move <- paste0("mv ", features, " ", Report_dir )
system(move)
star <- paste0(temp_output_dir, "/", sample_name, "Log.final.out")
move <- paste0("mv ", star, " ", Report_dir )
system(move)

star <- paste0(temp_output_dir, "/", sample_name, "Log.out")
move <- paste0("mv ", star, " ", Report_dir )
system(move)

star <- paste0(temp_output_dir, "/counts.txt.summary")
move <- paste0("mv ", star, " ", Report_dir )
system(move)

cat("Done. \n", file=log, append=TRUE)

clear <- paste0("rm  ", temp_output_dir, "/* 2> /dev/null" )
system(clear)

clear <- paste0("rm -r ", temp_output_dir, "/HUMAN_BAM_files")
system(clear)

clear <- paste0("rm -r ", temp_output_dir, "/EBV_BAM_file")
system(clear)



## -----------------------------------------------------------------------------------

end_time <- Sys.time()
Time_difference <- end_time - start_time
cat("Run Time ", Time_difference, "\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)
cat(" VIRNA_seq PIPELINE COMPLETE. \n", file=log, append=TRUE)
cat("----------------------------------------------\n", file=log, append=TRUE)

## -----------------------------------------------------------------------------------quit()





