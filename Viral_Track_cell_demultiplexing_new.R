
## ---------------------------
##
## Script name: VIRAL TRACK: PART TWO - Viral_demutliplexing 
## Function: Demultiplex viral reads and output per single cell. 
## Author: Pierre Bost (as used in Viral TRACK paper): Updated by Lauren Overend (LEO)
##
## Date Created: 2020 - June 
##
## Email: lauren.overend@oriel.ox.ac.uk
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(Matrix)
suppressMessages(library(Rsubread))
args <- commandArgs(trailingOnly = T)


#Loading parameters
Parameter_file_path = args[1]
#Parameter_file_path = '/well/immune-rep/users/kvi236/VIRUS/VIRAL_TRACK/Parameters.txt'
Parameters = read.table(Parameter_file_path,header = F,sep = "\t")
Parameters = as.character(Parameters$V1)
Parameters = strsplit(Parameters,split = "=",fixed = T)
Parameters_names = unlist(lapply(Parameters,function(x) {x[1]}))
Parameters_values = unlist(lapply(Parameters,function(x) {x[2]}))
names(Parameters_values) = Parameters_names
Output_directory = Parameters_values["Output_directory"] #Output directory
Name_run = make.names(Parameters_values["Name_run"]) #Name of the analytical run
Viral_annotation_file = as.character(Parameters_values["Viral_annotation_file"])

#Loading list files to process
Parameter_target_file = args[2]
#Parameter_target_file = '/well/immune-rep/users/kvi236/VIRUS/VIRAL_TRACK/Files_to_process_downsampled.txt'
File_to_process = read.table(Parameter_target_file,header = F,sep = "\t")
File_to_process = as.character(File_to_process$V1)
File_to_process = unique(File_to_process)

List_names = c()
for (k in File_to_process) {
  if (file.exists(k)) {
    
    name_target = base::strsplit(k,"/",fixed = T)
    name_target = name_target[[1]]
    l = length(name_target)
    name_target = name_target[l]
    name_target = gsub('/','',name_target)
    name_target = gsub('.fastq','',name_target) #Cleaning the name to get the original Amplification batch number
    
    is_gz_file = grepl(pattern = ".gz",name_target)
    
    if (is_gz_file) {
      name_target = gsub('.gz','',name_target) #Removing if necesseray the .gz 
    }
    List_names = c(List_names,name_target)
  }
}
List_target_path = paste(Output_directory,List_names,"/",sep = "/")


#Loading all viral fragments we detected in each plate
Identified_viral_fragments = c()
for (k in List_target_path) {
  name_target = names(k)
  report = read.table(paste(k,"/QC_filtered.txt",sep = ""),header = T)
  if (nrow(report)>0) {
    Identified_viral_fragments = c(Identified_viral_fragments,rownames(report))
    cat('Viral Fragments dectected \n')
  } else {
  cat('No Viral Fragments dectected \n')
  }
}
Identified_viral_fragments = unique(Identified_viral_fragments)
Identified_viral_fragments = Identified_viral_fragments[!is.na(Identified_viral_fragments)]



## -----------------------------------------------------------------
###Finding the path to the GTF : if none fonmd -> 'pseudo GTF' created
GTF_dir = paste(Output_directory, Name_run, sep = "")
dir.create(GTF_dir)
Path_GTF = paste(Output_directory, Name_run, "/Merged_GTF.gtf", sep="")
fileConn <- file(Path_GTF, open = "w+")

#Loading the viral annotation file 
Viral_annotation = read.delim(Viral_annotation_file)
Viral_annotation = Viral_annotation[Viral_annotation$Name_sequence!=" ",]

#Selecting the viral segments 
if (length(Identified_viral_fragments) > 0 ){
  cat("Viral fragments identified, genomes will be added to gtf")
  for (k in 1:length(Identified_viral_fragments)) {
    y <- strsplit(Identified_viral_fragments[k],split='|', fixed=TRUE)
    y = unlist(y)
    length =y[3]
    length <- strsplit(length,split='nt', fixed=TRUE)
    length = unlist(length)
    id <- y[2]
    z = paste(Identified_viral_fragments[k],'RefSeq', 'transcript', 1, length,1000,".",".", paste("gene_id ", Identified_viral_fragments[k] ,"_1;", sep =""), sep='\t')
    print(z)
    cat(z, file=fileConn, sep="\n")
    }
} else {
  cat("NO Identified Viral Fragments, only Human chromosomes will be added to GTF")
}
#adding human chromsomes to the gtf 
a = paste(1,'RefSeq', 'transcript', 1, 248956422, 1000,".",".", paste("gene_id ", "1" ,"_1;", sep =""), sep='\t')
b = paste(2,'RefSeq', 'transcript', 1, 242193529, 1000,".",".", paste("gene_id ", "2" ,"_1;", sep =""), sep='\t')
c = paste(3,'RefSeq', 'transcript', 1, 198295559, 1000,".",".", paste("gene_id ", "3" ,"_1;", sep =""), sep='\t')
d = paste(4,'RefSeq', 'transcript', 1, 190214555, 1000,".",".", paste("gene_id ", "4" ,"_1;", sep =""), sep='\t')
e = paste(5,'RefSeq', 'transcript', 1, 181538259, 1000,".",".", paste("gene_id ", "5" ,"_1;", sep =""), sep='\t')
f = paste(6,'RefSeq', 'transcript', 1, 170805979, 1000,".",".", paste("gene_id ", "6" ,"_1;", sep =""), sep='\t')
g = paste(7,'RefSeq', 'transcript', 1, 159345973, 1000,".",".", paste("gene_id ", "7" ,"_1;", sep =""), sep='\t')
h = paste(8,'RefSeq', 'transcript', 1, 145138636, 1000,".",".", paste("gene_id ", "8" ,"_1;", sep =""), sep='\t')
i = paste(9,'RefSeq', 'transcript', 1, 138394717, 1000,".",".", paste("gene_id ", "9" ,"_1;", sep =""), sep='\t')
j = paste(10,'RefSeq', 'transcript', 1, 133797422, 1000,".",".", paste("gene_id ", "10" ,"_1;", sep =""), sep='\t')
k = paste(11,'RefSeq', 'transcript', 1, 135086622, 1000,".",".", paste("gene_id ", "11" ,"_1;", sep =""), sep='\t')
l = paste(12,'RefSeq', 'transcript', 1, 133275309, 1000,".",".", paste("gene_id ", "12" ,"_1;", sep =""), sep='\t')
m = paste(13,'RefSeq', 'transcript', 1, 114364328, 1000,".",".", paste("gene_id ", "13" ,"_1;", sep =""), sep='\t')
n = paste(14,'RefSeq', 'transcript', 1, 107043718, 1000,".",".", paste("gene_id ", "14" ,"_1;", sep =""), sep='\t')
o = paste(15,'RefSeq', 'transcript', 1, 101991189, 1000,".",".", paste("gene_id ", "15" ,"_1;", sep =""), sep='\t')
p = paste(16,'RefSeq', 'transcript', 1, 90338345, 1000,".",".", paste("gene_id ", "16" ,"_1;", sep =""), sep='\t')
q = paste(17,'RefSeq', 'transcript', 1, 83257441, 1000,".",".", paste("gene_id ", "17" ,"_1;", sep =""), sep='\t')
r = paste(18,'RefSeq', 'transcript', 1, 80373285, 1000,".",".", paste("gene_id ", "18" ,"_1;", sep =""), sep='\t')
s = paste(19,'RefSeq', 'transcript', 1, 58617616, 1000,".",".", paste("gene_id ", "19" ,"_1;", sep =""), sep='\t')
t = paste(20,'RefSeq', 'transcript', 1, 64444167, 1000,".",".", paste("gene_id ", "20" ,"_1;", sep =""), sep='\t')
u = paste(21,'RefSeq', 'transcript', 1, 46709983, 1000,".",".", paste("gene_id ", "21" ,"_1;", sep =""), sep='\t')
v = paste(22,'RefSeq', 'transcript', 1, 50818468, 1000,".",".", paste("gene_id ", "22" ,"_1;", sep =""), sep='\t')
w = paste("X",'RefSeq', 'transcript', 1, 156040895, 1000,".",".", paste("gene_id ", "X" ,"_1;", sep =""), sep='\t')
x = paste("Y",'RefSeq', 'transcript', 1, 57227415, 1000,".",".", paste("gene_id ", "Y","_1;", sep =""), sep='\t')
y = paste("MT",'RefSeq', 'transcript', 1, 16569, 1000,".",".", paste("gene_id ", "MT","_1;", sep =""), sep='\t')

cat(a, file=fileConn, sep="\n")
cat(b, file=fileConn, sep="\n")
cat(c, file=fileConn, sep="\n")
cat(d, file=fileConn, sep="\n")
cat(e, file=fileConn, sep="\n")
cat(f, file=fileConn, sep="\n")
cat(g, file=fileConn, sep="\n")
cat(h, file=fileConn, sep="\n")
cat(i, file=fileConn, sep="\n")
cat(j, file=fileConn, sep="\n")
cat(k, file=fileConn, sep="\n")
cat(l, file=fileConn, sep="\n")
cat(m, file=fileConn, sep="\n")
cat(n, file=fileConn, sep="\n")
cat(o, file=fileConn, sep="\n")
cat(p, file=fileConn, sep="\n")
cat(q, file=fileConn, sep="\n")
cat(r, file=fileConn, sep="\n")
cat(s, file=fileConn, sep="\n")
cat(t, file=fileConn, sep="\n")
cat(u, file=fileConn, sep="\n")
cat(v, file=fileConn, sep="\n")
cat(w, file=fileConn, sep="\n")
cat(x, file=fileConn, sep="\n")
cat(y, file=fileConn, sep="\n")
close(fileConn)

## -----------------------------------------------------------------
### Counting by itself

for (path_temp in List_target_path) {
  cat("Demultiplexing reads from sample",names(List_target_path)[List_target_path==path_temp])
  cat("\n")
  
  #First aggregating all the reads from the detected viruses:
  List_bam_files =c()
  
  #Assuming viral reads are present: 
  for (segment_temp in Identified_viral_fragments) {
    List_bam_files = c(List_bam_files, paste(path_temp,"Viral_BAM_files/",segment_temp,".bam",sep = ""))
  }
  path_to_human <- paste0(path_temp, "HUMAN_BAM_files")
  List_bam_files <- c(List_bam_files, list.files(path_to_human, recursive = TRUE, full.name=TRUE))
  List_bam_files = paste("\'",List_bam_files,"\'",sep = "")
  
  #Check if viral reads are present: 
  if(is.null(List_bam_files)){
    cat("No Viral reads present: Moving on to next samples \n")
    next 
  }else{
    cat("Viral reads present: commensing demultiplexing. \n")
  }
  
  command_merge = base::paste("samtools merge ",path_temp,"Reads_to_demultiplex.bam -f ",paste(List_bam_files,collapse = " "),sep="")
  system(command_merge)
  
  #Assigning reads to transcripts using Rsubread Featurecounts
  featurecommand = paste("featureCounts -t transcript -M --primary -R BAM -g gene_id -a ", Path_GTF, " -o ", path_temp, "counts.txt ", path_temp, "Reads_to_demultiplex.bam", sep="")
  system(featurecommand)
  
  #We now have to order and index the BAM file
  command_sort =paste("samtools sort ",path_temp,"/Reads_to_demultiplex.bam.featureCounts.bam -o ",path_temp,"/Assigned_sorted.bam",sep = "")
  system(command_sort)
 
  command_index =paste("samtools index ",path_temp,"/Assigned_sorted.bam",sep = "")
  system(command_index)
 
 #Final command : Umi-tools command
 
 #Adding UMI_tools to the environment
 

 command_umi_tools = paste("umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XT --per-cell -I ",
                           path_temp,"/Assigned_sorted.bam  -S ",path_temp, "/Expression_table.tsv --wide-format-cell-counts",sep="")

 suppressMessages(system(command_umi_tools))

 # Saving matrix as Sparse Matrixform in an outs directory:
 MTX_dir = paste(path_temp,"viral_filtered_features/", sep = "")
 dir.create(MTX_dir)
 #Read non-sparse features:
 file = paste0(path_temp, "Expression_table.tsv")
 d <- read.table(file = file, sep = '\t', header = TRUE)
 rownames(d) <- d$gene 
 d$gene <- NULL
 d <- as.data.frame(t(d))
 viral_cols <- grep("refseq", colnames(d), value=TRUE)
 human_cols <- colnames(d)[!colnames(d) %in% viral_cols]
 d$percent_human_reads <- (rowSums(d[, c(human_cols)]))/ (rowSums(d))* 100
 if(length(viral_cols)>= 1){
  cat("Viruses PRESENT. Calculating group and individual viral percentages")
  if(length(viral_cols)==1){
    cat("one virus present")
    d$percent_viral_reads <- d[, c(viral_cols)]/ (rowSums(d[, c(viral_cols, human_cols)]))* 100
  } else {
    cat("Multiple Viruses Present")
    d$percent_viral_reads <- (rowSums(d[, c(viral_cols)]))/ (rowSums(d[, c(viral_cols, human_cols)]))* 100  
  } 
 } else {
   cat("No viral reads present: percent viral = 0 ") 
   d$percent_viral_reads <- "0"
 }
 # Calculate percentage for unique viruses 
 if(length(viral_cols)>= 1){
  cat("Calculating Percentage per Virus")
  for (i in viral_cols){
    y <- strsplit(i, split='|', fixed=TRUE)
    y = unlist(y)
    name=paste0("percent_viral_reads_", y[2], y[4])
    virus <- i
    d[, name] <- (d[, c(virus)])/ (rowSums(d[, c(viral_cols, human_cols)]))* 100
  }
 } 
 d <- (t(d))
 d <- Matrix(d, sparse = TRUE)
 writeMM(obj = d, file=paste0(MTX_dir, "viral_counts.mtx"))
 j <- colnames(d)
 s <- rownames(d)
 write.table(j, file=paste0(MTX_dir, 'barcodes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
 write.table(s, file=paste0(MTX_dir, 'genomes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

 ##create reports directory: 
 MTX_dir = paste(path_temp,"viral_filtered_features/Report/", sep = "")
 dir.create(MTX_dir)
 move <- paste0("cp ", path_temp, "counts.txt.summary ", MTX_dir)
 system(move)
 move <- paste0("cp ", path_temp, "QC_report.pdf ", MTX_dir)
 system(move)
 move <- paste0("cp ", path_temp, "QC_filtered.txt ", MTX_dir)
 system(move)
}



