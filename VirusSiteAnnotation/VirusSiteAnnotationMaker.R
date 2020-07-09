## ---------------------------
## Script name: VirusSiteAnnotationMaker
## Function: Generate VirusSite Annotation file for updated virus_site genomes.fa file
##           Compatible with ViralTrack
## Author: Lauren Overend (LEO). 
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
library(seqinr)
suppressMessages(library(optparse))
parser <- OptionParser()
option_list <- list( 
  make_option(c("-g", "--genomes"), action="store", type="character", default = "/gpfs2/well/immune-rep/users/kvi236/References/VIRUS_REFERENCE/genomes.fasta", help="Path to genomes .fa Downloaded from VirusSite [default]"),
  make_option(c("-o", "--output"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/References", help="Path to output directory [default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )


## Read in the VirusSite genomes.fa file:
l <- read.fasta(opt$genomes, seqtype = "DNA", as.string = TRUE, forceDNAtolower = TRUE,set.attributes = TRUE)
annotation <- unlist(getAnnot(l))
all_genomes <- names(l)
## Make an Empty Data Frame
df <- NULL 
## Fill with Relevant Columns
for (i in 1:length(all_genomes)){
  z <- all_genomes[i]
  z <- unlist(strsplit(z,"|",fixed = T))
  nc <- z[2]
  nucleo <- unlist(strsplit(z[3],split='nt', fixed=TRUE))
  x <- annotation[i]
  f <- unlist(strsplit(x, "|", fixed=TRUE))[4]
  q <- unlist(strsplit(f, ",", fixed=TRUE))[1]
  row <- c(nc, nucleo, q, f)
  df = rbind(df, row)
} 
colnames(df) <- c("Name_sequence", "Genome_length", "Virus_name", "Complete_segment_name")
corona <- c("NC_045512.2", "29903", "SARS_Cov_2", "Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome") 
df = rbind(df, corona)
## Save
file_name <- paste0(opt$output, "/Updated_VirusSite_Reference.txt")
write.table(df, file = paste0(opt$output, "/Updated_VirusSite_Reference.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)


