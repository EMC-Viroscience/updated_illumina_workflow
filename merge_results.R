#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(stringr)
library(tidyr)

#Load annotation data
data <- as.data.table(fread(snakemake@input[[1]]))
colnames(data) <- c("qseqid", "qlen", "sseqid", "slen", "approx_pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "sscinames", "stitle")

#count how many hits there are per qsqeqid;
data[,n_hits:=.N,.(qseqid)]

#for each annotation with multiple taxids (identical proteins from multiple organisms); create separate row in data table so that each staxid gets own row
data_sep <- as.data.table(separate_rows(data,staxids,convert=TRUE,sep=";"))

#Load taxonomic lineage information
taxid <- fread("/mnt/scratch2/DB/taxdump/dmp_files/rankedlineage.dmp")
taxid <- taxid[,c(1,3,5,7,9,11,13,15,17,19)]
colnames(taxid) <- c("tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom")

#Load the merged taxids
merged <- fread("/mnt/scratch2/DB/taxdump/dmp_files/merged.dmp")
merged <- merged[,c(1,3)]
colnames(merged) <- c("from","to")
merged <- merge.data.table(taxid, merged, by.x="tax_id", by.y="to")
merged[,tax_id:=from]
merged[,from:=NULL]

#Combine the regular and merged taxids for backwards compatibility
taxid <- rbind(taxid,merged)

# add taxonomic lineage information to annotation file
data_merged <- merge.data.table(data_sep, taxid, by.x="staxids", by.y="tax_id", all.x=T)

#fill in tax_name if species name is absent create column species_MV
data_merged[,species_MV:=species]
data_merged[species=="",species_MV:=tax_name]

#For each hit with multiple staxids, determine most common lineage on species level;
data_merged[,n_staxids:=.N,.(qseqid,sseqid)]   #n_staxids shows if there are multiple staxids for 1 hit, same protein several organisms

# for data with multiple staxids, do majority vote to determine most common lineage; per sseqid, count most common lineage
data_multiple <- data_merged[n_staxids>1] #only filter out the hits that have multiple staxids
data_multiple[,count:=.N,.(qseqid,sseqid,species_MV)] # count per hit (sseqid), how many counts of certain species is present (frequency would be count/n_staxids)
data_multiple <- data_multiple[order(count,decreasing=T),.SD,.(sseqid)] #order count number in decreasing order
data_multiple <- data_multiple[order(count,decreasing=T),.SD[1],.(sseqid)] #take only 1 staxid for hits with multiple staxids based on majority vote
data_multiple[,count:=NULL] #make columns compatible to later merge with data that have only 1 taxid

# extract data with 1 staxid
data_single <- data_merged[n_staxids==1]

# merge data with more than 1 staxid with data that have 1 taxid, now you have majority vote for the hits with multiple taxids + data on
data_total <- rbind(data_single,data_multiple)
data_total[,count:=.N,.(qseqid,species_MV)] # count per qseqid, how many counts of certain species is present (frequency would be count/n_hits )


#create two ways of filtering best hit; use second option
#1. filter based on majority vote
#data_major_vote <- data_total[order(count,decreasing=T),.SD[1],.(qseqid)] #based on majority vote
#data_major_vote_compare <- data_major_vote[,.(qseqid,staxids,bitscore,tax_name,species)]
#data_major_vote_compare[,method:="major_vote"]

#2. filter based on normalised bitscore --> best way to filter
#calculate normalized bitscore based on length contig
data_bit <- data_total[,contig_length:=str_extract(qseqid,"(?<=length_)[:digit:]+")]
data_bit <- data_bit[,bit_corrected:=(as.numeric(bitscore)/as.numeric(contig_length))]
data_bit <- data_bit[order(bit_corrected, decreasing=T),.SD[1],.(qseqid)]

#data_bit[,method:="normalised_bit"]
#data_compare_bit <- data_bit[,.(qseqid,staxids,bitscore,tax_name,species,method)]
#compare the majority voting way of lineage determination and filter based on normalized bitscore
#comparison <- rbind(data_major_vote_compare,data_compare_bit)
#EVENTUEEL eerst filteren op bitscore, slechte hits eruit halen en dan LCA; to be added

#Load the coverage file
coverage <- fread(snakemake@input[[2]])

#Calculate contig length and total coverage per contig and add to the annotation data
coverage <- coverage[,.(qseqid=.BY,contig_length_coverage_file=nrow(.SD),coverage=sum(V3)),V1][,.(qseqid,contig_length_coverage_file,coverage)]
data_merged <- as.data.table(merge.data.frame(data_bit, coverage, by="qseqid", all.x=T))  #here length is parsed twice ; so change that
data_merged[,average_coverage:=round((coverage/contig_length_coverage_file),2),.(qseqid)]
#data_merged[, average_coverage := round(average_coverage, 2)]
#data_merged[, n_hits := NULL]

#Load the readmapping data
readmap <- fread(snakemake@input[[3]])
colnames(readmap) <- c("read_id", "contig_id")

#Calculate how many reads mapped to each contig
readmap <- readmap[,.(n_reads=.N),contig_id]

#Add read mapping numbers to the data
data_merged <- merge.data.table(data_merged, readmap, by.x="qseqid", by.y="contig_id", all.x=T)
data_output <- data_merged[,.SD, .SDcols = !c('contig_length_coverage_file', 'species_MV','count')]

data_output[, contig_length := NULL]
setnames(data_output, "length", "hit_length")
setcolorder(data_output, c("qseqid", "qlen", "sseqid", "slen", "approx_pident", "hit_length", "evalue", "bitscore", "bit_corrected", "average_coverage", "sscinames", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom", "n_staxids", "n_hits", "n_reads", "coverage", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "stitle", "staxids"))


#Create the final annotation file with the merged annotation data
fwrite(data_output, snakemake@output[[1]], sep="\t", quote=F)
