#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
setwd(paste0(part_start,"prediction"))
library(bio3d)
library(dplyr)
library(ggplot2)

srart<-list.files(path = paste0(part_start,"start/topology/"))
start<-NULL
for (i in 1:length(srart)) {
  m<-strsplit( srart[i],split=".", fixed=T)[[1]][1]
  start<-c(start,m)
}
rm(srart)
seqname<-start[1]
if(!dir.exists(paste0( "conservative")))    { dir.create(paste0( "conservative"))}
for(seqname in start){
  seq<-read.fasta(paste0(part_start,"start/sequence/",seqname,".fasta"))
  aingment<-read.fasta(paste0(part_start,"start/alignment/",seqname,".fasta"))
  aingment<-seqaln(aingment)
  v_ali<-aingment$ali
  v_id<-aingment$id
  v_seqid<-seq$id
  num_id<-match(v_seqid,v_id)
  if(is.na(num_id)){print(paste0("wrong names of seqences\n",
                                 "check sequence id in start/sequence/",seqname,".fasta and ","start/alignment/",seqname,".fasta"))}
  v_ali<-t(v_ali)
  df_ali<-data.frame(v_ali)
  df_ali<-data.frame(lapply(df_ali, as.character), stringsAsFactors=FALSE)
#  colnames(df_ali)<-v_ids
  df_ali[df_ali=="-"]<-NA
  df_ali<-df_ali[!is.na(df_ali[,num_id]),]
  v_per<-ncol(df_ali)
  df_ali<-df_ali%>%mutate(number=1:nrow(df_ali))
  df_ali<-df_ali%>%mutate(percent=0)
#  df_ali
  for (i in 1:nrow(df_ali)) {
    a<-as.vector(t(df_ali[i,]))
    df_ali$percent[i]<-length(a[a%in%a[num_id]])/v_per*100
  }
  df_ali<-df_ali%>%select(number,colnames(df_ali)[num_id],percent)
  colnames(df_ali)<-c("number","amino","conservative")
  write.csv(df_ali,paste0("conservative/",seqname,".csv"),row.names = F)
}
