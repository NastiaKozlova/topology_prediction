#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
part<-(paste0(part_start,"prediction"))

setwd(part)
library(bio3d)
library(dplyr)
library(ggplot2)
library(Peptides)
library(cowplot)
threshold<-0.8
v_palette<-c("active center"="#fcfed6",
             "N-domain"="#c6fde2","C-domain"="#c3cafd","EMD4"="#FFC0CB",
             "topology"="#008F7A",
             "conservative"="#13E600",
             "hydrofobic"="#0F1DFA",
             "charge"="#FF9671")
if(!dir.exists(paste0( "plot/")))    { dir.create(paste0( "plot/"))}
if(!dir.exists(paste0( "plot/topology/")))    { dir.create(paste0( "plot/topology/"))}

if(!dir.exists(paste0( "prediction/")))    { dir.create(paste0( "prediction/"))}
if(!dir.exists(paste0( "prediction/topology_factors/")))  {dir.create(paste0("prediction/topology_factors/"))}
if(!dir.exists(paste0( "prediction/topology_topology/"))) {dir.create(paste0("prediction/topology_topology/"))}

if(!dir.exists(paste0( "topology/")))    { dir.create(paste0( "topology/"))}
if(!dir.exists(paste0( "topology/topology_topology/")))    { dir.create(paste0( "topology/topology_topology"))}
if(!dir.exists(paste0( "topology/topology_factors")))    { dir.create(paste0( "topology/topology_factors"))}
if(!dir.exists(paste0( "topology/topology_full/"))) {dir.create(paste0("topology/topology_full/"))}
srart<-list.files(path = paste0(part_start,"start/topology/"))
start<-NULL
for (i in 1:length(srart)) {
  m<-strsplit( srart[i],split=".", fixed=T)[[1]][1]
  start<-c(start,m)
}
rm(srart)
i<-1
v_tmd_len<-c(13:20)

for (i in 1:length(start)) {
  if(!dir.exists(paste0( "prediction/topology_factors/",start[i])))  {dir.create(paste0("prediction/topology_factors/",start[i]))}
  if(!dir.exists(paste0( "prediction/topology_topology/",start[i]))) {dir.create(paste0("prediction/topology_topology/",start[i]))}
  
  if(!dir.exists(paste0( "topology/topology_factors/",start[i])))  {dir.create(paste0("topology/topology_factors/",start[i]))}
  if(!dir.exists(paste0( "topology/topology_topology/",start[i]))) {dir.create(paste0("topology/topology_topology/",start[i]))}
  
  if(!dir.exists(paste0( "topology/topology_full/",start[i]))) {dir.create(paste0("topology/topology_full/",start[i]))}
  df_topology<-read.csv(paste0("charge/pH_7.4/",start[i],".csv"),stringsAsFactors = F)
  df_topology<-df_topology%>%filter(type=="Transmembrane")

  #find conservative amino acids SLC34A+flouder
  df_conservative<-read.csv(paste0("conservative/",start[i],".csv"),stringsAsFactors = F)
  #extract NaPi2b seqvense 
  seqv<-read.fasta(paste0(part_start,"start/sequence/",start[i],".fasta"))
  v_NaPi2b<-as.character(seqv$ali)
  df_type<-data.frame(matrix(ncol=2,nrow = length(v_NaPi2b)))
  colnames(df_type)<-c("amino","TMD")
  df_type$amino<-c(1:nrow(df_type))
  df_type$TMD<-0 

  for (j in 1:nrow(df_topology)) {
    df_type$TMD[df_type$amino>=df_topology$seq_beg[j]&df_type$amino<=df_topology$seq_end[j]]<-
      df_type$TMD[df_type$amino>=df_topology$seq_beg[j]&df_type$amino<=df_topology$seq_end[j]]+1
  }
  df_type<-df_type%>%mutate(TMD=TMD/length(unique(df_topology$program))*100)
  tmd_len<-v_tmd_len[1]
  for (tmd_len in v_tmd_len) {
    df_hydrophobicity<-read.csv(paste0("hydrophobicity/lenght_",tmd_len,"/",start[i],".csv"),stringsAsFactors = F)
    df_hydrophobicity<-df_hydrophobicity%>%mutate(conservative=NA)
    for (p in 1:nrow(df_hydrophobicity)) {
      df_hydrophobicity$conservative[p]<-mean(df_conservative$conservative[df_hydrophobicity$seq_beg[p]:df_hydrophobicity$seq_end[p]])
    }
#    df_combine<-full_join(df_conservative,df_hydrophobicity,by=c("number"="charge_x"))
#    df_combine<-full_join(df_combine,df_type,by=c("number"="amino"))
    df_hydrophobicity<-df_hydrophobicity%>%mutate(number=charge_x)
    df_hydrophobicity$charge_x<-NULL
    df_combine<-full_join(df_hydrophobicity,df_type,by=c("number"="amino"))
    df_combine<-df_combine%>%filter(!is.na(charge))
#    df_combine<-df_combine%>%mutate(conservative_norn=(conservative-min(conservative)))
#    df_combine<-df_combine%>%mutate(conservative_norn=(conservative_norn/max(conservative_norn)*100))

    df_combine<-df_combine%>%mutate(conservative_norn=conservative)
    df_combine<-df_combine%>%mutate(hydrofobic_norn=(-hydrofobic)) 
    df_combine<-df_combine%>%mutate(hydrofobic_norn=(hydrofobic_norn-min(hydrofobic_norn)))
    df_combine<-df_combine%>%mutate(hydrofobic_norn=(hydrofobic_norn/max(hydrofobic_norn)*100))

    df_combine<-df_combine%>%mutate(charge_norn=abs(charge))
    df_combine<-df_combine%>%mutate(charge_norn=(charge_norn-min(charge_norn)))
    df_combine<-df_combine%>%mutate(charge_norn=(1-charge_norn/max(charge_norn))*100)
    df_combine<-df_combine%>%mutate(topology_norm=TMD)
    df_combine<-df_combine%>%mutate(tmd_cons=(conservative_norn+ hydrofobic_norn+charge_norn)/3)
#    df_combine<-df_combine%>%mutate(tmd_cons=(topology_norm+ conservative_norn+ hydrofobic_norn+charge_norn)/4)
    write.csv(df_combine,paste0( "topology/topology_full/",start[i],"/",tmd_len,".csv"),row.names = F)
    topology_tresh<-quantile(df_combine$topology_norm,probs = threshold,na.rm = T)
    tmd_tresh<-quantile(df_combine$tmd_cons,probs = threshold,na.rm = T)
    df_combine_topology<-df_combine%>%filter(topology_norm>topology_tresh)

    df_combine_factors<-df_combine%>%filter(tmd_cons>tmd_tresh)

    
    if(nrow(df_combine_topology)>0){
    for (q in 1:(nrow(df_combine_topology)-1)) {
      if (df_combine_topology$seq_beg[q+1]<df_combine_topology$seq_end[q]) {
        df_combine_topology$seq_beg[q+1]<-df_combine_topology$seq_beg[q]
        df_combine_topology$seq_beg[q]<-NA
      }
    }
    df_combine_topology<-df_combine_topology%>%filter(!is.na(seq_beg))
    df_combine_topology<-df_combine_topology%>%mutate(leng=seq_end-seq_beg)

    write.csv(df_combine_topology,paste0("topology/topology_topology/",start[i],"/",tmd_len,".csv"),row.names = F)
    }
    if(nrow(df_combine_factors)>0){
    for (q in 1:(nrow(df_combine_factors)-1)) {
      if (df_combine_factors$seq_beg[q+1]<df_combine_factors$seq_end[q]) {
        df_combine_factors$seq_beg[q+1]<-df_combine_factors$seq_beg[q]
        df_combine_factors$seq_beg[q]<-NA
      }
    }
    df_combine_factors<-df_combine_factors%>%filter(!is.na(seq_beg))
    df_combine_factors<-df_combine_factors%>%mutate(leng=seq_end-seq_beg)
    write.csv(df_combine_factors,paste0("topology/topology_factors/",start[i],"/",tmd_len,".csv"),row.names = F)
    
    df_combine_topology<-df_combine%>%filter(topology_norm>topology_tresh)
    df_combine_factors<-df_combine%>%filter(tmd_cons>tmd_tresh)
    df_combine_topology<-df_combine_topology%>%mutate(seq_beg=number-20)
    df_combine_topology<-df_combine_topology%>%mutate(seq_end=number+20)
    df_combine_factors<-df_combine_topology%>%mutate(seq_beg=number-20)
    df_combine_factors<-df_combine_topology%>%mutate(seq_end=number+20)
    write.csv(df_combine_factors,paste0("prediction/topology_factors/",start[i],"/",tmd_len,".csv"),row.names = F)
    write.csv(df_combine_topology,paste0("prediction/topology_topology/",start[i],"/",tmd_len,".csv"),row.names = F)
    }
  }
}
