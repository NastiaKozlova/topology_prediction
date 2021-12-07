#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(Peptides)
library(ggplot2)

part<-paste0(part_start,"prediction/")
if (!dir.exists(part)){dir.create(part)}
v_hydro<-c("Casari", "Engelman")
#don't use thise topologeis
not_topologies<-c("merge","low_pH")
setwd(part)
seq_leng<-c(13:20)
pKscale_choose<-"Stryer"

v_pH<-c(seq(from=5,to= 7.5,by =0.5),5.8,7.4,c(8:14))

# name_interaction may be "strict interactions" (for aminoasids wicth interaction witch one aminoasids) 
# and "non-strict interactions" (for grupps of aminoasids wicth interaction witch every aminoasids in this grupp)  

srart<-list.files(path = paste0(part_start,"start/topology/"))
start<-NULL
for (i in 1:length(srart)) {
  m<-strsplit( srart[i],split=".", fixed=T)[[1]][1]
  start<-c(start,m)
}
rm(srart)
df_charactirise<-data.frame(matrix(nrow = length(start),ncol =9))

colnames(df_charactirise)<-c("pdbID","MD_num","TMD_num","chainID","seq","N_type","C_type","N_leng","C_leng")
df_charactirise$pdbID<-start
i<-1
if(!dir.exists(paste0( "charge")))    { dir.create(paste0( "charge"))}
i<-1

for (pH in v_pH) {
  if(!dir.exists(paste0( "charge/pH_",pH)))    { dir.create(paste0("charge/pH_",pH))}
}

#convert topology to different format
for (pH in v_pH) {
  for (i in 1:nrow(df_charactirise)) {
    #read seq
    seq<-read.fasta(paste0(part_start,"start/sequence/",(df_charactirise$pdbID[i]),".fasta"))
    v_seq<-as.vector(seq$ali)
    
    v_topology<-list.files(paste0(part_start,"start/topology/",(df_charactirise$pdbID[i]),"/seq_0/"))
    a<-c("OCTOPUS", "philius","PolyPhobius",  "SCAMPI_MSA",     "SPOCTOPUS","Topcons")
    v_topology<-a[a%in%v_topology]
    
    df_topology<-data.frame(matrix(ncol=4,nrow=length(v_seq)))
    colnames(df_topology)<-c("resno","resid","topology","program")
    df_topology$resid<-v_seq
    df_TEMP<-read.table(paste0(part_start,"start/topology/",df_charactirise$pdbID[i],"/seq_0/Topcons/topcons.top"),header = F)
    v_test<-strsplit(df_TEMP[1,1],split = "",fixed = T)[[1]]
    df_topology$topology<-v_test
    df_topology$program<-"Topcons"
    df_topology<-df_topology%>%mutate(resno=1:nrow(df_topology))
    v_topology<-v_topology[v_topology!="Topcons"]
    df_topology<-df_topology%>%mutate(seq_beg=resno)
    df_topology<-df_topology%>%mutate(seq_end=resno)
    j<-2
    for (j in 2:nrow(df_topology)) {
      if(df_topology$topology[j-1]==df_topology$topology[j]){
        df_topology$seq_beg[j]<-df_topology$seq_beg[j-1]
        df_topology$resno[j-1]<-NA
      }
    }
    df_topology<-df_topology%>%filter(!is.na(resno))
    for (p in 1:length(v_topology)) {
      
      df_topology_add<-data.frame(matrix(ncol=4,nrow=length(v_seq)))
      colnames(df_topology_add)<-c("resno","resid","topology","program")
      df_topology_add$resid<-v_seq
      df_TEMP<-read.table(paste0(part_start,"start/topology/",df_charactirise$pdbID[i],"/seq_0/Topcons/topcons.top"),header = F)
      v_test<-strsplit(df_TEMP[1,1],split = "",fixed = T)[[1]]
      df_topology_add$topology<-v_test
      df_topology_add$program<-v_topology[p]
      df_topology_add<-df_topology_add%>%mutate(resno=1:nrow(df_topology_add))
      df_topology_add<-df_topology_add%>%mutate(seq_beg=resno)
      df_topology_add<-df_topology_add%>%mutate(seq_end=resno)
      
      for (j in 2:nrow(df_topology_add)) {
        if(df_topology_add$topology[j-1]==df_topology_add$topology[j]){
          df_topology_add$seq_beg[j]<-df_topology_add$seq_beg[j-1]
          df_topology_add$resno[j-1]<-NA
        }
      }
      df_topology<-rbind(df_topology,df_topology_add)
    }
    df_topology<-df_topology%>%filter(!is.na(resno))
    df_topology$type[df_topology$topology=="i"]<-"Cytoplasmic"
    df_topology$type[df_topology$topology=="M"]<-"Transmembrane"
    df_topology$type[df_topology$topology=="o"]<-"Extracellular"
    df_topology<-df_topology%>%select(seq_beg,seq_end,type,program)
    df_topology<- df_topology%>%mutate(seq=NA)
    for (q in 1:nrow(df_topology)) {
      df_topology$seq[q]<-paste0(v_seq[df_topology$seq_beg[q]:df_topology$seq_end[q]],collapse = "")
    }
    
    df_topology<- df_topology%>%mutate(charge=round(charge(seq, pH = pH, pKscale = "Stryer") ,digits = 1))
    df_topology<- df_topology%>%mutate(hydrofobic=round(hydrophobicity(seq, scale = "Engelman"),digits = 1))
    
    df_topology$topology_number[df_topology$type=="Cytoplasmic"]<-1
    df_topology$topology_number[df_topology$type=="Transmembrane"]<-2
    df_topology$topology_number[df_topology$type=="Extracellular"]<-3
    df_topology$topology_number[df_topology$type=="Periplasmic"]<-3
    df_topology<-df_topology%>%mutate(charge_x=(seq_beg+seq_end)/2)
    df_topology<-df_topology%>%mutate(charge_y=topology_number+0.25)
    df_topology<-df_topology%>%mutate(y_hidro=topology_number-0.25)
    #sort not used topologies
    df_topology<-df_topology[!df_topology$program%in%not_topologies,]
    write.csv(df_topology,paste0("charge/pH_",pH,"/",df_charactirise$pdbID[i],".csv"),row.names = F)
  }
}
#calculate hydrophobicity and charge of fragments with length from 13 to 20
pH<-7.4
if(!dir.exists(paste0( "hydrophobicity")))    { dir.create(paste0( "hydrophobicity"))}
for (seq_len in seq_leng) {
  if(!dir.exists(paste0( "hydrophobicity/lenght_",seq_len)))    { dir.create(paste0( "hydrophobicity/lenght_",seq_len))}
  for (i in 1:nrow(df_charactirise)) {
    seq<-read.fasta(paste0(part_start,"start/sequence/",(df_charactirise$pdbID[i]),".fasta"))
    v_seq<-as.vector(seq$ali)
    df_topology<-data.frame(matrix(ncol=2,nrow = (length(v_seq)-seq_len+1)))
    colnames(df_topology)<-c("seq_beg","seq_end")
    df_topology$seq_beg<-1:nrow(df_topology)
    df_topology$seq_end<-c(seq_len:(nrow(df_topology)+seq_len-1))
    df_topology<-df_topology%>%mutate(seq=NA)
    for (j in 1:nrow(df_topology)) {
      df_topology$seq[j]<-paste0(v_seq[df_topology$seq_beg[j]:df_topology$seq_end[j]],collapse = "")
    }
    df_topology<- df_topology%>%mutate(charge=round(charge(seq, pH = pH, pKscale = "Stryer") ,digits = 1))
    df_topology<- df_topology%>%mutate(hydrofobic=round(hydrophobicity(seq, scale = "Engelman"),digits = 1))
    df_topology<-df_topology%>%mutate(charge_x=(seq_beg+seq_end)/2)
    write.csv(df_topology,paste0("hydrophobicity/lenght_",seq_len,"/",df_charactirise$pdbID[i],".csv"),row.names = F)
  }
}

if(!dir.exists(paste0( "plot")))    { dir.create(paste0( "plot"))}
if(!dir.exists(paste0( "plot/charge")))    { dir.create(paste0( "plot/charge"))}
x<-seq(from=0,to=5000,by=20)
#make plot protein topology with charges 
for (i in 1:nrow(df_charactirise)) {
  df_topology<-read.csv(paste0("charge/pH_7.4/",df_charactirise$pdbID[i],".csv"),stringsAsFactors =  F)
  df_topology<-df_topology%>%filter(program!="merge")
  
  df_topology<-df_topology%>%mutate(charge_y=topology_number)
  df_topology<-df_topology%>%mutate(topology_end=NA)
  df_topology$topology_end[df_topology$topology_number==1]<-1.05
  df_topology$topology_end[df_topology$topology_number==3]<-2.95
  
  df_topology$topology_end[df_topology$topology_number==2]<-3
  df_topology$topology_number[df_topology$topology_number==2]<-1.
  df_topology_ref<-df_topology
  df_topology<-df_topology_ref%>%filter(program!="low_pH")
  p0<-ggplot(data = df_topology)+
    labs(x="",y="",title =df_charactirise$pdbID[i] )+
    geom_rect(aes(xmin=seq_beg,xmax=seq_end,ymin=topology_number,ymax=topology_end,colour=type,fill=type))+
    geom_text(aes(x=charge_x,y=charge_y,label=charge))+
    scale_y_continuous(breaks = NULL,labels = NULL,limits = c(0.5,3.5))+
    scale_x_continuous(breaks = x,labels = x,limits=c(min(df_topology$seq_beg),max(df_topology$seq_end)))+
    facet_grid(program ~ .) + theme_bw()+ 
    guides(color = "none") + guides(fill = "none")
  
  ggsave(plot = p0, filename = paste0("plot/charge/",df_charactirise$pdbID[i],".png"),  width = 22, height = 30, units = c("cm"), dpi = 300 )
}
