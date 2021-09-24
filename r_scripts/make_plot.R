#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
part<-(paste0(part_start,"prediction"))
setwd(part)
library(bio3d)
library(dplyr)
library(ggplot2)
library(Peptides)
library(cowplot)
v_palette<-c("active center"="#fcfed6",
             "N-domain"="#c6fde2","C-domain"="#c3cafd","EMD4"="#FFC0CB",
             "topology"="#008F7A",
             "conservative"="#13E600",
             "hydrofobic"="#0F1DFA",
             "charge"="#FF9671",
             "QSSS"="#0000FF")
if(!dir.exists(paste0( "plot/")))    { dir.create(paste0( "plot/"))}
if(!dir.exists(paste0( "plot/topology/")))    { dir.create(paste0( "plot/topology/"))}
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
i<-1
tmd_len<-13
pH<-7.4
for (i in 1:length(start)) {
  seq<-read.fasta(paste0(part_start,"start/sequence/",start[i],".fasta"))
  if(file.exists(paste0(part_start,"start/mutations/",start[i],".csv"))){
    df_mutations<-read.csv(paste0(part_start,"start/mutations/",start[i],".csv"),stringsAsFactors = F)
    df_conservative<-read.csv(paste0(part_start,"prediction/conservative/",start[i],".csv"),stringsAsFactors = F)
    df_conservative<-full_join(df_conservative,df_mutations,by=c("number"="Amino_Acid"))
    df_conservative<-df_conservative%>%filter(!is.na(Sites))
    
#    df_conservative<-df_conservative%>%filter(conservative>90)
  }
  v_seq<-as.vector(seq$ali)
  protein_leng<-length(v_seq)
  for (tmd_len in v_tmd_len) {
    a<-file.exists(paste0( "topology/topology_full/",start[i],"/",tmd_len,".csv"))
    b<-file.exists(paste0( "topology/topology_topology/",start[i],"/",tmd_len,".csv"))
    c<-file.exists(paste0( "topology/topology_factors/",start[i],"/",tmd_len,".csv"))
    if (a&b&c){
      df_combine<-read.csv(paste0( "topology/topology_full/",start[i],"/",tmd_len,".csv"),stringsAsFactors = F)
      df_combine_topology<-read.csv(paste0("topology/topology_topology/",start[i],"/",tmd_len,".csv"),stringsAsFactors = F)
      df_combine_factors<-read.csv(paste0("topology/topology_factors/",start[i],"/",tmd_len,".csv"),stringsAsFactors = F)
      df_topology<-read.csv(paste0("charge/pH_",pH,"/",start[i],".csv"),stringsAsFactors = F)
      df_topology<-df_topology%>%filter(type=="Transmembrane")
      df_QSSS<-df_combine%>%mutate(motif=NA)
      for (j in 1:nrow(df_QSSS)) {
        df_QSSS$motif[j]<-gregexpr(df_QSSS$seq[j],pattern = "QSSS")[[1]]
      }
      df_QSSS<-df_QSSS%>%filter(motif>0)
      df_QSSS<-df_QSSS%>%mutate(motif=(motif+seq_beg))
      QSSS<-unique(df_QSSS$motif)
      v_mutations<-sort(unique(df_mutations$POSITION))
      p0_cf<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="charge full")+
        
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
       # geom_vline(xintercept=df_conservative$number)+
        geom_line(aes(x=number,y=charge,alpha=10))+
        geom_segment(aes(x=Inf,xend=-Inf,y=0,yend=0))+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        geom_vline(xintercept = QSSS)+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      
      p0_all<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="TMD nornalised")+
        
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_line(aes(x=number,y=tmd_cons,alpha=10,colour="hydrofobic"))+
        geom_line(aes(x=number,y=topology_norm,alpha=10,colour="topology"))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      
      p0_ch<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="charge")+
       # geom_vline(xintercept=df_conservative$number)+
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_line(aes(x=number,y=charge_norn,alpha=10))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      
      p0_h<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="hydrofobic")+
        
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_line(aes(x=number,y=hydrofobic_norn,alpha=10))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      
      p0_c<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="conservative")+
       # geom_vline(xintercept=df_conservative$number)+
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_line(aes(x=number,y=conservative_norn,alpha=10))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      
      p0_t<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="topology")+
       # geom_vline(xintercept=df_conservative$number)+
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_line(aes(x=number,y=topology_norm,alpha=10))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette,name="normalised parameters")+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+
        guides(alpha = "none")+
        guides(fill = "none")+
        #geom_segment(aes(x=seq_beg,xend=seq_end,y=0,yend=0),data=df_topology_Uniprot)+ 
        #geom_segment(aes(x=seq_beg,xend=seq_end,y=0,yend=0),data=df_topology_Uniprot)+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      p1<-ggplot(data = df_combine)+
        labs(x="Number amino acid",y="")+
       # geom_vline(xintercept=df_conservative$number)+
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_segment(aes(x=seq_beg,xend=seq_end,y="Topology\nfactors",yend="Topology\nfactors",alpha=10),data =df_combine_factors)+
        geom_segment(aes(x=seq_beg,xend=seq_end,y="Topology\non topology",yend="Topology\non topology", alpha=10),data =df_combine_topology )+
        scale_color_manual(values = v_palette)+scale_fill_manual(values = v_palette,name="Domain name")+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        guides(color = "none")+   guides(alpha = "none")+  #guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        theme_bw()+  theme(legend.position = "bottom")
      p2<-ggplot(data = df_topology)+
        labs(x="Number amino acid",y="")+
       # geom_vline(xintercept=df_conservative$number)+
        geom_rect(aes(xmin = 110, ymin = -Inf, xmax= 249, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 374, ymin = -Inf, xmax= 560, ymax = Inf,fill="active center",color="active center"))+
        geom_rect(aes(xmin = 1, ymin = -Inf, xmax= 100, ymax = Inf,fill="N-domain",color="N-domain"))+
        geom_rect(aes(xmin = 574, ymin = -Inf, xmax= protein_leng, ymax = Inf,fill="C-domain",color="C-domain"))+
        geom_rect(aes(xmin = 250, ymin = -Inf, xmax= 360, ymax = Inf,fill="EMD4",color="EMD4"))+
        geom_vline(xintercept = QSSS)+
       # geom_vline(xintercept=df_conservative$number)+
        scale_color_manual(values = v_palette)+scale_fill_manual(values = v_palette,name="Domain name")+
        guides(color = "none")+guides(alpha = "none")+  guides(fill = "none")+
        scale_x_continuous(limits = c(0,protein_leng),breaks = seq(from=0,to=protein_leng,by=20),labels =seq(from=0,to=protein_leng,by=20) )+
        geom_segment(aes(x=seq_beg,xend=seq_end,y=program,yend=program))+
        theme_bw()+  theme(legend.position = "bottom")
      p<-plot_grid(p0_cf,p0_c,p0_ch,p0_h,p0_t,p2,p1,ncol=1,align = "v",rel_heights = c(1.1,1.1,1.1,1.1,1.1,2.8,1.5,1,1))
      ggsave(plot = p, filename = paste0("plot/topology/",start[i],"_",tmd_len,".png"),  width = 30, height =30, units = c("cm"), dpi = 300)
    }
  }
}
