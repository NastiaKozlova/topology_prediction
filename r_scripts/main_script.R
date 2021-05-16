part_start<-'/home/nastia/projects/topology_prediction/'
#Put protein sequence into start/sequence/your_protein_name(),
#sequenses to wich you want to compare put into start/alignment/your_protein_name
#form Uniprot
#predicted topologies put into start/topology/protein_name.
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_data.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_consetvative_aminoasids.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/topology_prediction.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plot.R ",part_start),ignore.stdout=T,wait = T) 
