part_start<-'/home/nastia/mem/NaPi2b/topology/'
#put sequence of your protein into start folder as start/target.fasta

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_data.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_consetvative_aminoasids.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/topology_prediction.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plot.R ",part_start),ignore.stdout=T,wait = T) 
