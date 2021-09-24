part_start<-'path to your directory'
#Put protein sequence into start/sequence/your_protein_name(),
#sequenses to wich you want to compare put into start/alignment/your_protein_name
#form Uniprot
#predicted topologies put into start/topology/protein_name.
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_data.R ",part_start),ignore.stdout=T,wait = T) 
#charge and hydrophobicity of protein domains predicted by many programs are in prediction/charge
#charge and hydrophobicity of protein pieses of sequence are in prediction/hydrophobicity
#plot of charge of protein domains predicted by many programs are in prediction/plot/charge

#alignment are made by muscle
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_consetvative_aminoasids.R ",part_start),ignore.stdout=T,wait = T) 
#percentage of amino acids conservatism are output in prediction/conservative
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/topology_prediction.R ",part_start),ignore.stdout=T,wait = T) 
#topologies based on topologies are in topology/topology_topology
#topologies based on factors (charge, hydrophobicity and conservativity of pieces of protein) are in topology/topology_factors
#full uncut data about pieces are in topology/topology_full

#data to predict structure of protein based on topology are present in folders
#prediction/prediction/topology_factors
#prediction/prediction/topology_topology
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plot.R ",part_start),ignore.stdout=T,wait = T) 
#final topologies plots are in prediction/plot/topology