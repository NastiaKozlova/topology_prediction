# topology_prediction

Predict topology of transmembrane protein based on predicted topologies and other factors such as charge, hydrophobisity and conservatism.
Change part\_start to your working directory name
# Input
Put protein sequence into _*start/sequence/protein_name*_, sequenses to wich you
want to compare put into _*start/alignment/protein_name*_, predicted topologies put into 
_*start/topology/protein_name*_.

# Output 

charge and hydrophobicity of protein domains predicted by many programs are in prediction/charge
charge and hydrophobicity of protein pieses of sequence are in prediction/hydrophobicity
percentage of amino acids conservatism are output in prediction/conservative

# List of dependenses

Programs:

1. Muscle
2. R

R pacages:

1. dplyr (install.packages("dplyr"))
2. ggplot2 (install.packages("ggplot2"))
3. cowplot (install.packages("cowplot"))
4. Peptides (install.packages("Peptides"))