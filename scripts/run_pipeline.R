#Run tree generation pipeline

#User can edit these:
gene = "COI"
taxon = "simiiformes"

#Set some variables

search_term = paste(gene,"[gene] AND 100:1000[SLEN] AND ",taxon,"[ORGN]",sep="")
file_prefix = paste(gene,taxon,sep="_")
output_fasta = paste("data/",file_prefix,".fasta",sep="")

#Run the code to get a fasta
##First search
r_search = entrez_search(db="nucleotide",term=search_term,use_history = TRUE)
##We then process the returned results. First we download a summary of each record
r_summary = entrez_summary(db="nucleotide",web_history = r_search$web_history,retmax=500)
##Next we ask for what taxon each record belongs to
taxids = extract_from_esummary(r_summary,"taxid")
##Finally, we get a list of one sample per taxon - this makes life a little easier
uids = extract_from_esummary(r_summary,"uid")[match(unique(taxids),taxids)]
##Fetch the fastas
r_fastas = entrez_fetch(db="nucleotide",id=uids,rettype="fasta")
##And write out the fasta

write(r_fastas,file=output_fasta)


## Align
output_aligned = paste("results/",file_prefix,".aligned.fasta",sep="")
mafft_command = paste("mafft --auto ", output_fasta," > ", output_aligned,sep="")

system(mafft_command)

## IQTREE2

iqtree_command = paste("iqtree2 -s ",output_aligned," --prefix results/",file_prefix,sep="")

system(iqtree_command)

## Read in tree

tree_file = paste("results/",file_prefix,".treefile")
tree = read.tree(tree_file)
tree_plot = ggtree(tree)+geom_tiplab()

#Wrte out tree

tree_plot_file = paste("figures/",file_prefix,".pdf")
