require(rentrez)
require(ape)


#First we search for a group - here, a bunch of fungi, but you choose yours as appropriate

download_fasta = function(gene,taxon){
  search_term = paste(gene,"[gene] AND 100:1000[SLEN] AND ",taxon,"[ORGN]",sep="")
  r_search = entrez_search(db="nucleotide",term=search_term,use_history = TRUE)
  #We then process the returned results. First we download a summary of each record
  r_summary = entrez_summary(db="nucleotide",web_history = r_search$web_history,retmax=500)
  #Next we ask for what taxon each record belongs to
  taxids = extract_from_esummary(r_summary,"taxid")
  #Finally, we get a list of one sample per taxon - this makes life a little easier
  uids = extract_from_esummary(r_summary,"uid")[match(unique(taxids),taxids)]
  r_fastas = entrez_fetch(db="nucleotide",id=uids,rettype="fasta")
  file_out = paste("data/",taxon,"_",gene,".fasta",sep="")
  write(r_fastas,file=file_out)
}

system("mafft ")

system("iqtree2")