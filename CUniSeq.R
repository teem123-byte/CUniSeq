library(seqinr)
library(seqRFLP)
library(Biostrings)
library(msa)
library(DECIPHER)
library(plyr)

genes<-readDNAStringSet('/Users/hamiji/Downloads/targetDB_nt_cons.fasta')
genes_compare<-genes
CS<-data.frame()
for(j in 1:length(genes_compare)){
  gene<-1
  gene_2<-2
  for(i in 1:length(genes_compare)){
    if(gene_2<=length(genes_compare)){
      comparison<-paste(names(genes_compare[gene]), "/", names(genes_compare[gene_2]))
      MA<-msa(genes_compare[c(gene,gene_2)], method = "ClustalOmega")
      MA<-DNAStringSet(MA)
      CS_1<-ConsensusSequence(MA, threshold = 0.05)
      CS_1<-unlist(strsplit(as.character(CS_1),""))
      NWCBP<-c("R","W" ,"Y", "K", "D","S" ,"H" ,"M" ,"B", "V")
      CS_1<-replace(CS_1, which(CS_1 %in% NWCBP), "+")
      CS_1<-paste(CS_1, sep = "", collapse = "")
      CS_1<-cbind(comparison,CS_1)
      CS<-rbind(CS, CS_1)
    }
    gene_2<-gene_2+1
  }
  genes_compare<-genes_compare[-1]
}

#scanning function
UniqueConservedSequence<-function(CS, DistanceThreshold, MinSeqLength){
  seq<-list()
  for(i in 1:length(CS)){
    seqcnt<-1
    startingposition<-i
    currentposition<-i+1
    distance<-1
    for(j in 1:(length(CS)-currentposition)){
      if(CS[currentposition] %in% "+"){
        distance<-1
        currentposition<-currentposition+1
      }
      if(distance<=DistanceThreshold){
        distance<-distance+1
        currentposition<-currentposition+1
      }
      seq1<-paste(CS[startingposition:currentposition], collapse = "")
      seqlength<-nchar(seq1)
      seq1<-cbind(seq1, seqlength, startingposition, currentposition)
      seq[[seqcnt]]<-seq1
      print(seq[[seqcnt]])
      if(seqcnt !=1){
        currentpos_prev<-as.numeric(seq[[seqcnt-1]][4])
      }
      if(seqcnt==1){
        currentpos_prev<-0
      }
      if(currentpos_prev==currentposition) {break}
      seqcnt<-seqcnt+1
    }
  }
  seqf <- data.frame(matrix(unlist(seq), nrow=length(seq), byrow=T))
  colnames(seqf)<-c("sequence","sequence length", "starting position", "ending position")
  
  seqf<-seqf[which(as.numeric(seqf$`sequence length`) >= MinSeqLength),]
  seqf<-seqf[order(seqf$`sequence length`, decreasing = T),]
  return(seqf)
} 

#starting and ending positiions are wrong

teh<-CS

output_list<-list()
for(i in 1:nrow(teh)){
  input<-unlist(strsplit(teh$CS_1[i],""))
  v<-UniqueConservedSequence(input,10,100)
  output_list[[i]]<-v
}






