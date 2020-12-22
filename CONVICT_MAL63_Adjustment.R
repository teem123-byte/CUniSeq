library(openxlsx)

seqs<-read.xlsx("/Users/hamiji/Desktop/CONVICT/MAL63c9 and c2_ncbi.xlsx", colNames = F)
MalR<-readDNAMultipleAlignment('/Users/hamiji/Downloads/200703_MALR DNA alignment_Bergstrom and ncbi.fa')

seqs$X2<-sub(">","",seqs$X2)

c2<-seqs[which(seqs$X1=="Mal63c2"),]
c9<-seqs[which(seqs$X1=="Mal63c9"),]

rownames(c2)<-c2$X2
rownames(c9)<-c9$X2

c2$X2<-sub('.*sequence', '', c2$X2)
c2$X2<-sub('.*genome', '', c2$X2)
c2$X2<-sub('.*)', '', c2$X2)
c2$X2<-sub('.*NBRC', '', c2$X2)
c2$X2<-sub('.*CBS', '', c2$X2)
c2$X2<-sub('.*7', '', c2$X2)
c2$X2<-sub('.*cds', '', c2$X2)
c2$X2<-sub('.*III', '', c2$X2)
c2$X2<-sub('.*protein', '', c2$X2)
c2$X2<-sub('.*ScXI', '', c2$X2)
c2$X2<-sub('.*1948', '', c2$X2)
c2$X2<-sub('.*2003', '', c2$X2)
c2$X2<-sub('.*2031', '', c2$X2)
c2$X2<-sub('.*380', '', c2$X2)
c2$X2<-sub('.*539', '', c2$X2)
c2$X2<-sub('.*11', '', c2$X2)



c9$X2<-sub('.*sequence', '', c9$X2)
c9$X2<-sub('.*genome', '', c9$X2)
c9$X2<-sub('.*)', '', c9$X2)
c9$X2<-sub('.*NBRC', '', c9$X2)
c9$X2<-sub('.*CBS', '', c9$X2)
c9$X2<-sub('.*7', '', c9$X2)
c9$X2<-sub('.*cds', '', c9$X2)
c9$X2<-sub('.*III', '', c9$X2)
c9$X2<-sub('.*protein', '', c9$X2)
c9$X2<-sub('.*ScXI', '', c9$X2)
c9$X2<-sub('.*1948', '', c9$X2)
c9$X2<-sub('.*2003', '', c9$X2)
c9$X2<-sub('.*2031', '', c9$X2)
c9$X2<-sub('.*380', '', c9$X2)
c9$X2<-sub('.*539', '', c9$X2)
c9$X2<-sub('.*11', '', c9$X2)

#Ignore everything before this, I had to preprocess your data because I could 
#not read it in as a DNAstringset later. We should make sure that if we use this in the future
# we should just upload the data is a usable format.
library(seqRFLP)
library(Biostrings)
library(msa)
library(DECIPHER)
library(plyr)

c2<-DNAStringSet(c2$X2)
c9<-DNAStringSet(c9$X2)


c2MA<-msa(c2, method = "ClustalOmega")
c2MA<-DNAStringSet(c2MA)
c9MA<-msa(c9, method = "ClustalOmega")
c9MA<-DNAStringSet(c9MA)



c2cons<-ConsensusSequence(c2MA, threshold = 0.2)
c9cons<-ConsensusSequence(c9MA, threshold = 0.2)

c2c9<-DNAStringSet(c(c2cons,c9cons))

c2c9MA<-DNAStringSet(c2c9)
c2c9cons<-ConsensusSequence(c2c9MA, threshold = 0.01)
c2c9consdar<-unlist(strsplit(as.character(c2c9cons),""))
NWCBP<-c("R","W" ,"Y", "K", "D","S" ,"H" ,"M" ,"B", "V")
c2c9consdar<-replace(c2c9consdar, which(c2c9consdar %in% NWCBP), "+")
c2c9cons<-paste(c2c9consdar, sep = "", collapse = "")




MalR<-DNAStringSet(MalR)
MalR_MA<-ConsensusSequence(MalR, threshold = 0.4)
MalR_MA_consdar<-unlist(strsplit(as.character(MalR_MA),""))
MalR_MA_consdar<-replace(MalR_MA_consdar, which(MalR_MA_consdar %in% NWCBP), "+")
MalR_MA_cons<-paste(MalR_MA_consdar, sep = "", collapse = "")


CS<-DNAStringSet(c(c2c9cons,MalR_MA_cons))
CS<-ConsensusSequence(CS, threshold = 0.01)
CS_consdar<-unlist(strsplit(as.character(CS),""))
CS_consdar<-replace(CS_consdar, which(CS_consdar %in% NWCBP), "+")
CS<-paste(CS_consdar, sep = "", collapse = "")



inputseq<-strsplit(CS, "")[[1]]

#scanning function
UniqueConservedSequence<-function(CS, DistanceThreshold, MinSeqLength){
  sequences<-list()
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
      seq1<-paste(CS[startingposition:currentposition], collapse = " ")
      seqlength<-nchar(seq1)
      seq1<-cbind(seq1, seqlength, startingposition, currentposition)
      seq[[seqcnt]]<-seq1
      seqcnt<-seqcnt+1
    }
  }
  
  seqf <- data.frame(matrix(unlist(seq), nrow=length(seq), byrow=T))
  colnames(seqf)<-c("sequence","sequence length", "starting position", "ending position")
  
  MinSeqLength<-100
  seqf<-seqf[which(as.numeric(seqf$`sequence length`) >= minimumseqlength),]
  seqf<-seqf[order(seqf$`sequence length`, decreasing = T),]
  return(seqf)
} 

#starting and ending positiions are wrong
test<-UniqueConservedSequence(inputseq,5,100)
write.