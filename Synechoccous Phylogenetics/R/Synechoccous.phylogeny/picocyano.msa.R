




#Packages########

pacman::p_load(msa,seqinr,Biostrings,tidyverse,phangorn)




###Loading in 16S sequneces of picocyanobacteria from 16S amplicon sequencing####### 

#No type strain
pico.seqs <- Biostrings::readDNAStringSet("/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/Synechoccous Phylogenetics/R/Synechoccous.phylogeny/picocyano.rep.seqs_q2.v3.txt",
                                          format = "fasta",
                                          use.names = TRUE)


#with type strain
pico.seqs.type <- Biostrings::readDNAStringSet("/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/Synechoccous Phylogenetics/R/Synechoccous.phylogeny/picocyano.rep.seqs_q2.type_trim.txt",
                                          format = "fasta",
                                          use.names = TRUE)


#no type strains from blast 16S seq
a1 <- msa(pico.seqs, method = "Muscle")



print(a1, show="complete")




###With type strains


a1.type <- msa(pico.seqs.type, method = "Muscle")

print(a1.type, show="complete")#would need to trim if I use this...
###convert to seqinr format


a1.sr = msaConvert(a1,"seqinr::alignment")

#calculate identity disatnce
a1.dist = dist.alignment(a1.sr,"identity")
ai.dist.mat = as.matrix(a1.dist)#out groups highest percent



##convert to fasta file
A1_Dat = msaConvert(a1, type="phangorn::phyDat")
write.phyDat(A1_Dat, "picocyano.msa.v3.1.phy", format = "phylip")



##convert to fasta file..with type
A1_Dat.t = msaConvert(a1.type, type="phangorn::phyDat")
write.phyDat(A1_Dat.t, "picocyano.msa.v3.type.phy", format = "phylip")


