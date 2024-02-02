



###Packages####################33

##install/load packages
pacman::p_load(msa,genepop,seqinr,Biostrings,tidyverse)


####from msa document##############################
txt.shade <- system.file("tex", "texshade.sty", package="msa")
#[1] "C:/Users/trevo/AppData/Local/R/win-library/4.2/msa/tex/texshade.sty"


mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences


my.alignment <- msa(mySequences)
my.alignment



#######ammonia oxidizing microbe msa##############33


getwd()
aom.seqs <-readDNAStringSet("/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/aom.msa.alignment/amoA.msa.fasta",
                            format = "fasta")


aom.alignment <- msa(aom.seqs, method = "ClustalW")


print(aom.alignment, show="complete")#13 gaps in alignment

length(aom.alignment)
names(aom.alignment)

###determining frequency of base pairs

bp.freq <- alphabetFrequency(aom.alignment)

g = bp.freq[,3]
c = bp.freq[,2]


mean(g)#142
mean(c)#126

GC = (142+126)/644
GC#0.4161491
