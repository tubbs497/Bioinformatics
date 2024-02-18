###Packages####################33

##install/load packages
pacman::p_load(msa,genepop,seqinr,Biostrings,tidyverse,phangorn,UniprotR)



#Import a whole set of sequences in one fasta file and name the object dna.seqs
#current directory starts at R project folder


dna.seqs <-readDNAStringSet("input/sequences.fasta",
                            format = "fasta")

#align the sequences in the object dna.seqs with the clustalW method

alignment1 <- msa(dna.seqs, method = "ClustalW")

#pring entire alignement
print(alignment1, show = "complete") #Homo_sapiens_6 [20] has a gap and mismatch

#get alphabet frequency, including mutation notations

alFreq <- alphabetFrequency(alignment1) 
#[20] (HS6) contains one gap and has mismatchess with the rest of the sequences
#[20] A= 141 G= 184   The rest A=144 G=182



#Calcuate distance matrix between the aligned sequences
aln1 <- msaConvert(alignment1, type="seqinr::alignment")
d <- dist.alignment(aln1, "identity")
d# homo_sapien_6 is most different from the rest



#take a single sequences (non-mutation) and make it an object
hs1 <- dna.seqs[1]


##exporting fasta file for blast search
write.fasta(hs1, names = names(hs1), file.out = "output/hs1.fasta")

#Blast Results
#Homo sapiens hbb gene for beta globin, partial cds
#accession LC121775.1



#take the mutation  nt sequence and make it an object
hs6 <- dna.seqs[6]

#Translate the mutation sequence into the amino acid sequence with biostrings

hs6.aa <- Biostrings::translate(hs6)

print(hs6.aa)
#export amino acid sequence to fasta
write.fasta(hs6.aa, names = names(hs6.aa), file.out = "output/hs6_aa.fasta")

#Uniprot blast search results
#Acession A0A0J9YWK4    A0A0J9YWK4_HUMAN	
#Protien name Hemoglobin subunit beta
#Gene name HBB  # length 55 a.a


###create a character object with the accession  number to use GO functions

aa1 <- ("A0A0J9YWK4")
#Create GO object to get pathology/disease information
aa1.go = GetProteinGOInfo(aa1)

aa.1.path = GetPathology_Biotech(aa1.go)
#1] "Bad request. The resource you requested doesn't exist or There is a problem with your input."

#DISEASES#####################

#using the uniprot database, it was determined that the mutation of the second amino acid (Valine) to alanine would cause the disease of beta Thalassemia. If the V was converted to a Leucine it would cause Hb niigata.

print(hs6.aa)
#AAStringSet object of length 1:
#width seq                                                                             names               
#[1]   213 NLLPGAGRAGARAGHESQGRAIYCLHLLLTQLCSLATS...SCTVTSCTWILRTSG*VYGTLDVFFPLLFYG*VHVIGR Homo_sapiens_6

#the second amino acid is L so this variant of the gene may cause Hb niigata potentially causing diabetes.







