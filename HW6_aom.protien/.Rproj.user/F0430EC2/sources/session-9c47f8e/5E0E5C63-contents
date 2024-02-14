###Packages####################33
getwd()#"C:/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/HW6_aom.protien"

##install/load packages
pacman::p_load(UniprotR,protti,GenomicAlignments,seqinr,Biostrings,tidyverse,phangorn,r3dmol)


#amino acid sequence to fasta in previous excercise


#importing protein accession numbers


id <- read.csv("input/uniprot_accession.csv")


#creating character strings for each individual accession number

aa1 <- id$accession[1]
aa2 <- id$accession[2]
aa3 <- id$accession[3]
aa4 <- id$accession[4]
aa5 <- id$accession[5]

#obtaining gene ontology: aa1
aa1.go = GetProteinGOInfo(aa1)

PlotGoInfo(aa1.go)
#getting genes functional role
aa1.go$Gene.Ontology..molecular.function.
#[1] "ammonia monooxygenase activity [GO:0018597]"

aa1.go$Gene.Ontology..GO. #1] "membrane [GO:0016020]; ammonia monooxygenase activity [GO:0018597]"

#obtaining gene ontology: aa2
aa2.go = GetProteinGOInfo(aa2)

PlotGoInfo(aa2.go)

#obtaining gene ontology: aa3
aa3.go = GetProteinGOInfo(aa3)

PlotGoInfo(aa3.go)
#Error in `colnames<-`(`*tmp*`, value = c("Goterm", "Count")) : 
#  attempt to set 'colnames' on an object with less than two dimensions


aa4.go = GetProteinGOInfo(aa4)
aa5.go = GetProteinGOInfo(aa5)

#if want a publication ready graph
PlotGOAll(GOObj = aa2.go, Top = 10, directorypath = getwd(), width = 8, height = 5)




#GO information
aa1.go$Gene.Ontology..GO. #1] "membrane [GO:0016020]; ammonia monooxygenase activity [GO:0018597]"
aa5.go$Gene.Ontology..GO.#"membrane [GO:0016020]; monooxygenase activity [GO:0004497]"



#####Now using accession numbers from the hw

aa6 = "P0A799"
aa7 = "P08839"

#aa6
aa6.go = GetProteinGOInfo(aa6)

PlotGoInfo(aa6.go)

#publication ready plot
PlotGOAll(GOObj = aa6.go, Top = 10, directorypath = getwd(), width = 8, height = 5)

#aa7
aa7.go = GetProteinGOInfo(aa7)

PlotGoInfo(aa7.go)




aa.6.path = GetPathology_Biotech(aa6.go)
#Internet connection problem occurs and the function will return the original error
#URL using bad/illegal format or missing URLInternet connection problem occurs
#NULL


#use output of GetPathology_Biotech

#Get.diseases(Dataframe retrieved from UniprotR Function "GetPathology_Biotech")





##accessing structural information about the protein with protti package
#using amoA gene
aa1.df =as.data.frame(fetch_uniprot(aa1))


fetch_pdb(aa1.df$xref_pdb)#Error in fetch_pdb(aa1.df$xref_pdb) : No PDB IDs were provided.



#using Phosphoglycerate kinase

aa6.df =as.data.frame(fetch_uniprot(aa6))

summary(aa6.df)

fetch_pdb(aa6.df$xref_pdb)#no ID...has a colon after 1ZMR aa6.df$xref_pdb = 1ZMR;


#using example from hw

aa6.structure =fetch_pdb("1ZMR")


#get protein structure...using uniprot identifiers
aa6.prediction =fetch_alphafold_prediction(aa6)


aa6_pdb = c("1ZMR")

#view 3d structure
r3dmol() %>%
  m_add_model(data = aa6_pdb, format = "pdb") %>%
  m_zoom_to()
