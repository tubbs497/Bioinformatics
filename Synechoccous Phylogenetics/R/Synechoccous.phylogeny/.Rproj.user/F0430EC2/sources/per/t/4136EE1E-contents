

library(phyloseq)
library(tidyverse)
#load data sets generated from classification and denoising/error correction (dada2)
#Counts of ASVs
counts =read.csv("wc1.syn.epa.shelf.feature.table.csv")

#ASV Taxa information

taxa = read.csv("wc1.syn.epa.shelf.shelf.taxa.table.csv")


##Merging data sets by asv.id

merge = merge(taxa, counts, by = "ASV.ID")


write.csv(merge, file = '/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/Synechoccous Phylogenetics/R/Synechoccous.phylogeny/wc1.syn.asv_count_taxa.csv')



####Create Phyloseq object######################


# ASV Counts

count.data = read.csv ("wc1.syn.aligned.asv_count.csv",
                         row.names = 1)

attach(count.data)
#removing q2 ids
count.data = count.data[,-c(1:2)]

# Make matrix
count.data = as.matrix(count.data)

# Convert for phyloseq creation
ASV = otu_table(count.data,
                taxa_are_rows = TRUE)

#*******************************************Taxonomy****************************
#Aligned Taxa
taxa.data = read.csv ("wc1.syn.aligned.taxa.csv",
                      row.names = 1)

attach(taxa.data)
taxa.data = taxa.data[,-c(1:3)]

# Make a matrix
taxa.data = as.matrix(taxa.data)

# Convert for phyloseq creation
TAX = tax_table(taxa.data)

#***************************************Metadata********************************

metadata = read.csv ("wc1.syn.epa.metadata.csv",
                     row.names = 1)

attach(metadata)


# Convert for phyloseq creation
META = sample_data(metadata)

#***********************************MAKE PHYLOSEQ OBJECT************************
wc1_syn.phy = phyloseq(ASV,
                           TAX,
                           META)




############Validating Phyloseq obj.################
library(microViz)
wc1_syn.phy = tax_fix(wc1_syn.phy)

wc1_syn.phy = phyloseq_validate(wc1_syn.phy, remove_undetected = TRUE)


#Convert counts to rel.abunance
ps.ra = transform_sample_counts(wc1_syn.phy, function(OTU) OTU/sum(OTU) )


syn = subset_taxa(ps.ra, Genus %in% "g__Synechococcus_CC9902")

#none in freshwater samples... did not classify synechococcuus (reported in Hancock et al. 2024) in freshwater samples at genus level.
plot_bar(syn)+geom_bar(aes(fill=ASV), stat="identity", position="stack")



###Should be at genus level, but classification may not be that resolute

syn.f= subset_taxa(ps.ra, Family %in% "f__Synechococcaceae")
#Only in some of frank and lock monitoring

plot_bar(syn.f)
###use to determine outgroup sequences?
syn.o =  subset_taxa(ps.ra, Order %in% "o__Synechococcales")

plot_bar(syn.o)+geom_bar(aes( fill=Genus), stat="identity", position="stack")

##Cyanobium_PCC-6307 dominant in freshwater...in marine samples aswell..highest at N5


plot_bar(syn.o)+geom_bar(aes( fill=ASV), stat="identity", position="stack")



#Synechoccous is in cyanobiaceae family 
syn.cyanobium =subset_taxa(ps.ra, Family %in% "f__Cyanobiaceae")  


plot_bar(syn.cyanobium)+geom_bar(aes( fill=Genus), stat="identity", position="stack")

####use rep seqs from f__Cyanobiaceae
plot_bar(syn.cyanobium)+geom_bar(aes( fill=ASV), stat="identity", position="stack")


#export f__Cyanobiaceae table to excell to determine most prevalent ASVs
#counts
syn.asv = as.matrix(syn.cyanobium@otu_table)
syn.asv = as.data.frame(syn.asv)
#taxa
syn.tax = as.matrix(syn.cyanobium@tax_table)
syn.tax = as.data.frame(syn.tax)

syn.asv.tax = merge(syn.tax,syn.asv, "row.names")

row.names(syn.asv.tax) = syn.asv.tax[,1]

syn.asv.tax = syn.asv.tax %>% mutate(asv.id = Row.names)


#Adding q2 asv identifiers 

count.data1 = read.csv ("wc1.syn.aligned.asv_count.csv",
                       row.names = 1)

syn.asv.tax.merge = merge(count.data1, syn.asv.tax, "ASV")

write.csv(syn.asv.tax.merge, file = '/Users/trevo/OneDrive/Documents/GitHub/Bioinformatics/Synechoccous Phylogenetics/R/Synechoccous.phylogeny/wc1.cyanobiaceae_table.csv')





###looking at asvs in cyanobium (in freshwater and saltwater)
cyanobium = subset_taxa(ps.ra, Genus %in% "g__Cyanobium_PCC-6307")

plot_bar(cyanobium, "Sample.type")+geom_bar(aes(fill=ASV), stat="identity", position="stack")

####Trying nitrosomonas?


nitroso = subset_taxa(ps.ra, Family %in% "f__Nitrosomonadaceae")



plot_bar(nitroso)+geom_bar(aes( fill=Genus), stat="identity", position="stack")
