
#Install ggtree

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


BiocManager::install("treeio")



library(ggtree)
library(phangorn)
library(treeio)
library(tidytree)



#####Uploading tree##

#a formal class tree data
tree.rax <-read.raxml(file = "RAxML_bipartitionsBranchLabels.picocyano.msa.v3_1.tre")


#character string with asv names/tip labels
rax.tip <- as.phylo(tree.rax)$tip.label

#create "new" tree with only fSyn grouped together
tree.rax_c <- groupOTU(tree.rax, rax.tip[grep("fSyn", rax.tip)],
                       group_name = "sample.type")#

#Circular tree with paths  
ggtree(tree.rax_c,layout="circular",size = 1.5, aes(color=sample.type)) +  theme(legend.position='bottom') +
  scale_color_manual(values=c("lightblue", "lightgreen"), 
                     labels=c("Marine", "Freshwater"))+geom_tiplab(aes(angle=angle),size = 5.5, color="black")



#Use figtree for type strains....cant see association in current trees with ggtree





#####TEST from  ggtree book#############33



#nwk <- system.file("extdata", "sample.nwk", package="treeio")


#tree <- read.tree(nwk)

#No lables only structure of branchs
#ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()


#ggtree(tree, color="firebrick", size=2, linetype="dotted")



#ggtree(tree)
#ggtree(tree, ladderize=FALSE)#
#ggtree(tree, branch.length="none")#branch lengths all equal

#ggtree(tree) + 
#  geom_point(aes(shape=isTip, color=isTip), size=3)


#p <- ggtree(tree) + 
#  geom_nodepoint(color="#b5e521", alpha=1/4, size=10) #nodepoint sets colur of nodes
#p + geom_tippoint(color="#FDAC4F", shape=8, size=3)# geom_tippoint sets color of tippoints 


#p + geom_tiplab(size=3, color="purple")

#ggtree(tree) + geom_treescale()
#####uncolored#####
#ggtree(tree, layout="circular") + geom_tiplab(aes(angle=angle), color='blue')

#highlighting specific clades....only highlighting before the node not the whole clade like in the book
#ggtree(tree) + 
#  geom_hilight(node=21, fill="steelblue", alpha=.6, type="rect") +
#  geom_hilight(node=17, fill="darkgreen", alpha=.6) 

#ggtree(tree, layout="circular") + 
#  geom_hilight(node=21, fill="steelblue", alpha=.6) +
#  geom_hilight(node=23, fill="darkgreen", alpha=.6,type = "encircle")



#ggtree(tree) +
#  geom_balance(node=16, fill='steelblue', color='white', alpha=0.6, extend=1) +
#  geom_balance(node=19, fill='darkgreen', color='white', alpha=0.6, extend=1) 

#Warning messages:
#  1: Computation failed in `stat_balance()`.
#Caused by error in `map()`:
#  ℹ In index: 1.
#Caused by error in `offspring.tbl_tree_item()`:
#  ! could not find function "offspring.tbl_tree_item" 
#2: Computation failed in `stat_balance()`.
#Caused by error in `map()`:
#  ℹ In index: 1.
#Caused by error in `offspring.tbl_tree_item()`:
#  ! could not find function "offspring.tbl_tree_item"





#beast_file <- system.file("examples/MCC_FluA_H3.tree", 
#                          package="ggtree")
#beast_tree <- read.beast(beast_file)
#ggtree(beast_tree, mrsd="2013-01-01") + theme_tree2()


## separate the tree by host species
#tip <- as.phylo(beast_tree)$tip.label
#beast_tree <- groupOTU(beast_tree, tip[grep("Swine", tip)], 
#group_name = "host")


#p <- ggtree(beast_tree, aes(color=host)) + 
#  theme_classic() + theme(legend.position='none') +
#  scale_color_manual(values=c("blue", "red"), 
#                     labels=c("human", "swine")) 


#Clade labels

#set.seed(2015-12-21)
#tree <- rtree(30)
#p <- ggtree(tree) + xlim(NA, 8)#

#p + geom_cladelab(node=45, label="test label") +
#  geom_cladelab(node=34, label="another clade")



#Error in offspring.tbl_tree_item(.data = .data, .node = .node, tiponly = tiponly,  : 
#                                   could not find function "offspring.tbl_tree_item"

