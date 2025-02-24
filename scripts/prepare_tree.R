
require(ape)
t1 <- read.nexus("~/Library/CloudStorage/OneDrive-DurhamUniversity/Christophe_Jonathan_Shared_Folder/Titia phylogenetic phenotype (Chapter 4)/titia.mysnps-SNAPP-hydro_5_max_cov.trees 2.Anon")
to.keep = c("7050000010-NB0103","7050002710-ZANAa05","7050041410-RCJRa01","7050045480-ACSFc04","7050048620-CVa06","7050822250-CYa04","7050849600-TXRSa04","7050054670-BALB05")

t2 <-keep.tip(t2,to.keep)

#t2$tip.label
#[1] "7050000010-NB0103"  "7050002710-ZANAa05" "7050041410-RCJRa01" "7050045480-ACSFc04" "7050048620-CVa06"  
#[6] "7050054670-BALB05"  "7050822250-CYa04"   "7050849600-TXRSa04"


write.tree(t2,file="titia_tree_trimmed.tree")
