library(phytools)
data(anoletree)
tree_1 <- anoletree
pdf("tree_1.pdf",width=20,height=20)
plot(tree_1)
dev.off()
