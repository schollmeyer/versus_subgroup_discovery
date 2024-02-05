library(phytools)
data(anoletree)
tree_1 <- anoletree
pdf("tree_1.pdf",width=20,height=20)
plot(tree_1)
dev.off()

tree_2 <- force.ultrametric(tree_1)

pdf("tree_2.pdf",width=20,height=20)
plot(tree_2)
dev.off()



D1 <- cophenetic.phylo(tree_1)
D2 <- cophenetic.phylo(tree_2)

data(anole.data)
y <- anole.data[,1]
y <- y > median(y)
y <- oofos::compute_objective(data.frame(y=y),"y","TRUE")
table(y)
