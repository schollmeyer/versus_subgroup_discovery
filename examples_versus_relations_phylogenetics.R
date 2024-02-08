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
gamma <- 0.0018
D3 <- (gamma*D1+(1-gamma)*D2)


i <- 2
i<-6
data(anole.data)
y <- anole.data[,i]
lambda=1
y <- lambda*y +(1-lambda)*runif(length(y),min=min(y),max=max(y))
y <- y  > quantile(y,0.5)# & y < quantile(y,0.75)
y <- oofos::compute_objective(data.frame(y=y),"y","TRUE")
table(y)

Z <- oofos:::get_auto_conceptual_scaling(anole.data)

y <- ddandrda::compute_tukeys_depth(Z,Z)


y <- y  > quantile(y,0.25) & y < quantile(y,0.75)

y <- oofos::compute_objective(data.frame(y=y),"y","TRUE")
context_1 <- get_context_from_distance(D1,complemented=FALSE)
context_2 <- get_context_from_distance(D2,complemented=FALSE)
context_3 <- get_context_from_distance(D3,complemented=FALSE)


####

vcdims <- local_object_VCdims(t(context_1),outputflag = 0,timelimit=1000)

C <- vcdims$vcdims+0.25*vcdims$vccounts
plot(ecdf(C))

j <- which(C <=quantile(C,0.9))
context_1 <-context_1[,j]
####

L_1 <- oofos:::compute_concept_lattice(context_1)
L_2 <- oofos:::compute_concept_lattice(context_2)
L_3 <- oofos:::compute_concept_lattice(context_3)



vc_dimension_1 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_1))$objval
vc_dimension_2 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_2))$objval
vc_dimension_3 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_3))$objval

versus_discovery_1 <- oofos::optimize_on_context_extents(context_1,objective=y)
result_1 <- gurobi::gurobi(versus_discovery_1)
versus_discovery_1$objval <- result_1$objval
test_1 <- oofos::compute_extent_optim_test(versus_discovery_1,n_rep=10000)


versus_discovery_2 <- oofos::optimize_on_context_extents(context_2,objective=y)
result_2 <- gurobi::gurobi(versus_discovery_2)
versus_discovery_2$objval <- result_2$objval
test_2 <- oofos::compute_extent_optim_test(versus_discovery_2,n_rep=10000)


versus_discovery_3 <- oofos::optimize_on_context_extents(context_3,objective=y)
result_3 <- gurobi::gurobi(versus_discovery_3)
versus_discovery_3$objval <- result_3$objval
test_3 <- oofos::compute_extent_optim_test(versus_discovery_3,n_rep=10000)

## lokale objekt / Attribute vc dimensions/counts fuer regularisierung