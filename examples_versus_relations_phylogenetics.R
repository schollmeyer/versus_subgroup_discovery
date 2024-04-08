# relaxierter ultrametrik-gurobi algorihmus noch implementieren

# ls_fit_ultrametric {clue}
library(clue)
library(phytools) # for phylogenetic datasets
library(maotai) # for checking if matrix is a distance matrix : checkmetric()
#covariates:

data(anoletree)
tree_1 <- anoletree
tree_2 <- force.ultrametric(tree_1)
tree_1 <- ape::nj(D1)

D1 <- cophenetic.phylo(tree_1)
set.seed(1234567)
D1 <- D1+ rnorm(length(D1),sd=0.01)
diag(D1) <- 0
D1 <- D1 +t(D1)
checkmetric(D1)
D2 <- cophenetic.phylo(tree_2)

context_1 <- get_context_from_distance(D1,complemented=FALSE,eps=10^-10,eps2=10^-10,threshold=0.5)
estimate_size_mingen_k_new((context_1),k=7)/10^6

context_2 <- get_context_from_distance(D2,complemented=FALSE,eps=10^-10,eps2=10^-10,threshold=0)
# context_3 <- get_context_from_distance(D3,complemented=FALSE)


fc1 <- fcaR::FormalContext$new(t(context_1))
fc1$find_concepts()
lattice_1 <- t(as.matrix(fc1$concepts$intents()))
dim(lattice_1)
max(lattice_1%*%y)
n_rep <- 10000
h01 <- rep(0,n_rep)
for(k in (1:n_rep)){h01[k] <- max(lattice_1%*%sample(y))}
plot(ecdf(h01))
mean(h01 <=max(lattice_1%*%y))

fc2 <- fcaR::FormalContext$new(t(context_2))
fc2$find_concepts()
lattice_2 <- t(as.matrix(fc2$concepts$intents()))
dim(lattice_2)
max(lattice_2%*%y)
n_rep <- 10000
h02 <- rep(0,n_rep)
for(k in (1:n_rep)){h02[k] <- max(lattice_2%*%sample(y))}
plot(ecdf(h02))
mean(h02 <=max(lattice_2%*%y))
lattice_1 <- oofos:::compute_concept_lattice(context_1);dim(lattice_1$extents)
lattice_2 <- oofos:::compute_concept_lattice(context_2);dim(lattice_2$extents)
#lattice_3 <- oofos:::compute_concept_lattice(context_3);dim(lattice_3$extents)



vc_dimension_1 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_1))$objval
vc_dimension_2 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_2))$objval
#vc_dimension_3 <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_3))$objval

# target variable:

data(anole.data)
y <- anole.data[,2]
#lambda <- .05
#y <- lambda*y +(1-lambda)*mean(y)#runif(length(y),min=min(y),max=max(y))
y <- y  > quantile(y,0.5) #& y < quantile(y,0.53)
y <- oofos::compute_objective(data.frame(y=y),"y","TRUE")
table(y)



#pdf("tree_1.pdf",width=20,height=20)
#plot(tree_1)
#dev.off()



#pdf("tree_2.pdf",width=20,height=20)
#plot(tree_2)
#dev.off()




gamma <- 0.0018
D3 <- (gamma*D1+(1-gamma)*D2)


i <- 2




table(y)

Z <- oofos:::get_auto_conceptual_scaling(anole.data)

y <- ddandrda::compute_tukeys_depth(Z,Z)


y <- y  > quantile(y,0.25) & y < quantile(y,0.75)

y <- oofos::compute_objective(data.frame(y=y),"y","TRUE")


######
tree_2 <- regularize_tree(tree_1,6)
D2 <- cophenetic.phylo(tree_2)
context_2 <- get_context_from_distance(D2,complemented=FALSE)
L_2 <- oofos:::compute_concept_lattice(context_2)
dim(L_2$extents)

####

vcdims <- local_object_VCdims(t(context_1),outputflag = 0,timelimit=1000)

C <- vcdims$vcdims+0.25*vcdims$vccounts
plot(ecdf(C))

j <- which(C <=quantile(C,0.9))
context_1 <-context_1[,j]
####







versus_discovery_1 <- oofos::optimize_on_context_extents(context_1,objective=y)
result_1 <- gurobi::gurobi(versus_discovery_1)
versus_discovery_1$objval <- result_1$objval
test_1 <- oofos::compute_extent_optim_test(versus_discovery_1,n_rep=4)


versus_discovery_2 <- oofos::optimize_on_context_extents(context_2,objective=y)
result_2 <- gurobi::gurobi(versus_discovery_2)
versus_discovery_2$objval <- result_2$objval
test_2 <- oofos::compute_extent_optim_test(versus_discovery_2,n_rep=10000)


versus_discovery_3 <- oofos::optimize_on_context_extents(context_3,objective=y)
result_3 <- gurobi::gurobi(versus_discovery_3)
versus_discovery_3$objval <- result_3$objval
test_3 <- oofos::compute_extent_optim_test(versus_discovery_3,n_rep=10000)

## lokale objekt / Attribute vc dimensions/counts fuer regularisierung