

XX=X
L1=oofos:::compute_concept_lattice(X)

E=NULL


for(k in (1:10)){
  print(k)
  for( l in (1:10)[-k]){

for(m in seq_len(nrow(L1$extents))){

 # for(m2 in seq_len(nrow(L1$extents))){




  if(TRUE | (L1$extents[m,k]==0 & L1$extents[m2,k]==0 & L1$extents[m,l]==0 & L1$extents[m2,l]==0)){

  temp <- compute_versus_halfspace(X,k,l,extent1=L1$extents[m,],extent2=0*L1$extents[m2,])
  E=cbind(E,compute_object_closure(temp,X))
  }

  }}}

E=NULL;for(k in (1:20)){for( l in (1:20)[-k]){E=cbind(E,compute_object_closure(compute_versus_halfspace(X,k,l),X))}}



L2=oofos:::compute_concept_lattice(E)
L2 <- L1
L2$extents <- L2$extents[(1:100),]
is_sublattice(L2,L1)



is_sublattice <- function(lattice_1,lattice_2){
  n1 <- nrow(lattice_1$extents)
  n2 <- nrow(lattice_2$extents)

  for(k in (1:n1)){
    result <- FALSE
    for(l in (1:n2)){
      if(all(lattice_1$extents[k,] == lattice_2$extents[l,])){result <- TRUE}
    }
  if(!result){print(k);return(FALSE)}
  }
  return(TRUE)
}
####


# geometry


set.seed(1234567);x=rnorm(100*2);dim(x)=c(100,2);plot(x,pch=16);CI=convex.incidence(x)


x=as.matrix(expand.grid((1:20),(1:20)));plot(x)
CI=convex.incidence(x)
CT <- CI$context
set.seed(1234567)
E <- rep(0,100)#sample(c(0,1),prob=c(0.01,0.99),size=500,replace=TRUE)
E[(3:4)]=1

CT <- CI$context[(1:100),]
E <- compute_object_closure(E,CT)
for(k in (1:9900)){
I <- rep(0,ncol(CT))
k=587
I[k] <- 1
I[6000]<-1
E <- compute_phi(I,CT)
#if(abs(sum(E)-200)<=5){print(k)}}
idx <- sample(which(E==0),size=2)
G <- compute_versus_halfspace(CI$context,idx[1],idx[2],E,E)
sum(G)
plot(x,pch=16)


points(x[which(G==1),],pch=16,col="red")
points(x[idx[(1:2)],],col="green",pch=16)
points(x[which(E==1),],col="blue",pch=16)
points(x[which(G==1),],col="red")
}
