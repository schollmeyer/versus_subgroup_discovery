compute_random_dist <- function(n,dimension){
   x <- rnorm(n*dimension);dim(x) <- c(n,dimension)
   return(as.matrix(dist(x)))
  
}
ultra_dist <- function(distances){.0*distances + as.matrix(clue::ls_fit_ultrametric(distances))}
id <- function(distances){distances}
add_dist <- function(distances){ape::cophenetic.phylo(ape::nj(distances))}



get_distance_from_context <- function(context,bandwidth,method="manhattan",normalized=TRUE){
  result <- as.matrix(dist(context,method="manhattan"))
  if(normalized){result <- result/ncol(context)}
return(result)}

get_distance_from_extens <- function(extents){
  n_col <- ncol(extents)
  D <- array(0,c(n_col,n_col))
  for(k in seq_len(n_col)){
    for(l in seq_len(n_col)){
      D[k,l] <- sum(extents[,k ] & extents[,l]) /sum(extents[,k] |extents[,l])
    }
  }
  return(D)}  

get_halfspace_distance <- function(context){# same as get_distance_from_context?
  n_row <- nrow(context)
  D <- array(0,c(n_row,n_row))
  for(k in (1:n_row)){
    for(l in (1:n_row)){
      D[k,l] <- sum(context[k,] != context[l,]) 
      
    }}
return(D)}
  


 # n_objects <- nrow(context)
 # n_attributes <- ncol(context)
 # result <- array(0,c(n_objects,n_objects))
 # context <- t(context)
  #attribute_distances <- as.matrix(dist(t(context)))
 # for(k in (1:n_objects)){
 #   print(k)
 #   for(l in (1:n_objects)){
 #     #for(m in (1:n_attributes)){
 #     result[k,l] <- sum(abs(context[,k]-context[,l]))#dnorm(attribute_distances[m,],sd=bandwidth)* abs(contex[k,]-context[l,])))
      #}
 #   }

 # }
#  if(normalized){result <- result/ncol(context)}
#  return(result)

#}



check_three_point_condition <- function(dist_mat,eps=10^-6,lambda=1){
  m <- ncol(dist_mat)
  counterexamples <- array(0,c(m,m))
  for( k in (1:m)){
    print(k)
    for(l in (1:m)){
      # old version
      #counterexamples[k,l] <- lambda*sum(dist_mat[k,] > eps+ (pmax(dist_mat[k,l],dist_mat[l,])))
     counterexamples[k,l] <-  sum(dist_mat[,k] > eps+ (pmax(dist_mat[k,l],dist_mat[,l])))
     
      #counterexamples[k,l] <- counterexamples[k,l]+ (1-lambda)*sum(pmax(dist_mat[k,] -( eps+ (pmax(dist_mat[k,l],dist_mat[l,]))),0))



      #counterexamples[k,l] <- sum(dist_mat[k,l] > eps+ (pmax(dist_mat[k,],dist_mat[,l])))
    }



  }
  #print(sum(counterexamples>0))
return(list(counterexamples=counterexamples,result=all(counterexamples==0)))

}












get_context_from_distance <- function(dist_mat,threshold,complemented=FALSE,indexs=NULL,sampling_proportion=1,remove_duplicates=TRUE,set_seed=TRUE,seed=1234567,eps=10^-10,eps2=10^-10,lambda=1,counterexamples = check_three_point_condition(dist_mat,eps=eps2,lambda=lambda)$counterexamples){
  n_rows <- nrow(dist_mat)
  n_rows_sample <- ceiling(sampling_proportion*n_rows)
  if(set_seed){set.seed(seed)}
  sampled_indexs <- sample(seq_len(n_rows),size=n_rows_sample)
  if(sampling_proportion==1){sampled_indexs <- seq_len(n_rows)}
  if(!is.null(indexs)){sampled_indexs <-indexs;n_rows_sample <- length(indexs)}
  context <- array(FALSE,c(n_rows,n_rows_sample*(n_rows_sample-1)))
  ce <- rep(0,n_rows_sample*(n_rows_sample-1))
  #counterexamples <- check_three_point_condition(dist_mat[sampled_indexs,sampled_indexs],eps=eps2,lambda=lambda)$counterexamples
  t <- 1
  col_names <- NULL
  for(k in seq_len(n_rows_sample)){
    print(k)
     for(l in seq_len(n_rows_sample)[-k]){
       if(counterexamples[k,l] <= threshold){
         #print(counterexamples[k,l])
	     #print(c(k,l))
	     context[,t] <- (dist_mat[,sampled_indexs[k]] > dist_mat[,sampled_indexs[l]] +eps)
	     ce[t] <- counterexamples[k,l]
       #col_names <- c(col_names,paste(k,">",l,collaps=""))
       #print(col_names)
	      t <- t + 1
       }

	  }
  }
  context <- context[,(1:(t-1))]
colnames(context) <- col_names
if(complemented){context <- (cbind(context,1-context))}
if(remove_duplicates){
  non_duplicates <- which(!duplicated(t(context)))
  context <- context[,non_duplicates]
  ce <- ce[non_duplicates]
  
  }
return(list(context=context,counterexamples=ce))}

regularize_tree <- function(tree,lambda){

  h <- diag(vcv(tree))
  d <- max(h) - h
  ii <- sapply(1:Ntip(tree), function(x, y) which(y ==
                                                    x), y = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + d - lambda
  return(tree)}


local_object_VCdims=function(X,indexs=(1:dim(X)[1]),outputflag,timelimit,pool=FALSE,transpose=TRUE,additional.constraint=TRUE,threads=1){
  # Berechnet lokale Gegenstands-VC-dimension:
  # X : Kontext
  # indexs: Indizes derjenigen Punkte, für die die lokale Gegenstands-VC-Dimension berechnet werden soll
  # outputflag: Argument, das an urobi übergeben wird (zur Steuerung der Ausgabe während der Berechnung)
  # timelimit: Zeitlimit für die Berechnun einer einzelnen lokalen VC-Dimension
  # pool: Wenn True, dann werden alle Kontranominalskalen maximaler Kardinalität berechnet, ansonsten nur eine
  # Transpose: ob Kontext vorher transponiert werden soll: Bei Kontext mit mehr Zeilen als Spalten scheint mit transpose =TRUE die Berechnung schneller zu laufen, für mehr Spalten als Zeilen scheint transpose=FALSE oft schneller zu sein
  #additional.constraint: ob zusätzlicher Constraint (Anzahl Gegenstände der Kontranominalskala==Anzahl Merkmele der Kontranominalskala) mit implementiert werden soll (dadurch wird Berechnung oft leicht schneller)

  m=dim(X)[1]
  n=dim(X)[2]
  ans=list()
  vcdims=rep(0,m)
  vccounts=rep(0,m)
  for(k in indexs){
    if(transpose){
      temp=compute_extent_vc_dimension(t(X),additional_constraint=additional.constraint)
      temp$lb[k+n]=1
    }
    else{
      temp=extent.VC((X))
      temp$lb[k]=1
    }


    if(pool){a=gurobi(temp,list(outputflag=outputflag,timelimit=timelimit,PoolSolutions=100000000,PoolSearchMode=2,Poolgap=0.00001))}
    else{a=gurobi(temp,list(outputflag=outputflag,timelimit=timelimit,threads=threads))}
    a <<- a
    ans[[k]]=a
    vcdims[k]=a$objval
    vccounts[k]=length(a$pool)
    print(a$objval)
  }
  return(list(vcdims=vcdims,vccounts=vccounts,rest=ans))}



 F <- function(dist_mat, delta){

   m <- ncol(dist_mat)
   indexs <- seq_len(m^2);dim(indexs) <- c(m,m)

   A <- array(0,c(m^2+m^2,m^2+m^2))
   sense <- rep("",m^2+m^2)
   rhs <- rep(0,m^2+m^2)
   t <- 1

   t <- 1
   ## ultrametric property
   genconmax <- list()
   for(k in (1:m)){
     for(l in (1:m)){
         genconmax[[t]]<- list( resvar = indexs[k,l]+m^2, vars = pmax(indexs[k,],indexs[,l]))
          A[t,indexs[k,l]] <- 1
          A[t,indexs[k,l]+m^2] <- -1
          sense[t] <- "<"
          rhs[t] <- 0

         t <- t+1

     }}


   ## symmetry
   for(k in (1:m)){
     for(l in (1:m)){
       A[t,indexs[k,l]] <- 1
       A[t,indexs[l,k]] <- -1
       sense[t] <- "="
       rhs[t] <- 0
       t <- t+1

     }
   }

   model <- list(A=A,sense=sense,rhs=rhs,obj=NULL,genconmin=genconmax,lb=c(as.vector(dist_mat)-delta,rep(0,m^2)),ub=c(as.vector(dist_mat)+delta,rep(Inf,m^2)))
   result <- gurobi::gurobi(model)
   print(result$status)
   result <- result$x[(1:(m^2))];dim(result) <- c(m,m)
   return(result)

 }


 check_ultrametric_violations <- function(D){
   ans <- 0
   n_objects <- nrow(D)
   dim(D)=c(n_objects,n_objects)
   for(k in (1:n_objects)){
     for(l in (1:n_objects)){

       ans <- ans+ sum(pmax(0,D[k,l]-pmax(D[k,],D[,l]))^2)
     }}

   return(ans)}

 fit_ultrametric <- function(dist_mat,eps=rep(0,nrow(dist_mat)^2+nrow(dist_mat)^3), solve=FALSE,param=NULL,start_solution=FALSE,upper_bound=4*max(dist_mat),bounded=FALSE,eps_lower,eps_upper){
   x <- rep(0,nrow(dist_mat)^2+nrow(dist_mat)^3)
   print(dim(dist_mat))
   n_objects <- nrow(dist_mat)
   ub=rep(upper_bound,n_objects^2+n_objects^3)
   lb=rep(0,n_objects^2+n_objects^3)
   if(bounded){
     ub[seq_len(n_objects^2)] <- as.vector(dist_mat)+eps_upper
     lb[seq_len(n_objects^2)] <- pmax(0,as.vector(dist_mat)-eps_lower)
   }
   start <- NULL
   if(start_solution){

     D_heuristic <- as.matrix(clue::ls_fit_ultrametric(dist_mat))
     start=as.vector(D_heuristic)
     t <- 1

     start_2 <- rep(0,n_objects^3)
     for(k in (1:n_objects)){
       for(l in (1:n_objects)){
         for(m in (1:n_objects)){
           start_2[t] <- max(D_heuristic[k,m],D_heuristic[m,l])
           t <- t+1

         }}}

     start <- c(start,start_2)}

   genconmax=list()

   mat <- (1:n_objects^2)
   dim(mat) <- c(n_objects,n_objects)
   mat2 <- (1:n_objects^3)
   dim(mat2) <- c(n_objects,n_objects,n_objects)
   sense <- rep("",n_objects^2+n_objects^3)
   RHS <- rep(0,n_objects^2+n_objects^3)
   i <- j <- v <- rep(0,n_objects^2+n_objects^3)
   t <- 1
   obj <- rnorm(n_objects^2+n_objects^3)
   x[(1:n_objects^2)] <- as.vector(dist_mat)
   t <- 1
   for(k in (1:n_objects)){
     for(l in (1:n_objects)){
       for(m in (1:n_objects)){
         i[t] <- t
         j[t] <- mat[k,l]
         v[t] <- 1

         i[t+n_objects^3] <- t
         j[t+n_objects^3] <- mat2[k,m,l]+n_objects^2
         v[t+n_objects^3] <- -1
         #print(t-n_objects^2)
         genconmax[[t]] <- list(resvar = mat2[k,m,l]+n_objects^2,vars=c(mat[k,m],mat[m,l]))

         #TODO

         #print((dist_mat[k,l] - max(dist_mat[k,m],dist_mat[m,l])))
         RHS[t] <- (dist_mat[k,l] - (max(dist_mat[k,m],dist_mat[m,l])))

         x[t+n_objects^2] <- (max(dist_mat[k,m],dist_mat[m,l]))
         t <- t + 1

       }
       #tt <- tt+1
     }}

   t <- 1
   for(k in (1:n_objects)){
     for(l in (1:n_objects)){
       if(k==l){ub[t] <- 0}
       t <- t+1
     }
   }


   Q=slam::simple_triplet_matrix(i=(1:(n_objects^2+n_objects^3)),j=(1:(n_objects^2+n_objects^3)),v=c(rep(1,n_objects^2),rep(0,n_objects^3)))
   result <-(list(x=x,ub=ub,start=start,RHS2=RHS,rhs=eps,Q=Q,A=slam::simple_triplet_matrix(i, j, v, nrow=n_objects^2+n_objects^3),obj= (-2)*c(as.vector(dist_mat),rep(0,n_objects^3)),lb=lb,alpha=sum(dist_mat*dist_mat),genconmax=genconmax,sense=rep("<=",n_objects^2+n_objects^3 ) ))
   if(solve){
     result <- gurobi::gurobi(result,param=param)$x[seq_len(n_objects^2)]
     dim(result) <- c(n_objects,n_objects)
   }
  return(result)
  }





#compute_phi <- function(attribute_set, context){temp <- matrix(context[,attribute_set],ncol=sum(attribute_set));Rfast::rowAll(temp,parallel=FALSE)}
#compute_psi <- function(object_set, context){temp <- matrix(context[object_set,],nrow=sum(object_set));Rfast::colAll(temp,parallel=FALSE)}

compute_object_closure <- function(object_set,context){compute_phi(compute_psi(object_set,context),context)}

is_halfspace <- function(extent, index, context){
  extent_2 <- compute_object_closure(extent,context)
  if(extent[index]==1){return(FALSE)}
  if(any(extent_2!=extent)){return(FALSE)}
  for(k in which(extent==0)){
	extent_2 <- extent
	extent_2[k] <- 1
	extent_2 <- compute_object_closure(extent_2,context)
	if(any(extent_2[which(extent==0)]==1) & extent_2[index]==0){return(FALSE)}
  }

return(TRUE)

}



get_halfspace_context <- function(context){
 n_objects <- nrow(context)
 object_sets <- gtools::permutations(2,n_objects,repeats.allowed=TRUE)-1
 #object_sets <- oofos:::compute_concept_lattice(t(context))$intents
 dims=dim(object_sets)
 object_sets <- as.logical(object_sets)
 dim(object_sets) <- dims
 result <- NULL
 for(k in (1:nrow(object_sets))){
 print(k)
    for(index in which(object_sets[k,]==0)){
      temp <- is_halfspace(object_sets[k,], index,context)

	  if(temp==TRUE){ print(k);print(temp);result <- cbind(result,object_sets[k,])}
	}
 }
 return(result)}


compute_versus_halfspace <- function(context,index_1,index_2,extent1 = rep(0,nrow(context)),extent2 = rep(0,nrow(context))){
  n_row <- nrow(context)
  halfspace <- rep(0,n_row)

  for(k in (1:n_row)){
    extent_1 <- extent1;extent_1[c(k,index_2)] <- 1
    extent_1 <- compute_object_closure(extent_1,context)

    extent_2 <- extent2;extent_2[c(k,index_1)] <- 1
    extent_2 <- compute_object_closure(extent_2,context)
    if(extent_1[index_1]==1 & extent_2[index_2]==0){halfspace[k] <- 1}


  }
return(halfspace)
}

