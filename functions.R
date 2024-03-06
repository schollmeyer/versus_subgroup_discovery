check_three_point_condition <- function(dist_mat,eps=10^-6){
  m <- ncol(dist_mat)
  counterexamples <- array(0,c(m,m))
  for( k in (1:m)){
    for(l in (1:m)){
      counterexamples[k,l] <- sum(dist_mat[k,] > eps+ (pmax(dist_mat[k,l],dist_mat[l,])))
    
      
      #counterexamples[k,l] <- sum(dist_mat[k,l] > eps+ (pmax(dist_mat[k,],dist_mat[,l])))
    }
    
    
    
  }
  #print(sum(counterexamples>0))
return(counterexamples)  
  
}












get_context_from_distance <- function(dist_mat,threshold,complemented=TRUE,sampling_proportion=1,remove_duplicates=TRUE,set_seed=TRUE,seed=1234567,eps=10^-10,eps2=10^-10){
  counterexamples <- check_three_point_condition(dist_mat,eps=eps2)
  n_rows <- nrow(dist_mat)
  n_rows_sample <- ceiling(sampling_proportion*n_rows)
  if(set_seed){set.seed(seed)}
  sampled_indexs <- sample(seq_len(n_rows),size=n_rows_sample)
  if(sampling_proportion==1){sampled_indexs <- seq_len(n_rows)}
  context <- array(0,c(n_rows,n_rows_sample*(n_rows_sample-1)))
  t <- 1
  for(k in seq_len(n_rows_sample)){
     for(l in seq_len(n_rows_sample)[-k]){
       if(counterexamples[k,l] <= threshold){
         print(counterexamples[k,l])
	     #print(t)
	     context[,t] <- (dist_mat[,sampled_indexs[k]] > dist_mat[,sampled_indexs[l]] +eps)
       
	      t <- t + 1
       }
	   
	  }
	 }
	 
if(complemented){context <- (cbind(context,1-context))}
if(remove_duplicates){context <- t(unique(t(context)))}
return(context)}

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
