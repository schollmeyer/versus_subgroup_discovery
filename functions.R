get_context_from_distance <- function(dist_mat,complemented=TRUE,sampling_proportion=1,remove_duplicates=TRUE,set_seed=TRUE,seed=1234567,eps=10^-10){

  n_rows <- nrow(dist_mat)
  n_rows_sample <- ceiling(sampling_proportion*n_rows)
  if(set_seed){set.seed(seed)}
  sampled_indexs <- sample(seq_len(n_rows),size=n_rows_sample)
  if(sampling_proportion==1){sampled_indexs <- seq_len(n_rows)}
  context <- array(0,c(n_rows,n_rows_sample*(n_rows_sample-1)))
  t <- 1
  for(k in seq_len(n_rows_sample-1)){
     for(l in seq(k+1,n_rows_sample)){
	   print(t)
	   context[,t] <- (dist_mat[,sampled_indexs[k]] >= dist_mat[,sampled_indexs[l]] +eps)
	   t <- t + 1
	   
	  }
	 }
	 
if(complemented){context <- (cbind(context,1-context))}
if(remove_duplicates){context <- t(unique(t(context)))}
return(context)}
