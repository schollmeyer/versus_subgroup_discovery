##########################################
# Import data (code from hemberg-lab): https://github.com/hemberg-lab/scRNA.seq.datasets/blob/master/R/treutlein.R
#########################################

## set working directory accordingly !!!

# load additional functions for computing stylized betweenness:
#source("used_stylized_betweenness_functions.R")
# load additional functions needed for gene expression data (gene filter and scaling has to be applied):
source("../paperproject_starshaped-subgroup_discovery/additional_functions_for_gene_expression_data.R")

### DATA
d <- read.table("datasets/nature13173-s4.txt")
d <- t(d)
genes <- d[,1][5:nrow(d)]
# remove genes and bulk samples
d <- d[,2:(ncol(d) - 2)]
exprs_data <- as.data.frame(matrix(as.numeric(d[5:nrow(d),]), ncol = ncol(d)))
rownames(exprs_data) <- genes
colnames(exprs_data) <- d[1,]

### ANNOTATIONS
ann <- data.frame(cell_type1 = d[4,])
rownames(ann) <- d[1,]



############################################

# manual construction of x and y:
x <- t(exprs_data)
x <- gene_filter(x)
x <- log2(1+x)
dim(x)
x <- scaling(x)



# This context has VC dimension 80!!!

y <- ann$cell_type1
table(y)/length(y)

# y
#  AT1      AT2       BP ciliated    Clara
#  0.5125   0.1500   0.1625   0.0375   0.1375
set.seed(1234567)
idxs <- sample((1:80),size=40)
x_all <- x
y_all <- y
x <- x_all[idxs,]
y <- y_all[idxs]

# permute gene expression values for randomly selected genes
i=seq(1,6751,length.out=6100)
i=sample((1:6751),size=6100)
for(k in i){x[,k] <-sample(x[,k])}
context <- oofos:::get_auto_conceptual_scaling(x)
objective <- oofos:::compute_objective(data.frame(y=y),"y","AT1")



dist_mat_treutlein <- get_distance_from_context(context)
#saveRDS(dist_mat_treutlein,"dist_mat_treutlein.RDS")
#dist_mat_treutlein <- readRDS("dist_mat_treutlein.RDS")




D_plus <- dist_mat_treutlein[which(dist_mat_treutlein>0)]

eps <- c(0,.8*quantile(D_plus,seq(0,0.99,length.out=19)))

fitted_pseudoultrametrics <- list()
for(k in (1:20)){
ans <- fit_ultrametric(dist_mat_treutlein,eps=eps[k],upper_bound=4*max(dist_mat_treutlein)+2*eps[k],start_solution=TRUE)
gc()
fitted_pseudoultrametrics[[k]] <- try(gurobi::gurobi(ans,list(timelimit=60*60)))
}

D_ultra <- (fitted_pseudoultrametrics[[1]])$x[(1:1600)];dim(D_ultra) <- c(40,40)

D_ultra <- dist_mat_treutlein
check_ultrametric_violations(D_ultra)

tree_add <-  ape::nj(dist_mat_treutlein)
D_add <- ape::cophenetic.phylo(tree_add)

context_ultra <- get_context_from_distance(D_add,lambda=1,threshold=100,eps2=0,complemented=FALSE)
dim(context_ultra)
# 32
vc_dimension_ultra <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context_ultra))$objval

objective2 <- sample(objective)
discovery <- oofos::optimize_on_context_extents(context_ultra,objective=objective2)
result <- gurobi::gurobi(discovery)
discovery$objval <- result$objval

##
Z <- oofos:::get_non_stylized_betweenness(context_ultra)
discovery2 <- oofos::discover_starshaped_subgroups(stylized_betweenness=Z,objective=objective2,complexity_control=Inf)
discovery2$objval-result$objval


oofos::compute_extent_optim_test(discovery)

set.seed(1234567)
indexs <- sample((1:80),size=9)



D <- get_distance_from_context(context[,])

ans <- fit_ultrametric(D,eps=0,start_solution=TRUE)
B <- gurobi::gurobi(ans)
BB=BsaveRDS(gbsb,"results_treutlein_gbsb/absb.RDS")


# gbsb:
starshaped_discovery_gbsb <- oofos::discover_starshaped_subgroups(stylized_betweenness=gbsb,objective=objective,local_vc_dimension=8)
starshaped_discovery_gbsb$objval
test_gbsb <- oofos:::compute_starshaped_distr_test(starshaped_discovery_gbsb, n_rep=10000)
saveRDS(starshaped_discovery_gbsb,"results_treutlein/starshaped_discovery_gbsb")
saveRDS(test_gbsb,"results_treutlein/test_gbsb")


##### extensive analysis for gbsb
set.seed(1234567)
vc_dimensions <- seq(1,80,length.out=500)
p_values <- rep(0,length(vc_dimensions))
p_values_param <- rep(0,length(vc_dimensions))
objvalues <- array(0,c(500,100))
objval <- rep(0,100)
for(k in (1:500)){
  discovery <- oofos::discover_starshaped_subgroups(gbsb,objective=objective,complexity_control = vc_dimensions[k])
  test <- oofos::compute_starshaped_distr_test(discovery,n_rep=100)
  objvalues[k,] <- test$objvalues
  objval[k] <- discovery$objval
  p_values[k] <- test$p_value
  p_values_param[k] <-  (discovery$objval-mean(test$objvalues))/sd(test$objvalues)
  plot(vc_dimensions,(p_values_param))
}


