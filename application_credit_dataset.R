# set working directory
setwd("C:/Git/paperproject_starshaped-subgroup_discovery")
# needed libraries
library(foreign)
library(rsubgroup)
# load the data set credit ( https://archive.ics.uci.edu/dataset/144/statlog+german+credit+data)
# Thisis already the curated dataset, compare TODO
# The variable sex (column 9) will be later deleted because it is actually not
# the variable sex but instead a combination of sex and marital status
data(credit.data)
dat <- credit.data
dim(dat)
# [1] 1000   21

#generate formal context of covariates (sex (column 9) and credit status (column 21) are excluded)


whole_context <- oofos:::get_auto_conceptual_scaling(dat[,-c(9,21)])
# For the classical subgroup discoery we make a sample split approach because
# optimizing under H0 is computationally very expensive
set.seed(1234567)
indexs <- sample((1:1000),size =50)
context <- oofos:::get_auto_conceptual_scaling(dat[indexs,-c(9,21)])
training_context <- whole_context[indexs,]
test_context <- whole_context[-indexs,]


vc_dim <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context))

# >= 13

# generate objective vector for optimization. This corresponds
# (up to a multiplicative constant) to the
# Piatetsky-Shapiro quality function for the target value credit-state==good
objective <- oofos:::compute_objective(dat[indexs,],"class","good")

table(objective)

# objective
# -0.00628930817610063  0.00293255131964809
# 159                  341
#
# 341 good states and 159 bad states


############################################################################
# Classical Subgroup Discovery on a subset of size n=500
############################################################################


# generate gurobi model for classical subgroup discovery
model <- oofos:::optimize_on_context_extents(context=context,objective=objective,binary_variables="all" )
# For later computations it is more effective to treat the constraint matrix
# not as a simple triplet matrix
model$A <- as.matrix(model$A)
# add additional constraints that may help to speed up the optimization
model <- oofos:::add_attr_antiimplications(model)
#compute a preliminary result with timelimit of 10 minits
result_1 <- gurobi::gurobi(model,list(timelimit=60*10))

result_1$objval
# [1] 0.6388889

model$objval <- result_1$objval
test <- oofos::compute_extent_optim_test(model)


# p-value ca. 0.09

D <- as.matrix(dist(context))
#D <- as.matrix(clue::ls_fit_ultrametric(D))
context2 <- get_context_from_distance(D,threshold=20,lambda=0.6,complemented=FALSE)

context3 <- context2

p <- ncol(context2)
idxs <- NULL
for(k in (1:p)){
  d <- rep(0,nrow(context))
  for(l in seq_len(nrow(context))){
    d[l] <- mean(context2[,k]!= context[,l])
  }
  idxs <- c(idxs,which.min(d))
  }

idxs <- unique(idxs)

#context2 <- cbind(context2,context[,idxs])
context2 <- context[,idxs]

#vc_dim <- gurobi::gurobi(oofos::compute_extent_vc_dimension(context2))

context2 <- t(unique(t(context2)))
dim(context2)
model <- oofos:::optimize_on_context_extents(context=context2,objective=objective,binary_variables="all" )

model$A <- as.matrix(model$A)
# add additional constraints that may help to speed up the optimization
model <- oofos:::add_attr_antiimplications(model)
#compute a preliminary result with timelimit of 10 minits
result_2 <- gurobi::gurobi(model,list(timelimit=60*10))


model$objval <- result_2$objval
test <- oofos::compute_extent_optim_test(model)













result_1$runtime
# [1] 601.01

# save result
saveRDS(result_1,"results_credit_data/result_1.RDS")



# use semioptimal solution to tighten the search space via optimistic estimates
# (in the style of M. Boley and H. Grosskreutz. Non-redundant subgroup discovery us-
# ing a closure system.)
model_2 <- oofos:::add_sos_constraints(model,result_1$objval)
# add semioptimal solution as a start vector
model_2$start <- round(result_1$x,2)
result_2 <- gurobi::gurobi(model_2,list(PoolSolutions=20,PoolSearchMode=2))


result_2$objval

# [1] 0.3572364

#overall runtie in hours
(result_1$runtime  + result_2$runtime)/3600
# [1] 9.11968

saveRDS(result_2,"result_2.RDS")

############################################################################
# Smple splitting test
############################################################################

extent <- round(result_2$x[(1:500)],2)
intent <- round(result_2$x[-(1:500)],2)
extent_test <- oofos:::compute_phi(intent,test_context)
test_objective <- oofos:::compute_objective(dat[-indexs,],"class","good")
observed_test_statistic <- sum(test_objective*extent_test)
observed_test_statistic

n_rep <-1000000
h0 <- rep(0, n_rep)
set.seed(1234567)
for(k in (1:n_rep)){
  h0[k] <- sum(sample(test_objective)*extent_test)
}
plot(ecdf(h0))
mean(h0 >= observed_test_statistic)

#0

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#  END OF CURATED PART

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

############################################################################
# Classical Subgroup Discovery on a subset of size n=500
############################################################################
# Laufzeit ca 40 min
#unter H0: Laufzeit ca. > 83864s
# unter H0 mit gleichem Trick: 13040.84 s
# mit bestem SOS constraint: 3274s

result_2$objval
#[1] 0.1326288
set.seed(1234567)
indexs <- sample((1:1000),size=500)
sampled_context <- context[indexs,]

sampled_objective <- oofos:::compute_objective(dat[indexs,],"class","good")

model_sampled <- oofos:::optimize_on_context_extents(context=sampled_context,objective=sampled_objective,binary_variables="all" )
model_2_sampled <- oofos:::add_attr_antiimplications(model_sampled)
model_2_sampled <- oofos:::add_sos_constraints(model_2_sampled,0.278)
result_sampled <- gurobi::gurobi(model_2_sampled,list(timelimit=60*40))


############################################################################
# Starshaped Subgroup Discovery with angle based stylized betweenness
############################################################################

indexs_numerical_variables <- c(2,5,8,11,13,16,18)

absb <- get_absb(as.matrix(dat[,indexs_numerical_variables]))

saveRDS(absb,"results_credit_data/absb.RDS")

discovery_absb <- oofos:::discover_starshaped_subgroups(stylized_betweenness=absb,objective=objective,local_vc_dimension=100,params=list(outputflag=1))
