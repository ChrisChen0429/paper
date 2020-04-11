library(GDINA)
library(data.table)
library(truncdist)
library(MCMCpack)

temp <- sim10GDINA
x <- temp$simdat
Q <- temp$simQ
par <- temp$simItempar


saturated_gibbs <- function(x,Q,N_ITR){ # N_IRT:= number of iteration
	N <- dim(x)[1]
	J <- dim(x)[2]
	K <- dim(Q)[2]
  
	comp <- function(v1,v2){
	  # 1: p1 >= p2; -1: p1 < p2; 0: unconstrained.
	  if (sum(v1^v2) == length(v1)) return (1)
	    else if (sum(v2^v1) == length(v1)) return (-1)
	  else
	  return (0)
	}
	
	BinToDec <- function(x){
	  return (strtoi(do.call(paste, c(as.list(as.character(x)), sep = "")),base=2))
	}
	
	####  create the attribute profile parameter space (dichotomous) ### 
	attribute_level <- list()
	for (k in 1:K) attribute_level[[k]] <- 0:1 ## if ploytomous, change the level
	alpha_matrix <- as.matrix(do.call(expand.grid,attribute_level))
	#####################################################################
	
	constraint_matrix_list <- list()
	  # constraint_matrix_list is used for set up the monotonicity constraint
	p_vec_matrix <- matrix(NA,nrow=J,ncol=dim(alpha_matrix)[1])
	  # p_vec_matrix define for each item (row), what is the unique aggrement pattern each attribute profile belongs to
	
	for (itr in 1:J){
	  current_temp <- Q[itr,] * t(alpha_matrix) # aggreement of q vector with attribute profile 
	  v_temp <- apply(current_temp,2,BinToDec)
	        # v_temp: the index of attibute pattern for the current_temp
	  dt <- as.data.table(v_temp)[,list(list(.I)),by=v_temp]$V1
	        # dt: list of index for different group of idenfiable attribute pattern
	  M <- length(dt)
	        # number of the identifiable attribute profile group (same as length(unique(v_temp)))
	  
	  constraint_matrix <- matrix(NA,nrow=M,ncol=M)
	  for (i in 1:(M-1)){
	      for (j in (i+1):M){
	        constraint_matrix[i,j] <- comp(current_temp[,dt[[i]][1]],current_temp[,dt[[j]][1]])
	        constraint_matrix[j,i] <- -1 * constraint_matrix[i,j]
	      }	
	    }
	    
	  
	  p_vec <- rep(NA,length(v_temp))
	  for (i in 1:length(dt)){
	    p_vec[dt[[i]]] <- i
	  }
	  constraint_matrix_list[[itr]] <- constraint_matrix
	  p_vec_matrix[itr,] <- p_vec
	}
	


	
	#################### initialize the chain ###############
	#N_ITR <- 100
	
	##### marginal probability of attribute profile ####
	post_pi <- matrix(NA,nrow=N_ITR,ncol=dim(alpha_matrix)[1])
	post_pi[1,] <- rdirichlet(1,rep(1,8))
	###################################################
	
	##### examinee's attribute profile index #########
	post_c <- matrix(NA,nrow=N_ITR,ncol=N) 
	post_c[1,] <- sample(c(1:dim(alpha_matrix)[1]),N,replace=TRUE,post_pi[1,]) 
  #################################################
  
	##### correct response probability  with monotonicity #####
	post_p <- list() # list of item 
	for (i in 1:J){
		temp <- max(p_vec_matrix[i,]) # number of aggrement pattern 
		post_p[[i]] <- matrix(NA,nrow=N_ITR,ncol=temp)
	}
	
	for (i in 1:length(post_p)){
		for (j in 1:dim(post_p[[i]])[2]){
			lower <- -Inf
			upper <- Inf
			ind_min <- which(constraint_matrix_list[[i]][j,] == 1) # p1 >= p2
			ind_max <- which(constraint_matrix_list[[i]][j,] == -1) # p1 < p2
			if (length(ind_min) > 0)
				lower <- max(post_p[[i]][1,ind_min])
			if (length(ind_max) > 0)
				upper <- min(post_p[[i]][1,ind_max])
			if(is.na(lower)) lower <- -Inf
			if(is.na(upper)) upper <- Inf
			post_p[[i]][1,j] <- rtrunc(1,"beta",a=lower,b=upper,1,1)
		}
	}
	#########################################################

	n_cat <- apply(p_vec_matrix,1,max)

	for (n_itr in 2:N_ITR){
	
	######## get attribute profile (t -1) ##############
		dt <- apply(outer(post_c[n_itr-1,],c(1:dim(alpha_matrix)[1]),"=="),2,f <- function(input) which(input == 1))
		
	######## correct response probability  with monotonicity #####################
		for (j in 1:J){
			for (c in 1:n_cat[j]){
				ind_min <- which(constraint_matrix_list[[j]][c,] == 1)
				ind_max <- which(constraint_matrix_list[[j]][c,] == -1)
				lower <- -Inf
				upper <- Inf
				if(length(ind_min) > 0){
					for (m in ind_min){
						if (m < c){
							if (lower < post_p[[j]][n_itr,m])
								lower <- post_p[[j]][n_itr,m]
						}
						else{
							if (lower < post_p[[j]][n_itr-1,m])
								lower <- post_p[[j]][n_itr-1,m]
						}
					}
				}
				if(length(ind_max) > 0){
					for (m in ind_max){
						if (m < c){
							if (upper > post_p[[j]][n_itr,m])
								upper <- post_p[[j]][n_itr,m]
						}
						else{
							if (upper > post_p[[j]][n_itr-1,m])
								upper <- post_p[[j]][n_itr-1,m]
						}
					}
				}
				# if(is.na(lower)) lower = -Inf
				# if(is.na(upper)) upper = Inf
				# error handling for numerical issue with sampling from the truncated beta distribution where the intervals are in the low density regions.
				
				n_c <- 0
				x_sum_c <- 0
				ind <- which(p_vec_matrix[j,]==c)
				for (k in ind){
				  n_c <- n_c + length(dt[[k]])
				  x_sum_c <- x_sum_c + sum(x[dt[[k]],j])
				post_p[[j]][n_itr,c] <- tryCatch(rtrunc(1,"beta",a=lower,b=upper,1+x_sum_c,1+n_c-x_sum_c),error = function(e){rtrunc(1,"beta",a=lower,b=upper,1,1)})
				
				

				}
			}		
		}
    ################################################
		
		########## update class memberships ############################
		current_p_matrix <- matrix(NA,nrow=dim(p_vec_matrix)[1],ncol=dim(p_vec_matrix)[2])
		for (i in 1:dim(current_p_matrix)[1]){
			for (j in 1:dim(current_p_matrix)[2]){
				current_p_matrix[i,j] <- post_p[[i]][n_itr,p_vec_matrix[i,j]]
			}
		}
		for (i in 1:N){
			logp <- apply(x[i,] * log(current_p_matrix) + (1-x[i,]) * log(1-current_p_matrix),2,sum) + log(post_pi[n_itr-1,])
			prob_vec <- exp(logp)
			post_c[n_itr,i] <- sample(1:length(prob_vec),1,replace=TRUE,prob_vec)
		}

		#update prior probabilities for class assignment
		n_c <- apply(outer(post_c[n_itr,],c(1:dim(alpha_matrix)[1]),"=="),2,sum)
		post_pi[n_itr,] <- rdirichlet(1,n_c+1)
	}
	out <- list(post_p,post_c,post_pi,p_vec_matrix,constraint_matrix_list)
	return(out)
}


# temp <- saturated_gibbs(x,Q,100)

# #compute posterior means
# plot(post_p[[1]][,1]~dis)
# plot(post_p[[1]][,1] ~ post_p[[1]][,2])

# plot(post_p[[1]][,1]~c(1:3000),type="l")
# library(coda)
# item_1 <- mcmc(post_p[[1]],start=501,end=3000,thin=3)
# plot(item_1)



























library(CDM)
x <- data.ecpe$data
x <- as.matrix(x)
x <- x[,-1]
Q <- data.ecpe$q.matrix
Q <- as.matrix(Q)
start.time <- Sys.time()
out <- saturated_gibbs(x,Q,100)
save(out,file=("out.Rdata"))
end.time <- Sys.time()
end.time - start.time
item_1 <- mcmc(out[[1]][[1]],start=50,end=100,thin=3)
plot(item_1)
effectiveSize(item_1)

post_1 <- out[[1]][[1]][50:100,]
dim(post_1)
colMeans(post_1)

post_pi <- out[[3]][50:100,]
colMeans(post_pi)
mod2 <- CDM::gdina( data.ecpe$data[,-1], q.matrix= data.ecpe$q.matrix , linkfct="logit", rule="GDINA2")
summary(mod2)


#transform back to the log-linear model parameters
# mat_coef_2 <- matrix(c(1,0,0,0,
# 					   1,1,0,0,
# 					   1,0,1,0,
# 					   1,1,1,1),nrow=4,byrow=TRUE)
# mat_ceof_1 <- matrix(c(1,0,
# 					   1,1),nrow=2,byrow=TRUE)
# inv_mat_coef_2 <- solve(mat_coef_2)
# inv_mat_coef_1 <- solve(mat_ceof_1)

mat_coef <- matrix(c(1,0,0,0,0,0,0,0,
					 1,1,0,0,0,0,0,0,
					 1,0,1,0,0,0,0,0,
					 1,1,1,0,1,0,0,0,
					 1,0,0,1,0,0,0,0,
					 1,1,0,1,0,1,0,0,
					 1,0,1,1,0,0,1,0,
					 1,1,1,1,1,1,1,1),nrow=8,byrow=TRUE)
inv_mat_coef <- solve(mat_coef)

post_lambda <- list()
lambda_means <- NULL
lambda_sd <- NULL
for (item in 1:length(out[[1]])){
	matrix_logit_p <- apply(boot::logit(out[[1]][[item]][1001:5000,]),1,f <- function(vec) return(vec[out[[4]][item,]]))
	post_lambda[[item]] <- t(inv_mat_coef %*% matrix_logit_p)
	lambda_means <- round(rbind(lambda_means,colMeans(post_lambda[[item]])),2)
	lambda_sd <- round(rbind(lambda_sd,apply(post_lambda[[item]],2,sd)),2)
}
ind_column <- c(1,4,3,2,7,6,5,8)
lambda_means <- lambda_means[,ind_column]
lambda_sd <- lambda_sd[,ind_column]
lambda_means
lambda_sd

####################### making a latex table to report parameter estiamtes ###########
vec <- rep(NA,length(lambda_means))
for (i in 1:length(lambda_means)){
	if (lambda_means[i] == 0 & lambda_sd[i] == 0){
		vec[i] <- ""
	}
	else
		vec[i] <- paste0(lambda_means[i],"(",lambda_sd[i],")")
}
output_table <- matrix(vec,nrow=28,byrow=FALSE)
output_table <- output_table[,-8]
output_table <- cbind(c(1:28),output_table)
names_col <- c("item",expression(lambda[0]),expression(lambda[1]),expression(lambda[2]),expression(lambda[3]),
	expression(lambda[12]),expression(lambda[13]),expression(lambda[23]))
colnames(output_table) <-  names_col
stargazer::stargazer(output_table,digits=2)
##################################################

item <- 20
matrix_logit_p <- apply(boot::logit(out[[1]][[item]][1001:5000,]),1,f <- function(vec) return(vec[out[[4]][item,]]))
post_lambda <- t(inv_mat_coef %*% matrix_logit_p)
post_lambda <- post_lambda[,ind_column]
dim(post_lambda)
png("/users/xiangliu/desktop/bookchapter/item_20.png")
f <- MASS::kde2d(post_lambda[,2],post_lambda[,6],lims=c(-1.0,3.0,-1.0,2.0))
filled.contour(f,color.palette=gray.colors,xlab=expression(lambda[1]),ylab=expression(lambda[13]))
dev.off()

item <- 1
matrix_logit_p <- apply(boot::logit(out[[1]][[item]][1001:5000,]),1,f <- function(vec) return(vec[out[[4]][item,]]))
post_lambda <- t(inv_mat_coef %*% matrix_logit_p)
post_lambda <- post_lambda[,ind_column]
png("/users/xiangliu/desktop/bookchapter/item_1.png")
f <- MASS::kde2d(post_lambda[,2],post_lambda[,3],lims=c(0,2.0,0,2.0))
filled.contour(f,color.palette=gray.colors,xlab=expression(lambda[1]),ylab=expression(lambda[2]))
dev.off()

item <- 11
matrix_logit_p <- apply(boot::logit(out[[1]][[item]][1001:5000,]),1,f <- function(vec) return(vec[out[[4]][item,]]))
post_lambda <- t(inv_mat_coef %*% matrix_logit_p)
post_lambda <- post_lambda[,ind_column]
png("/users/xiangliu/desktop/bookchapter/item_11.png")
f <- MASS::kde2d(post_lambda[,2],post_lambda[,4],lims=c(0,3.0,0,2.0))
filled.contour(f,color.palette=gray.colors,xlab=expression(lambda[1]),ylab=expression(lambda[3]))
dev.off()




library(boot)


#colMeans(post_lambda[2000:5000,])
#plot(post_lambda[1000:5000,3] ~ post_lambda[1000:5000,2])
library(MASS)
ind <- c(1,3,7,11,12,16,17,20,21)
for (item in ind){
	n_req_attributes <- sum(Q[item,] != 0)
	if (n_req_attributes == 2){
		mat_inv <- inv_mat_coef_2
		post_lambda <- matrix(NA,nrow=5000,4)
	}
	if (n_req_attributes == 1){
		mat_inv <- inv_mat_coef_1
		post_lambda <- matrix(NA,nrow=5000,2)
	}
	for (i in 1:dim(post_lambda)[1]){
		post_lambda[i,] <- mat_inv %*% t(t(logit(out[[1]][[item]][i,])))
	}
	file_name <- paste0("/users/xiangliu/desktop/bookchapter/item_",item,".png")
	png(file_name)
	f <- kde2d(post_lambda[2000:5000,2],post_lambda[2000:5000,3])
	filled.contour(f,color.palette=gray.colors,xlab=expression(lambda[1]),ylab=expression(lambda[2]))
	dev.off()
}

vec_mean <- apply(post_lambda[1000:5000,],2,mean)
vec_sd <- apply(post_lambda[1000:5000,],2,sd)
table_lambda <- cbind(vec_mean,vec_sd)

post_class_size <- matrix(NA,nrow=5000,ncol=8)

for (i in 1:5000){
	post_class_size[i,] <- apply(outer(out[[2]][i,],c(1:8[1]),"=="), 2, f <- function(input) length(which(input == 1))) / dim(x)[1]
}
colMeans(post_class_size)
library(coda)
item_1 <- mcmc(post_class_size,start=501,end=3000,thin=3)
plot(item_1)

f <- kde2d(post_class_size[2000:5000,6],post_class_size[2000:5000,7])
filled.contour(f,color.palette=gray.colors)


ind <- c(3,5,6,7)
class_size_pca <- prcomp(post_class_size[2000:5000,ind])
plot(class_size_pca,type="l")
dat <- predict(class_size_pca)
f <- kde2d(dat[,1],dat[,2])
filled.contour(f,color.palette=gray.colors)




