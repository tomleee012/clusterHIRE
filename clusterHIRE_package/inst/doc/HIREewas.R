## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex()

## ----data_preparation1-----------------------------------------------------
###############################################################
#Generate the EWAS data 
###############################################################

set.seed(123)
###define a function to draw samples from a Dirichlet distribution
rDirichlet <- function(alpha_vec){
	num <- length(alpha_vec)
	temp <- rgamma(num, shape = alpha_vec, rate = 1)
	return(temp / sum(temp))
}

n <- 180     #number of samples 
n1 <- 60    #number of controls
n2 <- 120    #number of cases 

m <- 2000   #number of CpG sites
K <- 3       #underlying cell type number
	
###simulate methylation baseline profiles
#assume cell type 1 and cell type 2 are from the same lineage
#cell type 1
methy1 <- rbeta(m,3,6)
#cell type 2 
methy2 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(1:m, m/5) 
methy2[ind] <- rbeta(length(ind),3,6)

#cell type 3
methy3 <- rbeta(m,3,6)
mu <- cbind(methy1, methy2, methy3)

#number of covariates
p <- 2

###simulate covariates / phenotype (disease status and age)
X <- rbind(c(rep(0, n1),rep(1, n2)), runif(n, min=20, max=50))

###simulate phenotype effects
beta <- array(0, dim=c(m,K,p))

#control vs case
m_common <- 10
max_signal <- 0.15
min_signal <- 0.07

#we allow different signs and magnitudes
signs <- sample(c(-1,1), m_common*K, replace=TRUE)
beta[1:m_common,1:K,1] <- signs * runif(m_common*K, min=min_signal, max=max_signal)

m_seperate <- 10
signs <- sample(c(-1,1), m_seperate*2, replace=TRUE)
beta[m_common+(1:m_seperate),1:2,1] <- signs * 
					runif(m_seperate*2, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), m_seperate, replace=TRUE)
beta[m_common+m_seperate+(1:m_seperate),K,1] <- signs * 
					runif(m_seperate, min=min_signal, max=max_signal)

#age
base <- 20
m_common <- 10
max_signal <- 0.015
min_signal <- 0.007
signs <- sample(c(-1,1), m_common*K, replace=TRUE)
beta[base+1:m_common,1:K,2] <- signs * 
					runif(m_common*K, min=min_signal, max=max_signal)

m_seperate <- 10
signs <- sample(c(-1,1), m_seperate*2, replace=TRUE)
beta[base+m_common+(1:m_seperate),1:2,2] <- signs * 
					runif(m_seperate*2, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), m_seperate, replace=TRUE)
beta[base+m_common+m_seperate+(1:m_seperate),K,2] <- signs * 
					runif(m_seperate, min=min_signal, max=max_signal)

###generate the cellular compositions 
P <- sapply(1:n, function(i){
				if(X[1,i]==0){ #if control
					rDirichlet(c(4,4, 2+X[2,i]/10))
				}else{
					rDirichlet(c(4,4, 5+X[2,i]/10))
				}	
			})

###generate the observed methylation profiles 
Ometh <- NULL
for(i in 1:n){
	utmp <- t(sapply(1:m, function(j){
					tmp1 <- colSums(X[ ,i] * t(beta[j, , ]))
					rnorm(K,mean=mu[j, ]+tmp1,sd=0.01)
				}))
	tmp2 <- colSums(P[ ,i] * t(utmp))
	Ometh <- cbind(Ometh, tmp2 + rnorm(m, sd = 0.01))
}

#constrain methylation values between 0 and 1
Ometh[Ometh > 1] <- 1

Ometh[Ometh < 0] <- 0

## ----data-preparation2-----------------------------------------------------
#the class of the methylation matrix
class(Ometh)

#the values in the methylation matrix
head(Ometh[,1:6])

#the class of the covariate matrix
class(X)

#the values in the covariate matrix
X[ ,1:6]

## ----model1----------------------------------------------------------------
library(HIREewas)
ret_list <- HIRE(Ometh, X, num_celltype=K, tol=10^(-5), num_iter=1000, alpha=0.01)

## ----model2----------------------------------------------------------------
# the class of ret_list
class(ret_list)

#the estimated cellular compositions 
ret_list$P_t[ ,1:6]

#the estimated cell-type-specific methylation baseline profiles
head(ret_list$mu_t)

#the estimated phenotype effects
head(ret_list$beta_t)

#the penalized BIC value
ret_list$pBIC

#the estimated p-values to claim whether a CpG site is at risk
#in some cell type for a covariate

#p value matrix for case/control
head(ret_list$pvalues[ ,1:3])

#p value matrix for age
head(ret_list$pvalues[ ,4:6]) 

## --------------------------------------------------------------------------
#estimated cell compositions vs the truth
par(mfrow=c(1,3))
plot(ret_list$P_t[2, ], P[1, ], xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red")

plot(ret_list$P_t[1, ], P[2, ], xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red")

plot(ret_list$P_t[3, ], P[3, ], xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1, col="red")

## --------------------------------------------------------------------------
riskCpGpattern(ret_list$pvalues[1:100, K+c(2,1,3)], 
		main_title="Detected association pattern\n with age", hc_row_ind = FALSE)

