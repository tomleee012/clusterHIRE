setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE")
setwd("./HIREewas/src")
system("R CMD SHLIB HIRE.c")
dyn.load('HIRE.dll')

###################################################################################################################################
#function 1
###################################################################################################################################

HIRE <- function(Ometh, X, num_celltype, tol = 10^(-5), num_iter=1000, alpha=0.01){

	#generate samples from a Dirichlet distribution
	rDirichlet <- function(alpha_vec){
		num <- length(alpha_vec)
		temp <- rgamma(num, shape = alpha_vec, rate = 1)
		return(temp / sum(temp))
	}

	#initialize P_t

	CorDescent <- function(MethMatr, num_celltype, tol = 0.01, showIter = FALSE){
		err0 <- 0
		err1 <- 1000
		m <- nrow(MethMatr) #number of CpG sites
		n <- ncol(MethMatr) #number of samples
		K <- num_celltype
		if(m < n){
			stop("The CpG site number must be larger than the sample number!")
		}
		#initialize P_matr
		P_matr_t <- vapply(seq_len(n), function(i){ rDirichlet(rep(2,K)) }, FUN.VALUE=rep(-1,K))
		while(abs(err1 - err0) >= tol){
			err0 <- err1
			#update U_matr
			Dmat <- 2*P_matr_t%*%t(P_matr_t)
			Amat <- cbind(diag(rep(1,K)), diag(rep(-1,K)))
			bvec <- c(rep(0,K), rep(-1,K))
			U_matr_t <- t( vapply(seq_len(m), function(j){
							dvec <- 2*P_matr_t %*% as.numeric(MethMatr[j,])

							solu <- solve.QP(Dmat, dvec, Amat, bvec)
							solu$solution
					}, FUN.VALUE=rep(-1,K)) )
		
			#update P_matr
			Dmat <- 2*t(U_matr_t) %*% U_matr_t
			Amat <- cbind(matrix(1, K, K), diag(rep(1,K)))
			bvec <- c(rep(1, K), rep(0, K))
			P_matr_t <- vapply(seq_len(n), function(i){
						dvec <- 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
						solu <- solve.QP(Dmat, dvec, Amat, bvec, meq = K)
						solu$solution 
					}, FUN.VALUE=rep(-1,K))
					
			#calculate errors
			err1 <- sum((MethMatr - U_matr_t %*% P_matr_t)^2)
			if(showIter == TRUE){
				message("  ", err1, "\n")
			}
		}
	
		return(list(U=U_matr_t, P=P_matr_t))	
	}

	Initialize <- function(Ometh, num_celltype){
		K <- num_celltype
		sdrow <- apply(Ometh, 1, sd)
		ind <- order(sdrow, decreasing = TRUE)
		m <- nrow(Ometh)
		if(m <= 1000){
			num_cpg_for_init <- m
		}else{
			num_cpg_for_init <- max(3*ncol(Ometh), floor(m/10))
			if(num_cpg_for_init > m){
				num_cpg_for_init <- m
			}
		}

		Ometh_part <- Ometh[ind[seq_len(num_cpg_for_init)],] #select CpG sites with the most num_cpg_for_init variant methylation levels 

		result <- CorDescent(Ometh_part, num_celltype=K, tol = 0.1, showIter = FALSE)
		P_initial <- result$P

		mu_initial <- vapply(seq_len(m), function(j){
				if(K > 2){
					fit <- lm(Ometh[j,]~as.matrix(t(P_initial[-1, ])))
				}else{
					fit <- lm(Ometh[j,]~as.numeric(P_initial[-1, ]))
				}
				tmp <- as.numeric(summary(fit)$coeff[ ,1])
				tmp[-1] <- tmp[1] + tmp[-1]
				tmp
			}, FUN.VALUE = rep(-1,K) )
		return(list(P_initial, mu_initial))
	}

	EmEwasRcallC <- function(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter){
		args <- list("P_init"=as.numeric(P_t), "mu_init"=as.numeric(mu_t), "beta_init"=as.numeric(beta_t),
				"beta_init_dim"=as.integer(dim(beta_t)), "Ometh_r"=as.numeric(Ometh),
				"Ometh_r_dim"=as.integer(dim(Ometh)), "X_r"=as.numeric(X), "X_r_dim"=as.integer(dim(X)),
				"sig_sqTiss_init"=as.numeric(sig_sqTiss_t), "sig_sqErr_init"=as.numeric(sig_sqErr_t),
				"tol_r" = as.numeric(tol), "num_iter" = as.integer(num_iter))
		ret_list <- .Call("EmEwasRcallC", args)		
		return(ret_list)
	}


	m <- nrow(Ometh) #CpG site number
	n <- ncol(Ometh) #sample number
	p <- nrow(X)
	K <- num_celltype
	
	P_t <- matrix(-1, K, n)
	mu_t <- matrix(-1, m, K)
	beta_t <- array(-1, dim=c(m, K, p))
	sig_sqTiss_t <- matrix(-1, m, K)
	sig_sqErr_t <- rep(-1, m)

	init <- Initialize(Ometh, K)
	P_t <- init[[1]]
	mu_t <- t(init[[2]])
	message("  Initialization Done.\n")
	message("  Implementing EM algorithm... \n")
	ret_list <- EmEwasRcallC(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter)
	message("  Done! \n")
	
	message("  Calculating p-values...\n")
	tmp <- NULL
	for(ell in seq_len(p)){
		tmp <- cbind(tmp, X[ell, ]*t(ret_list$P_t))
	}
	x_matr <- cbind( tmp, t(ret_list$P_t)[, seq(2,K)])
	x_matr <- as.matrix(x_matr)

	pvalues <- t(vapply(seq_len(m), function(j){
					y_vec <- Ometh[j,]
					fit <- lm(y_vec~x_matr)
					summary(fit)$coef[seq(2, (1+p*K)),4]
				}, FUN.VALUE = rep(-1,p*K)))
	message("  Done!\n")
	ret_list[[7]] <- pvalues
	names(ret_list)[7] <- "pvalues"
	names(ret_list)[6] <- "pBIC"
	d <- (K-1)*n + m*(1+2*K+p*K) #number of parameters
	d0 <- sum(ret_list$pvalues > alpha/(p*m*K)) #number of zero parameters
	ret_list[[6]] <- ret_list[[6]] - log(n)*d + log(n)*(d-d0) 
	
	return(ret_list)
}


###################################################################################################################################
#function 2
###################################################################################################################################

riskCpGpattern <- function(pval_matr, main_title="Detected association pattern", hc_row_ind = FALSE){
	m <- nrow(pval_matr)
	K <- ncol(pval_matr)

	colors_func <- colorRampPalette(c("grey", "black"))
	colors <- colors_func(10)
	#for zero pvalues 
	pval_matr[pval_matr == 0] <- 10^(-100) 

	if(hc_row_ind == FALSE){
		heatmap.2(-log10(pval_matr), col = colors, scale = "none",
	 		key = TRUE, key.xlab = "-log10(p-value)",
	 		Colv=FALSE,Rowv=FALSE, density.info = "none", trace = "none", dendrogram="none", 
	 		ylab = paste0(m, " CpG sites"), xlab = "Cell types", margins = c(5,5),  
			main = main_title,labCol=paste0("Cell type ", seq_len(K)),
			labRow=FALSE, cexRow=0.9, srtCol=0, cexCol=1, adjCol=c(NA,1), colsep=seq_len(K))
	}else{
		heatmap.2(-log10(pval_matr), col = colors, scale = "none",
	 		key = TRUE, key.xlab = "-log10(p-value)",
	 		Colv=FALSE,Rowv=TRUE, density.info = "none", trace = "none", dendrogram="row", 
	 		ylab = paste0(m, " CpG sites"), xlab = "Cell types", margins = c(5,5),  
			main = main_title,labCol=paste0("Cell type ", seq_len(K)),
			labRow=FALSE, cexRow=0.9, srtCol=0, cexCol=1, adjCol=c(NA,1), colsep=seq_len(K))		
	}
}

