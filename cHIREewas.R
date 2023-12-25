library(quadprog)
library(reshape2)
library(ggplot2)
library(foreach)
library(doParallel)
library(reshape)
rm(list = ls())
setwd('~./')

##############################
# Implementation (fast)
##############################
N = 2000
M = 10000
K = 4
Q = 3
species = 'Humans'
ver = 4

# # change dir
# setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
# setwd(paste0('./Q', Q))

# change dir
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE")
subDir = paste0('N', N, 'K', K, 'Q', Q, 'V', ver)

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}
load('data.RData')
library(cHIRE2)
time0 = proc.time()
ret_list=cHIRE2::cHIRE(DGP$Ometh, DGP$Covariates, Init$EZ, Init$Pi_t, Init$Mu, Init$Beta,
                       Init$Gamma, Init$P_t, iter=2000, tol = 1e-20, core=20)
rm(Init)
print(proc.time() - time0)
save.image("cHIRE2_.RData")

##############################
# Implemenetation (Optimal)
##############################
N = 2000
M = 10000
K = 4
Q = 3
species = 'Humans'
ver = 4

# # change dir
# setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
# setwd(paste0('./Q', Q))

# change dir
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE")
subDir = paste0('N', N, 'K', K, 'Q', Q, 'V', ver)

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}
# implementation
library(cHIRE4)
for(q in 2){
  setwd(paste0('Q', q))
  for(k in 2:5){
    setwd(paste0('K', k))
    load('data.RData')
    time0 = proc.time()
    ret_list=cHIRE4::cHIRE(DGP$Ometh, DGP$Covariates[,-1], Init$EZ, Init$Pi_t, Init$Mu, Init$Beta,
                           Init$Gamma, Init$P_t, iter=1000, tol = 1e-20, core=20)
    print(proc.time() - time0)
    save.image("cHIRE4_.RData")
    cat("\014")
    rm(list = ls())
    gc()
    setwd('../')
  }
  setwd('../')
}

##############################
# self-written functions
##############################
# progress bar
progress = function(i,I,freq=10,t0=NA){
  if(is.na(t0)){
    if(i==2){cat("### Progress ###\n");cat(paste(rep("-",100),collapse=""));cat("\n")}
    if(i%%(I/100)<1){cat(">")}
  }else{
    if(i%%freq==0){
      cat("\n",paste(rep("=",60),collapse=""));cat("\n")
      cat(paste("   i=",i," at ",Sys.time(),"\n",sep=""))
      cat(paste("   Time Remaining: ",
                as.numeric((Sys.time()-t0)/i*(I-i),units="hours")," hours",sep=""))
      cat("\n",paste(rep("=",60),collapse=""))
    }
  }
}

# generate samples from a Dirichlet distribution
rDirichlet = function(alpha_vec, seed){
  set.seed(seed)
  num = length(alpha_vec)
  temp = rgamma(num, shape = alpha_vec, rate = 1)
  return(temp / sum(temp))
}

#############################
# Data Generating Process and initiation
#############################
MPCS.DGP = function(N, M, K, Q, species, prt = 0, ver, seed = 1){
  if(Q > K){
    message('Error! The number of celltypes must not be less than the number of cancer subtypes')
    return()
  }

  repeat{
    print(seed)
    set.seed(seed)

    cancer.subtype.prop = switch(Q, 1, c(.7, .3), c(.5, .3, .2), c(.5, .2, .2, .1),
                                 c(.4, .2, .15, .15, .1))

    message(paste0('The cancer proportions are ', cancer.subtype.prop, '\n'))

    cancer.identity = rep(1, N*cancer.subtype.prop[1])
    for(q in 2:Q)
      cancer.identity = c(cancer.identity, rep(q, N*cancer.subtype.prop[q]))

    disease_covariates = array(0, dim = c(Q,N))
    for(n in 1:N)
      disease_covariates[cancer.identity[n],n] = 1
    
    permute.cell.type.prop = list()
    permute.cell.type.prop[[1]] = rep(1, K)
    permute.cell.type.prop[[2]] = (1:K)^3
    permute.cell.type.prop[[3]] = (K:1)^3
    permute.cell.type.prop[[4]] = `if`(K %% 2 == 0, c( (1:(K/2))^3, ((K/2):1)^3 ),
                                       c( (1:floor(K/2))^3, (floor(K/2+1):1)^3 ))
    permute.cell.type.prop[[5]] = `if`(K %% 2 == 0, c( ((K/2):1)^3, (1:(K/2))^3 ),
                                       c( (ceiling(K/2):1)^3, (2:floor(K/2+1))^3 ))
    permute.cell.type.prop[[6]] = c(1, rep(3, K-1))
    permute.cell.type.prop[[7]] = c(rep(3, K-1), 1)
    
    ###simulate covariates / phenotype
    # gender
    num_female = N * .6
    genders = c(rep(0, num_female), rep(1, N - num_female))
    genders = sample(genders)
    # age
    ages = runif(N, min=20, max=70)
    # normalize the covariate age
    min_max_normal = function(arr){return((arr - min(arr)) / (max(arr) - min(arr)))}
    ages = min_max_normal(ages)
    # smokers' status
    num_cur = N * .2
    num_non = N - num_cur
    smokers = c(rep(0, num_non), rep(1, num_cur))
    smokers = sample(smokers)
    
    P = 3 # number of covariates
    Covariates = cbind(genders, ages, smokers)
    # Covariates = cbind(genders, 1 - genders)
    # Covariates = cbind(ages, rep(0, N))

    cell.type.prop = array(NA, dim = c(K,N))
    # for(q in 1:Q)
    #   permute.cell.type.prop[[q]] = sample(seq(1,K^Q), K, replace = T)
    cell.type.prop = vapply(seq_len(N),
           function(n){ 
             # temp = permute.cell.type.prop[[cancer.identity[n]]]
             temp = permute.cell.type.prop[[1]]
             k_tilde = sample(K, Q-1)
             for(k in k_tilde){
               temp[k] = temp[k] + 5 * cancer.identity[n]
             }
             
             temp[K-2] = temp[K-2] + 3 * Covariates[n,1]
             temp[K-1] = temp[K-1] + 4 * Covariates[n,2]
             temp[K] = temp[K] + 5 * Covariates[n,3]
             rDirichlet(temp, seed + n) 
             },
           FUN.VALUE=rep(-1,K))
    
    Beta = array(0, dim = c(M, K, P))
    M_common = min(M/20, 300);    M_seperate = min(M/20, 400)
    base = M_common + M_seperate*K;    max_sig = 0.05;    min_sig = 0.01
    for(p in 1:P){
      # max_signal = min(max_sig + 0.1 * (q-1), 1)
      # min_signal = min(min_sig + 0.1 * (q-1), 1)
      max_signal = max_sig
      min_signal = min_sig
      
      signs = sample(c(-1,1), M_common*K, replace=TRUE)
      ind = M - M_common*p+ 1:M_common
      Beta[ind,,p] = signs * runif(M_common*K, min=min_signal, max=max_signal)

      for(k in 1:K){
        signs = sample(c(-1,1), M_seperate, replace=TRUE)
        ind2 = M - base*p + M_common + M_seperate*(k-1) + (1:M_seperate)
        Beta[ind2,k,p] = signs * runif(M_seperate, min=min_signal, max=max_signal)
      }
      
    }
    
    for(p in 1:P){
      tmp = Beta[,,p]
      colnames(tmp) = c(paste0("Underlying Cell type ", 1:K))
      tmp = melt(tmp)
      colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
      p = ggplot(tmp, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue",
                             midpoint = 0,
                             mid = "#fcfcfc",
                             high = "red",
                             space="Lab") +
        theme(axis.text.x = element_text(angle = 90, face = "bold", size = 12),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        theme(axis.title = element_text(size=32,face="bold"),
              legend.text = element_text(size=20, face="bold"),
              legend.title = element_text(size = 20,face="bold"))
      print(p)
      Sys.sleep(1)
    }
    
    #Sparsity is common
    Gamma = array(0, dim=c(M,K,Q))
    M_common = min(M/20, 300);    M_seperate = min(M/20, 400)
    base = M_common + M_seperate*K;    max_sig = 0.3;    min_sig = 0.2
    
    for(q in 2:Q){
      # k_tilde = sort(sample(K,Q-1))
        
      # max_signal = min(max_sig + 0.1 * (q-1), 1)
      # min_signal = min(min_sig + 0.1 * (q-1), 1)
      max_signal = max_sig
      min_signal = min_sig

      signs = sample(c(-1,1), M_common*K, replace=TRUE)
      ind = base*(q-2) + 1:M_common
      Gamma[ind,,q] = signs * runif(M_common*K, min=min_signal, max=max_signal)
      # Gamma[ind,,q] = ifelse(Gamma[ind,,q] < -.3, -.3, Gamma[ind,,q])

      for(k in 1:K){
        signs = sample(c(-1,1), M_seperate, replace=TRUE)
        ind2 = base*(q-2) + M_common + M_seperate*(k-1) + (1:M_seperate)
        Gamma[ind2,k,q] = signs * runif(M_seperate, min=min_signal, max=max_signal)
        # Gamma[ind2,k,q] = ifelse(Gamma[ind2,k,q] < -.3, -.3, Gamma[ind2,k,q])
      }

    }

    for(q in 2:Q){
      tmp = Gamma[,,q]
      colnames(tmp) = c(paste0("Underlying Cell type ", 1:K))
      tmp = melt(tmp)
      colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
      p = ggplot(tmp, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue",
                             midpoint = 0,
                             mid = "#fcfcfc",
                             high = "red",
                             space="Lab") +
        theme(axis.text.x = element_text(angle = 90, face = "bold", size = 12),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) +
        theme(axis.title = element_text(size=32,face="bold"),
              legend.text = element_text(size=20, face="bold"),
              legend.title = element_text(size = 20,face="bold"))
      print(p)
      Sys.sleep(1)
    }

    #Baseline Methylation Profile (Between 0 and 1)
    Mu = array(0, dim = c(M, K))
    for(k in 1:K){
      Mu[,K-k+1] = rbeta(M, K + 1 - k, 2*K)
      # Mu[,k] = rbeta(M, 3, 6)
    }

    # if(K >= 3){
    #   Mu[,2] = Mu[,1] + rnorm(M, sd=0.01)
    #   ind = sample(1:M, M/5)
    #   Mu[ind,2] = rbeta(length(ind),3,6)
    # }
    #
    # if(K >= 5){
    #   Mu[,4] = Mu[,3] + rnorm(M, sd=0.01)
    #   ind = sample(1:M, M/5)
    #   Mu[ind,4] = rbeta(length(ind),3,6)
    # }

    tmp = Mu
    colnames(tmp) = c(paste0("Underlying Cell type ", 1:K))
    tmp = melt(tmp)
    colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
    p = ggplot(tmp, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue",
                           midpoint = 0,
                           mid = "#fcfcfc",
                           high = "red",
                           space="Lab") +
      theme(axis.text.x = element_text(angle = 90, face = "bold", size = 12),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      theme(axis.title = element_text(size=32,face="bold"),
            legend.text = element_text(size=20, face="bold"),
            legend.title = element_text(size = 20,face="bold"))
    print(p)
    Sys.sleep(1)

    # Noise
    Sigma_MK = array(.01, dim = c(M,K))
    Sigma_EpsilonM = rep(.01, M)

    # U = array(NA, dim = c(N,M,K)); Ometh = array(NA, dim = c(M,N))
    # for(n in 1:N){
    #   for(m in 1:M){
    #     for(k in 1:K){
    #       temp = 0
    #       for(p in 1:P){
    #         temp = temp + Beta[m,k,p] * Covariates[n,p]
    #       }
    #       U[n,m,k] = rnorm(1, Mu[m,k]+temp+Gamma[m,k,cancer.identity[n]], Sigma_MK[m,k])
    #     }
    #     Ometh[m,n] = rnorm(1, U[n,m,] %*% cell.type.prop[,n], Sigma_EpsilonM[m])
    #   }
    #   progress(n,N)
    # }
    
    Ometh = foreach(n = seq_len(N), .combine = 'cbind') %dopar% {
      temp = t(vapply(seq_len(M), function(m){
        vapply(seq_len(K), function(k){
          temp = 0
          for(p in 1:P){
            temp = temp + Beta[m,k,p] * Covariates[n,p]
          }
          rnorm(1, Mu[m,k]+temp+Gamma[m,k,cancer.identity[n]], Sigma_MK[m,k]) },
          FUN.VALUE=rep(-1,1)) }, FUN.VALUE=rep(-1,K)))
      temp2 = vapply(seq_len(M), function(m){
        rnorm(1, temp[m,] %*% cell.type.prop[,n], Sigma_EpsilonM[m])}, FUN.VALUE=rep(-1,1))
      temp2
    }

    print(paste0('The number of Ometh less than 0 or greater than 1: ', sum(Ometh<0 | Ometh>1)))
    print(paste0('The number of Ometh less than 0: ', sum(Ometh<0)))
    Ometh[Ometh > 1] = 1
    Ometh[Ometh < 0] = 0
    
    row.names(Ometh) = paste0('row', 1:M)
    colnames(Ometh) = paste0('col', 1:N)
    Ometh_melt = melt(Ometh)
    p = ggplot(Ometh_melt, aes(X1, X2)) + geom_tile(aes(fill = value)) + 
      scale_fill_gradient(low = "white", high = "red")
    print(p)
    accept = as.numeric(readline(prompt="Accept it or not: "))
    
    if( accept ){
      break
    }
    seed = seed + 1
  }

  z.nq = t(disease_covariates)
  z.nmq = array(0, dim = c(N,M,Q))
  for(n in 1:N){
    z.nmq[n, , cancer.identity[n]] = 1
  }

  Sigma_MK = array(1, dim = c(M,K))
  Sigma_EpsilonM = rep(1, M)

  if(prt){

    write.table(Ometh, file = paste0('observed_data_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Mu, file = paste0('mu_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Beta, file = paste0('beta_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Gamma, file = paste0('gamma_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(cancer.identity, file = paste0('label_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(cancer.subtype.prop, file = paste0('cancer_prop_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(cell.type.prop, file = paste0('celltype_prop_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Sigma_EpsilonM, file = paste0('epsilonM_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Sigma_MK, file = paste0('epsilonMK_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(disease_covariates, file = paste0('disease_covariates_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(z.nq, file = paste0('estep_znq_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(z.nmq, file = paste0('estep_znmq_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(c(N, M, K, Q), file = paste0('dim_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
  }

  return(list('Ometh' = Ometh, 'Covariates' = Covariates, 'Mu' = Mu, 'Beta' = Beta, 'Gamma' = Gamma, 
              'cancertype_prop' = cancer.subtype.prop, 'labels' = cancer.identity, 
              'disease_covariates' = disease_covariates, 'celltype_prop' = cell.type.prop,
              'Eps_M' = Sigma_EpsilonM, 'Eps_MK' = Sigma_MK, 'Ez' = z.nmq, 'EZ' = z.nq,
              'Initseed' = seed))
}

INIT = function(DGP, K, Q, P, species, prt = 0, ver, seed = 1){
  #initialize P_t
  Initialize = function(Ometh, num_celltype, seed){
    CorDescent = function(MethMatr, num_celltype, tol = 0.01, seed, showIter = FALSE){
      err0 = 0
      err1 = 1000
      m = nrow(MethMatr) #number of CpG sites
      n = ncol(MethMatr) #number of samples
      K = num_celltype
      if(m < n){
        stop("The CpG site number must be larger than the sample number!")
      }
      #initialize P_matr
      P_matr_t = vapply(seq_len(n), function(i){ rDirichlet(rep(2,K), seed + i) }, FUN.VALUE=rep(-1,K))
      while(abs(err1 - err0) >= tol){
        err0 = err1
        #update U_matr
        Dmat = 2*P_matr_t%*%t(P_matr_t)
        Amat = cbind(diag(rep(1,K)), diag(rep(-1,K)))
        bvec = c(rep(0,K), rep(-1,K))
        U_matr_t = t( vapply(seq_len(m), function(j){
          dvec = 2*P_matr_t %*% as.numeric(MethMatr[j,])
          
          solu = solve.QP(Dmat, dvec, Amat, bvec)
          solu$solution
        }, FUN.VALUE=rep(-1,K)) )
        
        #update P_matr
        Dmat = 2*t(U_matr_t) %*% U_matr_t
        Amat = cbind(matrix(1, K, K), diag(rep(1,K)))
        bvec = c(rep(1, K), rep(0, K))
        P_matr_t = vapply(seq_len(n), function(i){
          dvec = 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
          solu = solve.QP(Dmat, dvec, Amat, bvec, meq = K)
          solu$solution
        }, FUN.VALUE=rep(-1,K))
        
        #calculate errors
        err1 = sum((MethMatr - U_matr_t %*% P_matr_t)^2)
        if(showIter == TRUE){
          message("  ", err1, "\n")
        }
      }
      
      return(list(U=U_matr_t, P=P_matr_t))
    }
    
    K = num_celltype
    sdrow = apply(Ometh, 1, sd)
    ind = order(sdrow, decreasing = TRUE)
    m = nrow(Ometh)
    if(m <= 1000){
      num_cpg_for_init = m
    }else{
      num_cpg_for_init = max(3*ncol(Ometh), floor(m/10))
      if(num_cpg_for_init > m){
        num_cpg_for_init = m
      }
      num_cpg_for_init=m
    }
    
    #select CpG sites with the most num_cpg_for_init variant methylation levels
    Ometh_part = Ometh[ind[seq_len(num_cpg_for_init)],]
    
    result = CorDescent(Ometh_part, num_celltype=K, tol = 0.1, seed, showIter = FALSE)
    P_initial = result$P
    
    mu_initial = vapply(seq_len(m), function(j){
      if(K > 2){
        fit = lm(Ometh[j,]~as.matrix(t(P_initial[-1, ])))
      }else{
        fit = lm(Ometh[j,]~as.numeric(P_initial[-1, ]))
      }
      tmp = as.numeric(summary(fit)$coeff[ ,1])
      tmp[-1] = tmp[1] + tmp[-1]
      tmp
    }, FUN.VALUE = rep(-1,K) )
    return(list(P_initial, mu_initial))
  }
  
  N = dim(DGP$Ometh)[2]
  M = dim(DGP$Ometh)[1]
  # P = dim(DGP$Beta)[3]

  P_t = NULL
  while(is.null(P_t)){
    init = try( Initialize(DGP$Ometh, K, seed) )
    if(class(init)[1] != "try-error"){
      P_t = init[[1]]
      Mu_t = t(init[[2]])
    }
    seed = seed + 1
  }

  # print(cor(t(P_t), t(DGP$celltype_prop)))
  # print(cor(Mu_t, DGP$Mu))

  Beta_t = array(0, dim = c(M, K, P))
  Gamma_t = array(0, dim=c(M, K, Q))
  
  hat.CIden = as.vector(kmeans(t(DGP$Ometh), Q)$cluster)
  print(table(hat.CIden, DGP$labels))
  ind = rep(NA, Q)
  for(q in 1:Q){
    ind[q] = as.numeric(readline(paste0("Enter the ", q, " th disease subtype's label: ")))
  }
  
  # ind = as.numeric(names(sort(table(hat.CIden), decreasing =  T)))
  for(n in 1:length(hat.CIden)){
    hat.CIden[n] = which(hat.CIden[n] == ind)
  }
  print(table(hat.CIden, DGP$labels))
  
  Pi_t = as.vector(table(hat.CIden)) / sum(as.vector(table(hat.CIden)))
  # Pi_t = sort(Pi_t, decreasing = T)

  Ez = array(0, dim = c(N, M, Q))
  EZ = array(0, dim = c(N, Q))
  for(n in 1:N){
    EZ[n, hat.CIden[n]] = 1
    Ez[n, , hat.CIden[n]] = 1
  }
  
  # m = M / 3
  # for(q in 2:Q){
  #   Ometh_q = DGP$Ometh[,hat.CIden == q]
  #   Gamma_t[,,q] = Mu_t - apply(Ometh_q, 1, mean)
  #   for(k in 1:K){
  #     ind = order(abs(Gamma_t[,k,q]), decreasing = T)[1:m]
  #     Gamma_t[-ind,k,q] = 0
  #   }
  # }

  # ind = which.max(apply(EZ[1:(N*Pi_t[1]),], 2, sum))
  # temp = EZ[,ind]
  # EZ[,ind] = EZ[,1]
  # EZ[,1] = temp
  # print(apply(EZ, 2, sum))

  if(prt){
    write.table(Mu_t, file = paste0('mu_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Beta_t, file = paste0('beta_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Gamma_t, file = paste0('gamma_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Pi_t, file = paste0('cancer_prop_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(Ez, file = paste0('estep_znmq_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(EZ, file = paste0('estep_znq_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
    write.table(c(N, M, K, Q), file = paste0('dim_', species,"_v", ver,'.txt'),
                row.names = FALSE, col.names = FALSE)
  }
  return(list('P_t' = P_t, 'Pi_t' = Pi_t, 'Mu' = Mu_t, 'Beta' = Beta_t, 
              'Gamma' = Gamma_t, 'Ez' = Ez, 'EZ' = EZ))
}

#############################
# visualization algo
#############################
library(scales)
require(gridExtra)
library(dplyr)
library(tidyr)
library(caret)
# Mu
heatmap.mu = function(DGP, ret_list, ret_list2){
  N = dim(DGP$Ometh)[2]
  M = dim(DGP$Ometh)[1]
  K = dim(ret_list$P_t)[1]
  Q = length(ret_list$Pi_t)
  
  Mu = DGP$Mu
  est_Mu = ret_list$Mu
  est_Mu2 = ret_list2$Mu
  
  cancer.type.cor = cor(Mu, est_Mu)
  row.names(cancer.type.cor) = paste0('Underlying cell type ', 1:K)
  colnames(cancer.type.cor) = paste0('Estimated cell type ', 1:K)
  print(cancer.type.cor)

  ind = 1:K
  for(k in 1:K)
    ind[k] = as.numeric(readline(prompt=paste0("Enter underlying cell type vs estimated cell type ", k, ": ")))

  Mu = Mu[,ind]
  tmp = cbind(Mu, est_Mu, est_Mu2)
  colnames(tmp) = c(paste0('cell type ', ind), 
                    paste0('cell type ', 1:K, ' '), 
                    paste0('cell type ', 1:K, '  '))
  # colnames(tmp) = c(ind, paste0(1:K, 'M'), paste0(1:K, 'H'))
  tmp = melt(tmp)
  # tmp$Var2 = as.character(tmp$Var2)
  # for(k in 1:K){
  #   tmp$Var2 = ifelse(tmp$Var2 == paste0(k, 'M'), paste0(' ', k), 
  #                     ifelse(tmp$Var2 == paste0(k, 'H'), paste0('  ', k), tmp$Var2))
  # }
  # tmp$Var2 = as.factor(tmp$Var2)
  colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
  class.names = c('truth', 'clusterHIRE', 'HIRE')
  
  for(i in 1:3){
    tmp2 = tmp[(1 + M*K*(i-1)):(M*K*i),]
    
    png(paste0('N', N, 'K', K, 'Q', Q, '_', class.names[i], '_mu.png'),
        width = 1920, height = 1080)
    p = ggplot(tmp2, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue",
                           midpoint = 0,
                           mid = "#fcfcfc",
                           high = "red",
                           space="Lab", limits = range(tmp$Val)) +
      theme(axis.text.x = element_text(angle = 0, face = "bold", size = 36),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.title.x = element_blank()) +
      theme(axis.title = element_text(size=32,face="bold"),
            legend.text = element_text(size=20, face="bold"),
            legend.title = element_text(size = 20,face="bold")) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0))
    print(p)
    dev.off()
  }
}

# Beta
heatmap.beta = function(DGP, ret_list, ret_list2){
  N = dim(DGP$Ometh)[2]
  K = dim(ret_list$P_t)[1]
  P = dim(ret_list$Beta)[3]
  
  Beta = DGP$Beta[,,1]
  est_Beta = ret_list$Beta[,,1]
  est_Beta2 = ret_list2$Beta[,,1]
  for(p in 1:P){
    Beta = cbind(Beta, DGP$Beta[,,p])
    est_Beta = cbind(est_Beta, ret_list$Beta[,,p])
    est_Beta2 = cbind(est_Beta2, ret_list2$Beta[,,p])
  }
  
  for(p in 1:P){
    p3 = K*p + 1:K
    for(p1 in 1:P){
      p2 = K*p1 + 1:K
      phenotype.cor = cor(Beta[,p3], est_Beta[,p2])
      row.names(phenotype.cor) = paste0('Underlying cell type ', 1:K)
      colnames(phenotype.cor) = paste0('Estimated cell type ', 1:K)
      message(paste0('The underlying phenotype ', p, ' vs estimated phenotype ', p1, '  \n'))
      print(phenotype.cor)
    }
  }
  
  ind = rep(NA, K)
  for(p in 1:P){
    p3 = K*p + 1:K
    p1 = as.numeric(readline(prompt=paste0("Enter underlying phenotype vs estimated phenotype ", p, ": ")))
    p2 = K*p1 + 1:K
    phenotype.cor = cor(Beta[,p3], est_Beta[,p2])
    row.names(phenotype.cor) = paste0('Underlying cell type ', 1:K)
    colnames(phenotype.cor) = paste0('Estimated cell type ', 1:K)
    print(phenotype.cor)
    
    if(is.na(ind[1])){
      for(k in 1:K)
        ind[k] = 
          as.numeric(readline(prompt=paste0("Enter underlying cell type vs estimated cell type ", k, ": ")))
    }
    
    p3 = K*p + ind
    tmp = cbind(Beta[,p3], est_Beta[,p2], est_Beta2[,p2])
    colnames(tmp) = c(paste0('cell type ', ind), 
                      paste0('cell type ', 1:K, ' '), 
                      paste0('cell type ', 1:K, '  '))
    # colnames(tmp) = c(ind, paste0(1:K, 'M'), paste0(1:K, 'H'))
    tmp = melt(tmp)
    # tmp$Var2 = as.character(tmp$Var2)
    # for(k in 1:K){
    #   tmp$Var2 = ifelse(tmp$Var2 == paste0(k, 'M'), paste0(' ', k), 
    #                     ifelse(tmp$Var2 == paste0(k, 'H'), paste0('  ', k), tmp$Var2))
    # }
    # tmp$Var2 = as.factor(tmp$Var2)
    colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
    class.names = c('truth', 'clusterHIRE', 'HIRE')
    
    for(i in 1:3){
      tmp2 = tmp[(1 + M*K*(i-1)):(M*K*i),]
      
      png(paste0('N', N, 'K', K, 'Q', Q, '_', class.names[i], '_', p, '_pheno.png'),
          width = 1920, height = 1080)
      pp = ggplot(tmp2, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue",
                             midpoint = 0,
                             mid = "#fcfcfc",
                             high = "red",
                             space="Lab", limits = range(tmp$Val)) +
        theme(axis.text.x = element_text(angle = 0, face = "bold", size = 36),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), axis.title.x = element_blank()) +
        theme(axis.title = element_text(size=32,face="bold"),
              legend.text = element_text(size=20, face="bold"),
              legend.title = element_text(size = 20,face="bold")) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
      print(pp)
      dev.off()
    }
  }
}

# Gamma
heatmap.gamma = function(DGP, ret_list, ret_list2){
  N = dim(DGP$Ometh)[2]
  K = dim(ret_list$P_t)[1]
  Q = length(ret_list$Pi_t)
  
  Gamma = DGP$Gamma[,,1]
  est_Gamma = ret_list$Gamma[,,1]
  est_Gamma2 = ret_list2$Gamma[,,1]
  for(q in 2:Q){
    Gamma = cbind(Gamma, DGP$Gamma[,,q])
    est_Gamma = cbind(est_Gamma, ret_list$Gamma[,,q])
    est_Gamma2 = cbind(est_Gamma2, ret_list2$Gamma[,,q])
  }
  
  for(q in 2:Q){
    q3 = K*(q - 1) + 1:K
    for(q1 in 2:Q){
      q2 = K*(q1 - 1) + 1:K
      cancer.type.cor = cor(Gamma[,q3], est_Gamma[,q2])
      row.names(cancer.type.cor) = paste0('Underlying cell type ', 1:K)
      colnames(cancer.type.cor) = paste0('Estimated cell type ', 1:K)
      message(paste0('The underlying cancer subtype ', q, ' vs estimated cancer subtype ', q1, '  \n'))
      print(cancer.type.cor)
    }
  }

  ind = rep(NA, K)
  for(q in 2:Q){
    q3 = K*(q - 1) + 1:K
    q1 = as.numeric(readline(prompt=paste0("Enter underlying cancer type vs estimated cancer type ", q, ": ")))
    q2 = K*(q1 - 1) + 1:K
    cancer.type.cor = cor(Gamma[,q3], est_Gamma[,q2])
    row.names(cancer.type.cor) = paste0('Underlying cell type ', 1:K)
    colnames(cancer.type.cor) = paste0('Estimated cell type ', 1:K)
    print(cancer.type.cor)

    if(is.na(ind[1])){
      for(k in 1:K)
        ind[k] = 
          as.numeric(readline(prompt=paste0("Enter underlying cell type vs estimated cell type ", k, ": ")))
    }
    
    q3 = K*(q - 1) + ind
    tmp = cbind(Gamma[,q3], est_Gamma[,q2], est_Gamma2[,q2])
    colnames(tmp) = c(paste0('cell type ', ind), 
                      paste0('cell type ', 1:K, ' '), 
                      paste0('cell type ', 1:K, '  '))
    # colnames(tmp) = c(ind, paste0(1:K, 'M'), paste0(1:K, 'H'))
    tmp = melt(tmp)
    # tmp$Var2 = as.character(tmp$Var2)
    # for(k in 1:K){
    #   tmp$Var2 = ifelse(tmp$Var2 == paste0(k, 'M'), paste0(' ', k), 
    #                     ifelse(tmp$Var2 == paste0(k, 'H'), paste0('  ', k), tmp$Var2))
    # }
    # tmp$Var2 = as.factor(tmp$Var2)
    colnames(tmp) = c('CpG_sites', 'Cell_types', 'Val')
    class.names = c('truth', 'clusterHIRE', 'HIRE')
    
    for(i in 1:3){
      tmp2 = tmp[(1 + M*K*(i-1)):(M*K*i),]
      
      png(paste0('N', N, 'K', K, 'Q', Q, '_', class.names[i], '_', q, '.png'),
          width = 1920, height = 1080)
      p = ggplot(tmp2, aes(x = Cell_types, y = CpG_sites, fill = Val)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue",
                             midpoint = 0,
                             mid = "#fcfcfc",
                             high = "red",
                             space="Lab", limits = range(tmp$Val)) +
        theme(axis.text.x = element_text(angle = 0, face = "bold", size = 36),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), axis.title.x = element_blank()) +
        theme(axis.title = element_text(size=32,face="bold"),
              legend.text = element_text(size=20, face="bold"),
              legend.title = element_text(size = 20,face="bold")) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
      print(p)
      dev.off()
    }
  }
}

# estimated cell compositions vs the truth
scatter.plot = function(DGP, ret_list){
  N = dim(DGP$Ometh)[2]
  K = dim(ret_list$P_t)[1]
  Q = length(ret_list$Pi_t)
  
  true_P = DGP$celltype_prop
  est_P = ret_list$P_t
  
  cell.type.cor = cor(t(true_P), t(est_P))
  row.names(cell.type.cor) = paste0('Underlying cell type ', 1:K)
  colnames(cell.type.cor) = paste0('Estimated cell type ', 1:K)
  print(cell.type.cor)

  png(paste0('N', N, 'K', K, 'Q', Q, '_celltype_prop.png'),
      width = 1920, height = 1080)
  par(mar = c(15,15,5,2) + 0.1) ## default is c(5,4,4,2) + 0.1
  par(mfrow=c(1,K))
  for (k in 1:K) {
    k1 = as.numeric(readline(prompt=paste0("Enter underlying cell type on the y-axis vs estimated cell type ", k, ": ")))
    plot(est_P[k, ], true_P[k1, ], xlim=c(0,1), ylim=c(0,1),
         main='', xlab=paste0('Estimated cell type ', k),
         ylab=paste0('Underlying cell type ', k1),
         cex.main = 8, cex.lab = 5) + abline(a=0, b=1, col="red")
  }
  dev.off()
}

# estimated cancer compositions vs the truth
pie.cancer.prop = function(DGP, ret_list){
  N = dim(DGP$Ometh)[2]
  K = dim(ret_list$P_t)[1]
  Q = length(ret_list$Pi_t)
  
  true_Pi = DGP$cancertype_prop
  est_Pi = ret_list$Pi_t
  
  cancer.prop = data.frame('disease_subtype' = factor(1:Q),
                           'prop' = sort(true_Pi, decreasing = T))
  # cancer.prop$lab.ypos = cumsum(cancer.prop$prop) - pos*cancer.prop$prop

  estimated.prop = data.frame('disease_subtype' = factor(1:Q),
                              'prop' = round(sort(est_Pi, decreasing = T), 3))
  # estimated.prop$lab.ypos = cumsum(estimated.prop$prop) - pos*estimated.prop$prop

  plot1 = ggplot(cancer.prop, aes(x = "", y = prop, fill = disease_subtype)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(label = prop), position = position_stack(vjust = 0.5),
              size = 20, color = "white") +
    scale_fill_manual(values = hue_pal()(Q)) +
    theme_void() +
    labs(title="Underlying proportion") +
    theme(legend.key.size = unit(3, 'cm'), legend.position="bottom",
          legend.title = element_text(size=40), legend.text = element_text(size=34),
          plot.title = element_text(hjust=0.5, size=32))

  plot2 = ggplot(estimated.prop, aes(x = "", y = prop, fill = disease_subtype)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(label = prop), position = position_stack(vjust = 0.5),
              size = 20, color = "white") +
    scale_fill_manual(values = hue_pal()(Q)) +
    theme_void() +
    labs(title="Estimated proportion") +
    theme(legend.key.size = unit(3, 'cm'), legend.position="bottom",
          legend.title = element_text(size=40), legend.text = element_text(size=34),
          plot.title = element_text(hjust=0.5, size=32))

  png(paste0('N', N, 'K', K, 'Q', Q, '_subtype_prop.png'),
      width = 1920, height = 1080)
  grid.arrange(plot1, plot2, ncol=2)

  dev.off()
}

# confusion matrix of indiviual labels
cm.cancer = function(DGP, ret_list, strings){
  N = dim(DGP$Ometh)[2]
  K = dim(ret_list$P_t)[1]
  Q = length(ret_list$Pi_t)
  
  true_C = DGP$labels
  est_C = apply(ret_list$E_Z, 1, which.max)
  
  est.num = est_C + 1000
  for(q in 1:Q){
    ind = as.numeric(names(sort(table(est.num), decreasing = T))[q])
    est.num[which(est.num == ind)] = q
  }

  result = data.frame(confusionMatrix(data = factor(est.num),
                                      reference = factor(true_C))$table)
  result$Reference = factor(result$Reference, levels = Q:1)
  result = result %>% group_by(Reference) %>% mutate(prop = Freq/sum(Freq))
  if(strings == 'gender'){
    result$Prediction = ifelse(result$Prediction == 1, 'Female', 'Male')
    result$Reference = ifelse(result$Reference == 1, 'Female', 'Male')
  }
  png(paste0('N', N, 'K', K, 'Q', Q, '_subtype_heatmap.png'),
      width = 1920, height = 1080)
  p = ggplot(result, aes(x=Prediction, y=Reference, fill=prop)) + geom_tile()  +
    labs(title=paste0('Heatmap of underlying ',strings, ' vs estimated ',strings, ' among individuals')) +
    theme(text = element_text(size=28), legend.key.size = unit(2, 'cm'), legend.position="right",
          legend.title = element_text(size=28), legend.text = element_text(size=20),
          plot.title = element_text(hjust=0.5, size=32)) +
    scale_fill_gradient(low = "white", high = "red") +
    guides(fill=guide_legend(title="Proportion"))
  print(p)
  dev.off()
}

#############################
# HIRE and MPCS Algo
#############################
HIRE = function(Ometh, X, mu_t, gamma_t, P_t, iter=1000, core=32, reshape=0){

  HIRERcallC = function(Ometh, X, P_t, mu_t, gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core){
    args = list("P_init"=as.numeric(P_t), "Mu_init"=as.numeric(mu_t),
                "Gamma_init"=as.numeric(gamma_t), "Gamma_init_dim"=as.integer(dim(gamma_t)),
                "Ometh_r"=as.numeric(Ometh), "Ometh_r_dim"=as.integer(dim(Ometh)),
                "X_r"=as.numeric(X), "Eps_MK_init"=as.numeric(Eps_MK_t),
                "Eps_M_init"=as.numeric(Eps_M_t), "num_iter" = as.integer(iter),
                "num_core" = as.integer(core))
    ret_list = .Call("HIRERcallC", args)
    return(ret_list)
  }

  M = nrow(Ometh) #CpG site number
  N = ncol(Ometh) #sample number
  Q = dim(gamma_t)[3]
  K = dim(gamma_t)[2]

  Eps_MK_t = matrix(1, M, K)
  Eps_M_t = rep(1, M)

  if(reshape){
    Gamma_t = array(0, dim=c(M, K, Q))
    for(q in 1:Q){
      Gamma_t[,,q] = gamma_t[,(1 + K*(q-1)):(K*q)]
    }
  } else{
    Gamma_t = gamma_t
  }

  message("  Implementing EM algorithm... \n")
  ret_list = HIRERcallC(Ometh, X, P_t, mu_t, Gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core)
  message("  Done! \n")

  return(ret_list)
}

MPCS = function(Ometh, EZ, Pi_t, mu_t, gamma_t, P_t, iter=1000, core=32, reshape=0){

  MPCSRcallC = function(Ometh, EZ, Pi_t, P_t, mu_t, gamma_t, Eps_MK_t, Eps_M_t, iter, core){
    args = list("Pi_init"=as.numeric(Pi_t), "P_init"=as.numeric(P_t), "Mu_init"=as.numeric(mu_t),
                "Gamma_init"=as.numeric(gamma_t), "Gamma_init_dim"=as.integer(dim(gamma_t)),
                "Ometh_r"=as.numeric(Ometh), "Ometh_r_dim"=as.integer(dim(Ometh)),
                "Z_r"=as.numeric(EZ), "Eps_MK_init"=as.numeric(Eps_MK_t),
                "Eps_M_init"=as.numeric(Eps_M_t), "num_iter" = as.integer(iter),
                'num_core' = as.integer(core))
    ret_list = .Call("MPCSRcallC", args)
    return(ret_list)
  }

  M = nrow(Ometh) #CpG site number
  N = ncol(Ometh) #sample number
  K = dim(gamma_t)[2]
  Q = dim(gamma_t)[3]

  Eps_MK_t = matrix(1, M, K)
  Eps_M_t = rep(1, M)

  if(reshape){
    Gamma_t = array(0, dim=c(M, K, Q))
    for(q in 1:Q)
      Gamma_t[,,q] = gamma_t[,(1 + K*(q-1)):(K*q)]
  } else{
    Gamma_t = gamma_t
  }

  # message("  Initialization Done.\n")

  message("  Implementing EM algorithm... \n")
  ret_list = MPCSRcallC(Ometh, EZ, Pi_t, P_t, mu_t, Gamma_t, Eps_MK_t, Eps_M_t, iter, core)
  message("  Done! \n")

  return(ret_list)
}

#############################
# Parameters
#############################
N = 2000
M = 10000
K = 4
Q = 3
species = 'Humans'
ver = 4

# # change dir
# setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
# setwd(paste0('./Q', Q))

# change dir
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE")
subDir = paste0('N', N, 'K', K, 'Q', Q, 'V', ver)

if (file.exists(subDir)){
  setwd(subDir)
} else {
  dir.create(subDir)
  setwd(subDir)
}

#############################
# generate data
#############################
myCluster = makeCluster(detectCores(), type = "PSOCK")
registerDoParallel(myCluster)
DGP = MPCS.DGP(N, M, K, Q, species, 0, ver, 123)
stopCluster(myCluster)

Init = INIT(DGP, K, Q, dim(DGP$Beta)[3], species, 0, ver, DGP$Initseed)
rm(list=setdiff(ls(), c('DGP', 'Init')))
save.image("data.RData")

load('data.RData')
print(cor(t(Init$P_t), t(DGP$celltype_prop)))
print(cor(Init$Mu, DGP$Mu))
print(cor(Init$Gamma[,,2], DGP$Gamma[,,2]))

#############################
# debug the MPCS algo
#############################
# HIRE algo
# setwd('C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/HIREewas/src')
# system("R CMD SHLIB HIRE.c")
# dyn.load('HIRE.dll')
#
# time0 = proc.time()
# ret_list = HIRE(DGP$Ometh, DGP$disease_covariates, DGP$Mu, DGP$Gamma, DGP$celltype_prop, DGP$Eps_MK,
#                 DGP$Eps_M, K, num_iter=100)
# print(proc.time() - time0)

vers = 1
setwd('C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/HIREewas/src')
system(paste0('R CMD SHLIB HIRE_v', vers, '.c'))
dyn.load(paste0('HIRE_v', vers, '.dll'))

library(HIREcewas)
time0 = proc.time()
ret_list = HIREcewas::HIRE(DGP$Ometh, DGP$disease_covariates, Init$Mu, Init$Gamma, Init$P_t,
                            iter=1000, tol = 1e-10, core = 32)
print(proc.time() - time0)

setwd('C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/HIREewas/src')
system(paste0('R CMD SHLIB MPCS_v', 6, '.c'))
dyn.load(paste0('MPCS_v', 6, '.dll'))

library(cHIRE4)
time0 = proc.time()
# ret_list=cHIRE4::cHIRE(DGP$Ometh, DGP$Covariates, Init$EZ, Init$Pi_t, Init$Mu, Init$Beta,
#                        Init$Gamma, Init$P_t, iter=1000, tol = 1e-20, core=20)
ret_list=cHIRE4::cHIRE(DGP$Ometh, DGP$Covariates, t(DGP$disease_covariates), DGP$cancertype_prop, 
                       Init$Mu, Init$Beta, Init$Gamma, DGP$celltype_prop, iter=1000, tol = 1e-20, core=20)
# ret_list=cHIRE4::cHIRE(DGP$Ometh, DGP$Covariates, t(DGP$disease_covariates), DGP$cancertype_prop, 
#                        DGP$Mu, DGP$Beta, DGP$Gamma, DGP$celltype_prop, iter=100, tol = 1e-20, core=20)
rm(Init)
print(proc.time() - time0)
# save.image("cHIRE4.RData")
save.image("cHIRE4_.RData")

# change dir
# rm(list=setdiff(ls(), c('DGP', 'Init', 'ret_list', 'ret2_list', 'ret3_list')))
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE")
setwd(paste0('N', N, 'K', K, 'Q', Q, 'V', ver))


#############################
# pBIC
#############################
K = 2:7; Q = 2:5
BIC = array(NA, dim = c(length(K), length(Q)), dimnames = list(K, Q))
# BIC = setNames(rep(NA, length(Q)), Q)

BIC_cal = function(DGP, ret_list){
  N = dim(DGP$Ometh)[2]
  M = dim(DGP$Ometh)[1]
  K = dim(ret_list$Gamma)[2]
  Q = dim(ret_list$Gamma)[3]
  P = dim(ret_list$Beta)[3]
  ploglike = 0
  
  for(m in 1:M) {
    for(n in 1:N) {
      s1 = 0
      s3 = 0
      s4 = 0
      sum = 0
      for (k in 1:K) {
        s1 = s1 + ret_list$P_t[k,n] * ret_list$Mu[m,k]
        s3 = s3 + ret_list$P_t[k,n] * ret_list$P_t[k,n] * ret_list$Eps_MK[m,k]
        for(p in 1:P){
          s4 = s4 + ret_list$P_t[k,n] * ret_list$Beta[m,k,p]
        }
      }
      
      for (q in 1:Q) {
        s2 = 0
        for (k in 1:K) {
          s2 = s2 + ret_list$P_t[k,n] * ret_list$Gamma[m,k,q]
        }
        sum = sum + ret_list$Pi_t[q] * 
          dnorm(DGP$Ometh[m,n], s1 + s2 + s4, sqrt(s3 + ret_list$Eps_M[m]))
      }
      ploglike = ploglike + log(sum)
    }
  }
  return (-2.0 * ploglike + log(N) * (Q - 1 + (K - 1) * N + M * (1 + 2 * K + (Q-1) * K + P*K)))
}

pBIC_cal = function(DGP, ret_list, alpha = .01){
  N = dim(DGP$Ometh)[2]
  M = dim(DGP$Ometh)[1]
  K = dim(ret_list$Gamma)[2]
  Q = dim(ret_list$Gamma)[3]

  # # message("  Calculating p-values...\n")
  # tmp = NULL
  # for(q in 2:Q){
  #   tmp = cbind(tmp, ret_list$E_Z[,q]*t(ret_list$P_t))
  # }
  # # for(q in 2:Q){
  # #   tmp = cbind( tmp, ret_list$E_Z[,q]*t(ret_list$P_t)[, seq(2,K)])
  # # }
  # # x_matr = cbind( tmp, t(ret_list$P_t)[, seq(2,K)])
  # x_matr = as.matrix(tmp)
  # 
  # pvalues = t(vapply(seq_len(M), function(m){
  #   fit = lm(DGP$Ometh[m,]~x_matr)
  #   summary(fit)$coef[seq(2, (1+(Q-1)*K)),4]
  # }, FUN.VALUE = rep(-1,(Q-1)*K)))
  # # message("  Done!\n")
  # 
  # ret_list[[9]] = pvalues
  # names(ret_list)[9] = "pvalues"
  # names(ret_list)[8] = "pBIC"
  
  d = (Q-1) + (K-1)*N + M*(1+2*K+Q*K) #number of parameters
  d1 = (Q-1) + (K-1)*N + M*(1+2*K+(Q-1)*K) #number of parameters
  # d0 = sum(ret_list$pvalues > alpha/((Q-1)*M*K)) #number of zero parameters
  ret_list[[8]] = ret_list[[8]] - log(N)*d + log(N)*d1
  # ret_list[[8]] = ret_list[[8]] - log(N)*d + log(N)*(d1-d0)
  
  return(ret_list[[8]])
}

AICcal = function(DGP, ret_list, alpha = .01){
  N = dim(DGP$Ometh)[2]
  M = dim(DGP$Ometh)[1]
  K = dim(ret_list$Gamma)[2]
  Q = dim(ret_list$Gamma)[3]
  
  num.para = (Q-1) + (K-1)*N + M*(1+2*K+(Q-1)*K) #number of parameters
  loglike = (ret_list[[8]] - log(N)*num.para) * -.5
  
  # aic = -2 / N * loglike + 2 * num.para / N
  aic = 2*num.para - 2 * loglike
  return(aic)
}

for(jj in Q){ # Q
  for(ii in K){ # K
    message(ii, ' ', jj)
    # change dir
    setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/GSE42861V2/")
    setwd(paste0('./Q', jj, './K', ii))
    # setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE")
    # setwd(paste0('N1000', 'K', ii, 'Q', jj, 'V0'))
    load("MPCSewas4.RData")
    # rm(list=setdiff(ls(), c('DGP', 'Init', 'ret_list', 'ii', 'jj')))
    # save.image("MPCSewas4.RData")
    # print(ret_list$BIC)
    # print(AICcal(DGP, ret_list, .01))
    BIC[which(ii == K), which(jj == Q)] = AICcal(DGP, ret_list, .01)
  }
}

# options(scipen = -1)
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/GSE42861V2/")
setwd("./Q2/K6/")
# setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE")
# setwd("./N1000K4Q3V0/")

## reshape data (tidy/tall form)
library(tidyverse)
dat2 = BIC %>% tbl_df() %>% rownames_to_column('K') %>% gather(Q, value, -K) %>%
  mutate(
    K = factor(K, levels=1:nrow(BIC)),
    Q = factor(gsub("V", "", Q), levels=1+1:ncol(BIC))
  )
dat2$K = as.factor(as.numeric(dat2$K) + 1)
## plot data
png("AIC.png", width = 1920, height = 1080)
ggplot(dat2, aes(K, Q)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1)),  size=12) +
  scale_fill_gradient(low = "red", high = "white") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 32),
        axis.text.y = element_text(face = "bold", size = 32),
        axis.ticks.y = element_blank()) +
  theme(axis.title = element_text(size=32,face="bold"),
        legend.position="none")
dev.off()

# par(mar = c(7,10,7,10) + 0.1) ## default is c(5,4,4,2) + 0.1
# plot(names(BIC[3,]), BIC[3,], type="n", xlab = 'Number of cancer subtypes', ylab = "AIC",
#      main = 'AICs among different number of cancer subtypes with K = 4', cex.lab = 4,
#      las = 0, cex.main = 3, yaxt = "n", xaxt = "n")
# axis(1, at = 1:7, cex.axis = 5)
# axis(2, cex.axis = 1.5)
# points(names(BIC[3,]), BIC[3,], type="b", pch=19, cex=2)
# abline(v = names(which.min(BIC[3,])), h = min(BIC[3,]), col = 'red')
# # mtext(paste0('Q = ', names(which.min(BIC[3,]))), side=1, line=-2, at=5, col = 'red', cex = 5)
# dev.off()

#############################
# visualization
#############################
load('HIREcewasV.RData')
load("cHIRE4_.RData")

# comparison
cHIRE.data = ret_list
hire.data = ret_list

scatter.plot(DGP, cHIRE.data)
scatter.plot(DGP, hire.data)
heatmap.mu(DGP, cHIRE.data, hire.data)
heatmap.beta(DGP, cHIRE.data, hire.data)
heatmap.gamma(DGP, cHIRE.data, hire.data)
pie.cancer.prop(DGP, cHIRE.data)
DGP$labels = DGP$disease_covariates[3,] + 1
cm.cancer(DGP, cHIRE.data, 'disease')
table(ifelse(cHIRE.data$E_Z > .5, 1, 0)[,1], DGP$disease_covariates[1,])
dev.off()

graph.plot = function(num_samples){
  for(n in num_samples){
    for(k in 3:6){
      for(q in 2:5){
        if(k >= q){

          setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE")
          setwd(paste0('N', n, 'K', k, 'Q', q, 'V', 0))
          load('data.RData')
          load("MPCSewas0.RData")
          data = ret_list
          true_Gamma = DGP$Gamma[,,1]
          est_Gamma = data$Gamma[,,1]
          for(q1 in 2:q){
            true_Gamma = cbind(true_Gamma, DGP$Gamma[,,q1])
            est_Gamma = cbind(est_Gamma, data$Gamma[,,q1])
          }

          scatter.plot(DGP$celltype_prop, data$P_t, n, k, q)
          heatmap.mu(DGP$Mu, data$Mu, n, k, q)
          heatmap.gamma(true_Gamma, est_Gamma, n, k, q)
          pie.cancer.prop(DGP$cancertype_prop, data$Pi_t, n, k, q)
          cm.cancer(DGP$labels, apply(data$E_Z, 1, which.max), n, k, q)
          cat('\014')
        }
      }
    }
  }
}
graph.plot(c(1500))

#############################
# Real data analysis Luo
#############################
library(sva)
library(meffil)
library(data.table)
library(rvest)
library(RSelenium)
# GSE42861
# change dir
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
load('data.RData')
data = fread('GSE42861_processed_methylation_matrix.txt', data.table = F)
CpG.names = unlist(read.table('./Luo/RA/cpg_names_after_quality_control_10K.txt'))
# data = read.table('processed_methylation_matrix_2.txt', header = T, row.names = 1)
row.names(data) = data[,1]
data = data[,-1]
data = data[,-dim(data)[2]]
data = data[,-c(11, 167)] # remove samples GSM1051535 and GSM1051691
data = data[-which(apply(data, 1, mean) < .2 | apply(data, 1, mean) > .8),]
print(data[1:6,1:6])

batch.names = unique(gsub("_.*","",colnames(data)))
batches = rep(NA, ncol(data))
for(n in 1:ncol(data)){
  batches[n] = grep(gsub("_.*","",colnames(data)[n]), batch.names)
}

disease_covariates = DGP$disease_covariates[-1,]
disease_covariates = rbind(1, disease_covariates)
data = ComBat(dat = data, batch = batches, mod = t(disease_covariates))
print(data[1:6,1:6])

mVG = meffil.most.variable.cpgs(as.matrix(data), n = 10000, autosomal = F)
print(sum(mVG %in% CpG.names))
data = data[row.names(data) %in% mVG,]

# disease_covariates
URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
URLs = paste0(URL, read.table('GSE42861.txt')$V3)

setwd('C:/Program Files/Google/Chrome Gamma/Application')
# java -Dwebdriver.chrome.driver=chromedriver.exe -jar selenium-server-standalone-3.12.0.jar -port 2444

setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
remDr = remoteDriver(remoteServerAddr = "localhost", port = 2430, browserName = "chrome")
remDr$open()
Sys.sleep(3)
webElems = rep(NA, dim(DGP$Ometh)[2])

for(n in 1:dim(DGP$Ometh)[2]){
  remDr$navigate(URLs[n])
  Sys.sleep(5)
  temp = remDr$getPageSource()[[1]] %>% read_html() %>% html_table()
  webElems[n] = temp[[16]]$X2[9]
  Sys.sleep(3)
}

webElems = strsplit(webElems, " ")

# remove samples 11 and 167
DGP$Ometh = DGP$Ometh[,-c(11, 167)]
webElems = webElems[-c(11, 167)]

disease_covariates = array(0, dim = c(6, dim(DGP$Ometh)[2]))
for(n in 1:dim(DGP$Ometh)[2]){
  if(webElems[[n]][1] == 'Patient')
    disease_covariates[1,n] = 1

  disease_covariates[2,n] = webElems[[n]][2]

  if(webElems[[n]][3] == 'm')
    disease_covariates[3,n] = 1

  if(webElems[[n]][4] == 'ex')
    disease_covariates[4,n] = 1
  else if(webElems[[n]][4] == 'current')
    disease_covariates[5,n] = 1
  else if(webElems[[n]][4] == 'occasional')
    disease_covariates[6,n] = 1
}

disease_covariates = matrix(as.numeric(disease_covariates), nrow = 6)

DGP = list()
DGP$Ometh = as.matrix(data)
DGP$disease_covariates_cluster = as.matrix(disease_covariates)[-1,]
rm(list=setdiff(ls(), c('DGP')))
save.image("data.RData")

setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/clusterHIRE/GSE42861")
load('data.RData')
DGP$covariates = DGP$covariates[-1,] # remove the labels RA
# min max normalize on Age
DGP$covariates[1,] = 
  (DGP$covariates[1,] - min(DGP$covariates[1,])) / (max(DGP$covariates[1,]) - min(DGP$covariates[1,]))

for(q in 2:5){
  subDir = paste0('Q', q)
  if (file.exists(subDir)){
    setwd(subDir)
  } else {
    dir.create(subDir)
    setwd(subDir)
  }

  for(k in 2:9){
    subDir = paste0('K', k)
    if (file.exists(subDir)){
      setwd(subDir)
    } else {
      dir.create(subDir)
      setwd(subDir)
    }

    print(dim(DGP$Covariates)[2])
    Init = INIT(DGP, k, q, dim(DGP$Covariates)[2], 'Humans', 0, 0, 123)
    rm(list=setdiff(ls(), c('DGP', 'Init', 'rDirichlet', 'INIT', 'k', 'q')))
    save.image("data.RData")
    setwd('../')
  }
  setwd('../')
}

# implementation
library(cHIRE4)
for(q in 2){
  setwd(paste0('Q', q))
  for(k in 2:9){
    setwd(paste0('K', k))
    load('data.RData')
    time0 = proc.time()
    ret_list=cHIRE4::cHIRE(DGP$Ometh, DGP$Covariates, Init$EZ, Init$Pi_t, Init$Mu, Init$Beta,
                           Init$Gamma, Init$P_t, iter=2000, tol = 1e-20, core=20)
    print(proc.time() - time0)
    save.image("cHIRE4_.RData")
    cat("\014")
    rm(list = ls())
    gc()
    setwd('../')
  }
  setwd('../')
}


### GSE77716
setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/GSE77716")

# data
data = as.data.frame(fread('GSE77716_Matrix_processed.tsv'))
row.names(data) = data[,1]
data = data[,-1]
data = data[,-seq(2, ncol(data), 2)]
data = data[,order(colnames(data))]
data = data[,-238] # remove sample GSM2057284
data = data[-which(apply(data, 1, mean) < .2 | apply(data, 1, mean) > .8),]
print(data[1:6,1:6])

mVG = meffil.most.variable.cpgs(as.matrix(data), n = 10000, autosomal = F)
data = data[row.names(data) %in% mVG,]

# disease_covariates
URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
URLs = paste0(URL, read.table('GSE77716.txt')$V3)

setwd('C:/Program Files/Google/Chrome Gamma/Application')
# java -Dwebdriver.chrome.driver=chromedriver.exe -jar selenium-server-standalone-3.12.0.jar -port 2420

setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/GSE77716")
library(RSelenium)

remDr = remoteDriver(remoteServerAddr = "localhost", port = 2420, browserName = "chrome")
remDr$open()
Sys.sleep(3)
webElems = rep(NA, length(URLs))

for(n in 1:length(URLs)){
  remDr$navigate(URLs[n])
  Sys.sleep(5)
  temp = remDr$getPageSource()[[1]] %>% read_html() %>% html_table()
  webElems[n] = temp[[16]]$X2[9]
  Sys.sleep(2)
}

webElems = gsub('sample.*ethnicity: ', ' ', webElems)
webElems = gsub('.*: ', '', webElems)
webElems = webElems[-238]
webElems = gsub('ale ', 'ale:', webElems)
webElems = strsplit(webElems, ":")

disease_covariates = array(0, dim = c(4, length(webElems)))
for(n in 1:length(webElems)){
  if(webElems[[n]][1] == 'Male')
    disease_covariates[1,n] = 1

  if(webElems[[n]][2] == 'Mixed Latino')
    disease_covariates[2,n] = 1
  else if(webElems[[n]][2] == 'Puerto Rican')
    disease_covariates[3,n] = 1
  else if(webElems[[n]][2] == 'Other Latino')
    disease_covariates[4,n] = 1
}

DGP = list()
DGP$Ometh = as.matrix(data)
DGP$disease_covariates = as.matrix(disease_covariates)
rm(data); gc()

for(q in 2:4){
  subDir = paste0('Q', q)
  if (file.exists(subDir)){
    setwd(subDir)
  } else {
    dir.create(subDir)
    setwd(subDir)
  }

  for(k in 2:7){
    subDir = paste0('K', k)
    if (file.exists(subDir)){
      setwd(subDir)
    } else {
      dir.create(subDir)
      setwd(subDir)
    }

    Init = INIT(DGP, k, q, 'Humans', 0, 0, 123)
    # rm(list=setdiff(ls(), c('DGP', 'Init')))
    save.image("data.RData")
    setwd('../')
  }
  setwd('../')
}

#############################
# Real data analysis GDC
#############################
library(data.table)
library(sva)
library(meffil)
library(rvest)
library(RSelenium)

setwd('C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/GDC/Breast/')

dir.list = unlist(read.table('FolderList.txt'))

setwd('raw_data/')
setwd(paste0('./', dir.list[1]))
file.name = list.files("./", full.names = TRUE, pattern=".txt")
data = as.data.frame(fread(file.name))
data = na.omit(data)
gc()
setwd('../')

for(n in 2:length(dir.list)){
  setwd(paste0('./', dir.list[n]))
  file.name = list.files("./", full.names = TRUE, pattern=".txt")
  if(length(file.name) != 1)
    file.name = file.name[-which(file.name == "./annotations.txt")]
  one.data = as.data.frame(fread(file.name))
  one.data = na.omit(one.data)
  data = merge(data, one.data, by = 'V1')
  rm(one.data); gc()
  print(dim(data)[1])
  setwd('../')
}

# write.table(data, 'raw_breast_data.txt', quote = F, row.names = F, col.names = F)

data = as.data.frame(fread('breast_data.txt'))
row.names(data) = data[,1]
data = data[,-1]
data = data[-which(apply(data, 1, mean) < .2 | apply(data, 1, mean) > .8),]

mVG = meffil.most.variable.cpgs(as.matrix(data), n = 10000, autosomal = F)
data = data[row.names(data) %in% mVG,]

DGP = list()
DGP$Ometh = as.matrix(data)
# DGP$disease_covariates = as.matrix(disease_covariates)
rm(data); gc()

for(q in 2:4){
  subDir = paste0('Q', q)
  if (file.exists(subDir)){
    setwd(subDir)
  } else {
    dir.create(subDir)
    setwd(subDir)
  }

  for(k in 3:10){
    subDir = paste0('K', k)
    if (file.exists(subDir)){
      setwd(subDir)
    } else {
      dir.create(subDir)
      setwd(subDir)
    }

    Init = INIT(DGP, k, q, 'Humans', 0, 0, 123)
    # rm(list=setdiff(ls(), c('DGP', 'Init')))
    save.image("data.RData")
    setwd('../')
  }
  setwd('../')
}


#############################
# Luo's simulation
#############################
luo_simulation = function(K, Q){

  set.seed(123)

  N = 180     #number of samples
  n1 = 60    #number of controls
  n2 = 120    #number of cases
  M = 2000   #number of CpG sites
  # K = 3       #underlying cell type number

  ###simulate methylation baseline profiles
  #assume cell type 1 and cell type 2 are from the same lineage
  #cell type 1
  methy1 = rbeta(M,3,6)
  #cell type 2
  methy2 = methy1 + rnorm(M, sd=0.01)
  ind = sample(1:M, M/5)
  methy2[ind] = rbeta(length(ind),3,6)

  #cell type 3
  methy3 = rbeta(M,3,6)
  Mu = cbind(methy1, methy2, methy3)

  #number of disease_covariates
  # Q = 2

  ###simulate disease_covariates / phenotype (disease status and age)
  disease_covariates = rbind(c(rep(0, n1),rep(1, n2)), runif(N, min=20, max=50))

  ###simulate phenotype effects
  Gamma = array(0, dim=c(M,K,Q))

  #control vs case
  m_common = 10
  max_signal = 0.15
  min_signal = 0.07

  #we allow different signs and magnitudes
  signs = sample(c(-1,1), m_common*K, replace=TRUE)
  Gamma[1:m_common,1:K,1] = signs * runif(m_common*K, min=min_signal, max=max_signal)

  m_seperate = 10
  signs = sample(c(-1,1), m_seperate*2, replace=TRUE)
  Gamma[m_common+(1:m_seperate),1:2,1] = signs *
    runif(m_seperate*2, min=min_signal, max=max_signal)

  signs = sample(c(-1,1), m_seperate, replace=TRUE)
  Gamma[m_common+m_seperate+(1:m_seperate),K,1] = signs *
    runif(m_seperate, min=min_signal, max=max_signal)

  #age
  base = 20
  m_common = 10
  max_signal = 0.015
  min_signal = 0.007
  signs = sample(c(-1,1), m_common*K, replace=TRUE)
  Gamma[base+1:m_common,1:K,2] = signs *
    runif(m_common*K, min=min_signal, max=max_signal)

  m_seperate = 10
  signs = sample(c(-1,1), m_seperate*2, replace=TRUE)
  Gamma[base+m_common+(1:m_seperate),1:2,2] = signs *
    runif(m_seperate*2, min=min_signal, max=max_signal)

  signs = sample(c(-1,1), m_seperate, replace=TRUE)
  Gamma[base+m_common+m_seperate+(1:m_seperate),K,2] = signs *
    runif(m_seperate, min=min_signal, max=max_signal)

  ###generate the cellular compositions
  P = sapply(1:N, function(i){
    if(disease_covariates[1,i]==0){ #if control
      rDirichlet(c(4,4, 2+disease_covariates[2,i]/10), 1)
    }else{
      rDirichlet(c(4,4, 5+disease_covariates[2,i]/10), 1)
    }
  })

  ###generate the observed methylation profiles
  Ometh = NULL
  for(i in 1:N){
    utmp = t(sapply(1:M, function(j){
      tmp1 = colSums(disease_covariates[ ,i] * t(Gamma[j, , ]))
      rnorm(K,mean=Mu[j, ]+tmp1,sd=0.01)
    }))
    tmp2 = colSums(P[ ,i] * t(utmp))
    Ometh = cbind(Ometh, tmp2 + rnorm(M, sd = 0.01))
  }

  #constrain methylation values between 0 and 1
  Ometh[Ometh > 1] = 1
  Ometh[Ometh < 0] = 0

  return(list('Ometh' = Ometh, 'disease_covariates' = disease_covariates))
}

library(HIREcewas)

K = 3; Q = 2
DGP = luo_simulation(K, Q)
N = dim(DGP$Ometh)[2]; M = dim(DGP$Ometh)[1]
Init = INIT(DGP, 3, 2, 'Humans', 0, 0, 123)
time0 = proc.time()
ret_list = HIREcewas::HIRE(DGP$Ometh, DGP$disease_covariates, Init$Mu, Init$Gamma, Init$P_t, iter=1000)
print(proc.time() - time0)

message("  Calculating p-values...\n")
tmp = NULL
for(ell in seq_len(Q)){
  tmp = cbind(tmp, DGP$disease_covariates[ell, ]*t(ret_list$P_t))
}
x_matr = cbind( tmp, t(ret_list$P_t)[, seq(2,K)])
x_matr = as.matrix(x_matr)

pvalues = t(vapply(seq_len(M), function(m){
  y_vec = DGP$Ometh[m,]
  fit = lm(y_vec~x_matr)
  summary(fit)$coef[seq(2, (1+Q*K)),4]
}, FUN.VALUE = rep(-1,Q*K)))
message("  Done!\n")
ret_list[[7]] = pvalues
names(ret_list)[7] = "pvalues"
names(ret_list)[6] = "pBIC"
d = (K-1)*N + M*(1+2*K+Q*K) #number of parameters
d0 = sum(ret_list$pvalues > alpha/(Q*M*K)) #number of zero parameters
ret_list[[6]] = ret_list[[6]] - log(N) * d0
















