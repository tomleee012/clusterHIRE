cHIRE = function(Ometh, Covariates, EZ, Pi_t, mu_t, beta_t, gamma_t, P_t,
                       tol = 1e-8, iter=1000, core=32, reshape=0){

  clusterHIRERcallC =
    function(Ometh, Covariates, EZ, Pi_t, P_t, mu_t, beta_t, gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core){
    args = list("Pi_init"=as.numeric(Pi_t), "P_init"=as.numeric(P_t), "Mu_init"=as.numeric(mu_t),
                "Beta_init" = as.numeric(beta_t), "Beta_init_dim" = as.integer(dim(beta_t)),
                "Gamma_init"=as.numeric(gamma_t), "Gamma_init_dim"=as.integer(dim(gamma_t)),
                "Ometh_r"=as.numeric(Ometh), "Ometh_r_dim"=as.integer(dim(Ometh)),
                "Covariate_r" = as.numeric(Covariates),
                "Z_r"=as.numeric(EZ), "Eps_MK_init"=as.numeric(Eps_MK_t),
                "Eps_M_init"=as.numeric(Eps_M_t), "tol_r" = as.numeric(tol),
                "num_iter" = as.integer(iter), 'num_core' = as.integer(core))
    ret_list = .Call("clusterHIRERcallC", args)
    return(ret_list)
  }

  if(dim(Ometh)[2] > dim(Ometh)[1]){
    Ometh = t(Ometh)
  }
  
  M = dim(Ometh)[1] #CpG site number
  N = dim(Ometh)[2] #sample number
  K = dim(gamma_t)[2]
  Q = dim(gamma_t)[3]
  P = dim(beta_t)[3]
  
  if(dim(Covariates)[1] != N){
    Covariates = t(Covariates)
  }
  if(dim(EZ)[1] != N){
    EZ = t(EZ)
  }
  if(dim(P_t)[2] != N){
    P_t = t(P_t)
  }

  Eps_MK_t = matrix(.0001, M, K)
  Eps_M_t = rep(.0001, M)

  if(reshape){
    Beta_t = array(0, dim = c(M, K, P))
    Gamma_t = array(0, dim=c(M, K, Q))
    for(p in 1:P)
      Beta_t[,,p] = beta_t[,(1 + K*(p-1)):(K*p)]
    for(q in 1:Q)
      Gamma_t[,,q] = gamma_t[,(1 + K*(q-1)):(K*q)]
  } else{
    Gamma_t = gamma_t
    Beta_t = beta_t
  }

  # message("  Initialization Done.\n")

  message("  Implementing EM algorithm... \n")
  ret_list = clusterHIRERcallC(Ometh, Covariates, EZ, Pi_t, P_t, mu_t, Beta_t,
                               Gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core)
  message("  Done! \n")

  return(ret_list)
}
