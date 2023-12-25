HIRE = function(Ometh, X, mu_t, gamma_t, P_t, tol = 10^(-5), iter=1000, core=32, start = 1, reshape=0){
  
  HIRERcallC = function(Ometh, X, P_t, mu_t, gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core, start){
    args = list("P_init"=as.numeric(P_t), "Mu_init"=as.numeric(mu_t), 
                "Gamma_init"=as.numeric(gamma_t), "Gamma_init_dim"=as.integer(dim(gamma_t)), 
                "Ometh_r"=as.numeric(Ometh), "Ometh_r_dim"=as.integer(dim(Ometh)), 
                "X_r"=as.numeric(X), "Eps_MK_init"=as.numeric(Eps_MK_t), 
                "Eps_M_init"=as.numeric(Eps_M_t), "tol_r" = as.numeric(tol), 
                "num_iter" = as.integer(iter),  "num_core" = as.integer(core),
                'start_r' = as.integer(start))
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
  ret_list = HIRERcallC(Ometh, X, P_t, mu_t, Gamma_t, Eps_MK_t, Eps_M_t, tol, iter, core, start)
  message("  Done! \n")
  
  return(ret_list)
}