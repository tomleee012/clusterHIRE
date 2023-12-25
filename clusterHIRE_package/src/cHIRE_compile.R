library(cHIRE0)

load('data.RData')
#############################
# Parameters
#############################
nc = 16

time0 = proc.time()
ret_list=cHIRE0::cHIRE(DGP$Ometh, DGP$Covariates, Init$EZ, Init$Pi_t, Init$Mu, Init$Beta,
                       Init$Gamma, Init$P_t, iter=1000, tol = 1e-20, core=nc)
print(proc.time() - time0)

rm(list=setdiff(ls(), c('ret_list')))
save.image('cHIRE0_.RData')