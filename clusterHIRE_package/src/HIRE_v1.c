#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

# include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "HIRE_v1.h"


extern SEXP HIRERcallC(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"HIRERcallC",          (DL_FUNC) &HIRERcallC,          1},
    {NULL, NULL, 0}
};

void R_init_HIREewas(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

SEXP getListElement (SEXP list, char *str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for(i = 0; i < length(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}

	return elmt;
}

SEXP HIRERcallC(SEXP args) {
	int nProtected = 0;
	double **P_t,  **Mu, ***Gamma, *Eps_M, **Eps_MK, **Ometh, **X_scaled;
	double BIC=-1;
	double *BICpointer;
	double tol;
	BICpointer = &BIC;
	int K, N, M, Q, num_iteration, nc, start;
	
	SEXP list, P_init, Mu_init, Gamma_init, Gamma_init_dim, Eps_MK_init, Eps_M_init, Ometh_r, Ometh_r_dim, X_r, tol_r, num_iter, num_core, start_r;
	
	//input variables from R
	list = args;
	P_init = getListElement(list, "P_init");
	Mu_init = getListElement(list, "Mu_init");
	Gamma_init = getListElement(list, "Gamma_init");
	Gamma_init_dim = getListElement(list, "Gamma_init_dim");
	Ometh_r = getListElement(list, "Ometh_r");
	Ometh_r_dim = getListElement(list, "Ometh_r_dim");
	X_r = getListElement(list, "X_r");
	Eps_MK_init = getListElement(list, "Eps_MK_init");
	Eps_M_init = getListElement(list, "Eps_M_init");
	tol_r = getListElement(list, "tol_r");
	num_iter = getListElement(list, "num_iter");
	num_core = getListElement(list, "num_core");
	start_r = getListElement(list, "start_r");

	M = INTEGER(Ometh_r_dim)[0];  // CpG site number
	N = INTEGER(Ometh_r_dim)[1];  // sample number
	K = INTEGER(Gamma_init_dim)[1];// cell number
	Q = INTEGER(Gamma_init_dim)[2];      // covariate number
	tol = REAL(tol_r)[0]; 		  // tolerance
	num_iteration = INTEGER(num_iter)[0]; // maximum iteration number
	nc = INTEGER(num_core)[0]; // number of cores
	start = INTEGER(start_r)[0];

	//declare and initialize variables in C
	Ometh = make2Darray(M, N);
	X_scaled = make2Darray(Q, N);
	P_t = make2Darray(K, N);
	Mu = make2Darray(M, K);
	Gamma = make3Darray(M, K, Q);
	Eps_MK = make2Darray(M, K);
	Eps_M = (double *) malloc(M*sizeof(double));

	int n, m, k, q;
	
	for(n=0; n<N; n++){
		for(m=0; m<M; m++){
			Ometh[m][n] = REAL(Ometh_r)[m+M*n];
		}
		
		for(q=0; q<Q; q++){
			X_scaled[q][n] = REAL(X_r)[q+Q*n];
		}
		
		for(k=0; k<K; k++){
			P_t[k][n] = REAL(P_init)[k+K*n];
		}
	}
	
	for(m=0; m<M; m++){
		Eps_M[m] = REAL(Eps_M_init)[m];
		for(k=0; k<K; k++){
			Mu[m][k] = REAL(Mu_init)[m + M * k];
			Eps_MK[m][k] = REAL(Eps_MK_init)[m+M*k];
			for(q=0; q<Q; q++){
				Gamma[m][k][q] = REAL(Gamma_init)[m+M*k+M*K*q];
			}
		}
	}
	
	// Set the number of cores for parallel
	omp_set_num_threads(nc); 
	EmEwas(P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh, X_scaled, BICpointer, tol, num_iteration, start);

	//return P_t
	SEXP ret_P_t, ret_P_t_dim;
	PROTECT(ret_P_t = allocVector(REALSXP, K*N));
	++nProtected;
	for(k=0; k < K; k++){
		for(n=0; n < N; n++){
			REAL(ret_P_t)[k + K*n] = P_t[k][n];
		}
	}
	PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
	++nProtected;	
	INTEGER(ret_P_t_dim)[0] = K;
	INTEGER(ret_P_t_dim)[1] = N;
	setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);
	
	//return Mu
	SEXP ret_mu_t, ret_mu_t_dim;
	PROTECT(ret_mu_t = allocVector(REALSXP, M*K));
	++nProtected;
	for(m=0; m < M; m++){
		for(k=0; k < K; k++){
			REAL(ret_mu_t)[m + M*k] = Mu[m][k];
		}
	}
	PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
	++nProtected;	
	INTEGER(ret_mu_t_dim)[0] = M;
	INTEGER(ret_mu_t_dim)[1] = K;
	setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);
	
	//return Gamma
	SEXP ret_gamma_t, ret_gamma_t_dim;
	PROTECT(ret_gamma_t = allocVector(REALSXP, M*K*Q));
	++nProtected;
	for(m=0; m < M; m++){
		for(k=0; k < K; k++){
			for(q=0; q < Q; q++){
				REAL(ret_gamma_t)[m + M*k + M*K*q] = Gamma[m][k][q];
			}
		}
	}
	PROTECT(ret_gamma_t_dim = allocVector(INTSXP, 3));
	++nProtected;	
	INTEGER(ret_gamma_t_dim)[0] = M;
	INTEGER(ret_gamma_t_dim)[1] = K;
	INTEGER(ret_gamma_t_dim)[2] = Q;
	setAttrib(ret_gamma_t, R_DimSymbol, ret_gamma_t_dim);
	
	//return Eps_M
	SEXP ret_epsm_t;
	PROTECT(ret_epsm_t = allocVector(REALSXP, M));
	++nProtected;
	for(m=0; m < M; m++){
		REAL(ret_epsm_t)[m] = Eps_M[m];
	}

	//return Eps_MK
	SEXP ret_epsmk_t, ret_epsmk_t_dim;
	PROTECT(ret_epsmk_t = allocVector(REALSXP, M*K));
	++nProtected;
	for(m=0; m < M; m++){
		for(k=0; k<K; k++){
			REAL(ret_epsmk_t)[m+M*k] = Eps_MK[m][k];
		}
	}
	PROTECT(ret_epsmk_t_dim = allocVector(INTSXP, 2));
	++nProtected;	
	INTEGER(ret_epsmk_t_dim)[0] = M;
	INTEGER(ret_epsmk_t_dim)[1] = K;
	setAttrib(ret_epsmk_t, R_DimSymbol, ret_epsmk_t_dim);	
	
	//return BIC
	SEXP ret_BIC;
	PROTECT(ret_BIC = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_BIC)[0] = BIC;

	//return list
	SEXP retList; 
	PROTECT(retList = allocVector(VECSXP, 6));
	++nProtected;
	
	SEXP names; // components names in the return list
	PROTECT(names = allocVector(STRSXP, 6));
	++nProtected;
	
	SET_STRING_ELT(names, 0, mkChar("P_t"));
	SET_STRING_ELT(names, 1, mkChar("Mu"));
	SET_STRING_ELT(names, 2, mkChar("Gamma"));
	SET_STRING_ELT(names, 3, mkChar("Eps_M"));
	SET_STRING_ELT(names, 4, mkChar("Eps_MK"));
	SET_STRING_ELT(names, 5, mkChar("BIC"));
		
	//add elements to the return list
	SET_VECTOR_ELT(retList, 0, ret_P_t);
	SET_VECTOR_ELT(retList, 1, ret_mu_t);
	SET_VECTOR_ELT(retList, 2, ret_gamma_t);
	SET_VECTOR_ELT(retList, 3, ret_epsm_t);
	SET_VECTOR_ELT(retList, 4, ret_epsmk_t);
	SET_VECTOR_ELT(retList, 5, ret_BIC);
	setAttrib(retList, R_NamesSymbol, names);
	
	UNPROTECT(nProtected);
	delet2Darray(P_t, K, N);
	delet2Darray(Mu, M, K);
	delet2Darray(Ometh, M, N);
	delet2Darray(X_scaled, Q, N);
	delet2Darray(Eps_MK, M, K);
	free(Eps_M);
	delet3Darray(Gamma, M, K, Q);

	return retList;		     							
}