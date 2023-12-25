#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

# include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cHIRE_v6.h"

extern SEXP clusterHIRERcallC(SEXP);

static const R_CallMethodDef CallEntries[] = {
	{"clusterHIRERcallC",          (DL_FUNC)&clusterHIRERcallC,          1},
	{NULL, NULL, 0}
};

void R_init_cHIREewas(DllInfo* dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}

SEXP getListElement(SEXP list, char* str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	short int i;

	for (i = 0; i < length(list); i++) {
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}

	return elmt;
}

SEXP clusterHIRERcallC(SEXP args) {
	short int nProtected = 0;
	float** P_t, ** Mu, *** Beta, *** Gamma, * Eps_M, ** Eps_MK, ** Ometh, ** Covariates, ** E_Z, * Pi_t;
	double BIC = -1;
	double* BICpointer;
	float tol;
	BICpointer = &BIC;
	short int K, N, M, Q, P, num_iteration, nc = 32;

	SEXP list, Pi_init, P_init, Mu_init, Beta_init, Beta_init_dim, Gamma_init, Gamma_init_dim, \
		Eps_MK_init, Eps_M_init, Ometh_r, Ometh_r_dim, Covariate_r, Z_r, tol_r, num_iter, num_core;

	//input variables from R
	list = args;
	Pi_init = getListElement(list, "Pi_init");
	P_init = getListElement(list, "P_init");
	Mu_init = getListElement(list, "Mu_init");
	Beta_init = getListElement(list, "Beta_init");
	Beta_init_dim = getListElement(list, "Beta_init_dim");
	Gamma_init = getListElement(list, "Gamma_init");
	Gamma_init_dim = getListElement(list, "Gamma_init_dim");
	Ometh_r = getListElement(list, "Ometh_r");
	Covariate_r = getListElement(list, "Covariate_r");
	Ometh_r_dim = getListElement(list, "Ometh_r_dim");
	Z_r = getListElement(list, "Z_r");
	Eps_MK_init = getListElement(list, "Eps_MK_init");
	Eps_M_init = getListElement(list, "Eps_M_init");
	tol_r = getListElement(list, "tol_r");
	num_iter = getListElement(list, "num_iter");
	num_core = getListElement(list, "num_core");

	M = INTEGER(Ometh_r_dim)[0];  // CpG site number
	N = INTEGER(Ometh_r_dim)[1];  // sample number
	K = INTEGER(Gamma_init_dim)[1];// cell number
	Q = INTEGER(Gamma_init_dim)[2];      // disease subtype number
	P = INTEGER(Beta_init_dim)[2];      // covariate number
	tol = REAL(tol_r)[0]; 		  // tolerance
	num_iteration = INTEGER(num_iter)[0]; // maximum iteration number
	nc = INTEGER(num_core)[0]; // number of cores

	//declare and initialize variables in C
	Ometh = make2Darray(M, N);
	Covariates = make2Darray(N, P);
	E_Z = make2Darray(N, Q);
	P_t = make2Darray(K, N);
	Mu = make2Darray(M, K);
	Beta = make3Darray(M, K, P);
	Gamma = make3Darray(M, K, Q);
	Eps_MK = make2Darray(M, K);
	Eps_M = (float*)malloc(M * sizeof(float));
	Pi_t = (float*)malloc(Q * sizeof(float));

	short int n, m, k, q, p;

	for (q = 0; q < Q; q++) {
		Pi_t[q] = REAL(Pi_init)[q];
	}

	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++) {
			Ometh[m][n] = REAL(Ometh_r)[m + M * n];
		}

		for (p = 0; p < P; p++) {
			Covariates[n][p] = REAL(Covariate_r)[n + N * p];
		}

		for (q = 0; q < Q; q++) {
			E_Z[n][q] = REAL(Z_r)[n + N * q];

		}

		for (k = 0; k < K; k++) {
			P_t[k][n] = REAL(P_init)[k + K * n];
		}
	}

	for (m = 0; m < M; m++) {
		Eps_M[m] = REAL(Eps_M_init)[m];
		for (k = 0; k < K; k++) {
			Mu[m][k] = REAL(Mu_init)[m + M * k];
			Eps_MK[m][k] = REAL(Eps_MK_init)[m + M * k];
			for (q = 0; q < Q; q++) {
				Gamma[m][k][q] = REAL(Gamma_init)[m + M * k + M * K * q];
			}
			for (p = 0; p < P; p++) {
				Beta[m][k][p] = REAL(Beta_init)[m + M * k + M * K * p];
			}
		}
	}

	// Set the number of cores for parallel
	omp_set_num_threads(nc);
	clusterHIREAlgo(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
		Ometh, Covariates, E_Z, BICpointer, tol, num_iteration);

	//return Pi_t
	SEXP ret_pi_t;
	PROTECT(ret_pi_t = allocVector(REALSXP, Q));
	++nProtected;
	for (q = 0; q < Q; q++) {
		REAL(ret_pi_t)[q] = Pi_t[q];
	}

	//return P_t
	SEXP ret_P_t, ret_P_t_dim;
	PROTECT(ret_P_t = allocVector(REALSXP, K * N));
	++nProtected;
	for (k = 0; k < K; k++) {
		for (n = 0; n < N; n++) {
			REAL(ret_P_t)[k + K * n] = P_t[k][n];
		}
	}
	PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_P_t_dim)[0] = K;
	INTEGER(ret_P_t_dim)[1] = N;
	setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);

	//return Mu
	SEXP ret_mu_t, ret_mu_t_dim;
	PROTECT(ret_mu_t = allocVector(REALSXP, M * K));
	++nProtected;
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			REAL(ret_mu_t)[m + M * k] = Mu[m][k];
		}
	}
	PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_mu_t_dim)[0] = M;
	INTEGER(ret_mu_t_dim)[1] = K;
	setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);

	//return Beta
	SEXP ret_beta_t, ret_beta_t_dim;
	PROTECT(ret_beta_t = allocVector(REALSXP, M * K * P));
	++nProtected;
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			for (p = 0; p < P; p++) {
				REAL(ret_beta_t)[m + M * k + M * K * p] = Beta[m][k][p];
			}
		}
	}
	PROTECT(ret_beta_t_dim = allocVector(INTSXP, 3));
	++nProtected;
	INTEGER(ret_beta_t_dim)[0] = M;
	INTEGER(ret_beta_t_dim)[1] = K;
	INTEGER(ret_beta_t_dim)[2] = P;
	setAttrib(ret_beta_t, R_DimSymbol, ret_beta_t_dim);

	//return Gamma
	SEXP ret_gamma_t, ret_gamma_t_dim;
	PROTECT(ret_gamma_t = allocVector(REALSXP, M * K * Q));
	++nProtected;
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			for (q = 0; q < Q; q++) {
				REAL(ret_gamma_t)[m + M * k + M * K * q] = Gamma[m][k][q];
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
	for (m = 0; m < M; m++) {
		REAL(ret_epsm_t)[m] = Eps_M[m];
	}

	//return Eps_MK
	SEXP ret_epsmk_t, ret_epsmk_t_dim;
	PROTECT(ret_epsmk_t = allocVector(REALSXP, M * K));
	++nProtected;
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			REAL(ret_epsmk_t)[m + M * k] = Eps_MK[m][k];
		}
	}
	PROTECT(ret_epsmk_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_epsmk_t_dim)[0] = M;
	INTEGER(ret_epsmk_t_dim)[1] = K;
	setAttrib(ret_epsmk_t, R_DimSymbol, ret_epsmk_t_dim);

	//return E_Z
	SEXP ret_eZ_t, ret_eZ_t_dim;
	PROTECT(ret_eZ_t = allocVector(REALSXP, N * Q));
	++nProtected;
	for (n = 0; n < N; n++) {
		for (q = 0; q < Q; q++) {
			REAL(ret_eZ_t)[n + N * q] = E_Z[n][q];
		}
	}
	PROTECT(ret_eZ_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_eZ_t_dim)[0] = N;
	INTEGER(ret_eZ_t_dim)[1] = Q;
	setAttrib(ret_eZ_t, R_DimSymbol, ret_eZ_t_dim);

	//return BIC
	SEXP ret_BIC;
	PROTECT(ret_BIC = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_BIC)[0] = BIC;

	//return list
	SEXP retList;
	PROTECT(retList = allocVector(VECSXP, 9));
	++nProtected;

	SEXP names; // components names in the return list
	PROTECT(names = allocVector(STRSXP, 9));
	++nProtected;

	SET_STRING_ELT(names, 0, mkChar("Pi_t"));
	SET_STRING_ELT(names, 1, mkChar("P_t"));
	SET_STRING_ELT(names, 2, mkChar("Mu"));
	SET_STRING_ELT(names, 3, mkChar("Beta"));
	SET_STRING_ELT(names, 4, mkChar("Gamma"));
	SET_STRING_ELT(names, 5, mkChar("Eps_M"));
	SET_STRING_ELT(names, 6, mkChar("Eps_MK"));
	SET_STRING_ELT(names, 7, mkChar("E_Z"));
	SET_STRING_ELT(names, 8, mkChar("BIC"));

	//add elements to the return list
	SET_VECTOR_ELT(retList, 0, ret_pi_t);
	SET_VECTOR_ELT(retList, 1, ret_P_t);
	SET_VECTOR_ELT(retList, 2, ret_mu_t);
	SET_VECTOR_ELT(retList, 3, ret_beta_t);
	SET_VECTOR_ELT(retList, 4, ret_gamma_t);
	SET_VECTOR_ELT(retList, 5, ret_epsm_t);
	SET_VECTOR_ELT(retList, 6, ret_epsmk_t);
	SET_VECTOR_ELT(retList, 7, ret_eZ_t);
	SET_VECTOR_ELT(retList, 8, ret_BIC);
	setAttrib(retList, R_NamesSymbol, names);

	UNPROTECT(nProtected);
	delet2Darray(E_Z, N, Q);
	delet2Darray(P_t, K, N);
	delet2Darray(Mu, M, K);
	delet3Darray(Beta, M, K, P);
	delet3Darray(Gamma, M, K, Q);
	delet2Darray(Ometh, M, N);
	delet2Darray(Covariates, N, P);
	delet2Darray(Eps_MK, M, K);
	free(Eps_M);
	free(Pi_t);

	return retList;
}