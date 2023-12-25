#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MAXX(a, b) (((a) > (b)) ? (a) : (0)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define ABS(a, b)  (((a) > (b)) ? (a-b) : (b-a)) 
#define SIGN(a)  (((a) > (0)) ? (1) : (-1)) 

//========================================================================================================================
//Define multiple dimensional arrays (Memory is not continuous here)
//========================================================================================================================

double**** make4Darray(int a1, int a2, int a3, int a4) {
	double**** tmp;
	tmp = (double****)malloc(a1 * sizeof(double***));
	for (int n = 0; n < a1; n++) {
		tmp[n] = (double***)malloc(a2 * sizeof(double**));
		for (int m = 0; m < a2; m++) {
			tmp[n][m] = (double**)malloc(a3 * sizeof(double*));
			for (int k = 0; k < a3; k++) {
				tmp[n][m][k] = (double*)malloc(a4 * sizeof(double));
			}
		}
	}
	return tmp;
}

void delet4Darray(double**** tmp, int a1, int a2, int a3, int a4) {
	for (int n = 0; n < a1; n++) {
		for (int m = 0; m < a2; m++) {
			for (int k = 0; k < a3; k++) {
				free(tmp[n][m][k]);
			}
			free(tmp[n][m]);
		}
		free(tmp[n]);
	}
	free(tmp);
}

double*** make3Darray(int a1, int a2, int a3) {
	double*** tmp;
	tmp = (double***)malloc(a1 * sizeof(double**));
	for (int n = 0; n < a1; n++) {
		tmp[n] = (double**)malloc(a2 * sizeof(double*));
		for (int m = 0; m < a2; m++) {
			tmp[n][m] = (double*)malloc(a3 * sizeof(double));
		}
	}
	return tmp;
}

void delet3Darray(double*** tmp, int a1, int a2, int a3) {
	for (int n = 0; n < a1; n++) {
		for (int m = 0; m < a2; m++) {
			free(tmp[n][m]);
		}
		free(tmp[n]);
	}
	free(tmp);
}

double** make2Darray(int a1, int a2) {
	double** tmp;
	tmp = (double**)malloc(a1 * sizeof(double*));
	for (int n = 0; n < a1; n++) {
		tmp[n] = (double*)malloc(a2 * sizeof(double));
	}
	return tmp;
}

void delet2Darray(double** tmp, int a1, int a2) {
	for (int n = 0; n < a1; n++) {
		free(tmp[n]);
	}
	free(tmp);
}

//========================================================================================================================
//calculate the inverse of a matrix (the code is from http://www.sourcecodesworld.com/source/show.asp?ScriptID=1086)
//========================================================================================================================
void inverse(double** A, int K, double** I) {
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < K; j++) {
			if (i == j)
				I[i][j] = 1;
			else
				I[i][j] = 0;
		}
	}

	/*---------------LoGiC starts here------------------*/		//procedure to make the matrix A to unit matrix
	for (int k = 0; k < K; k++) {
		double temp = A[k][k];
		for (int j = 0; j < K; j++) {
			A[k][j] /= temp;
			I[k][j] /= temp;
		}
		for (int i = 0; i < K; i++) {
			temp = A[i][k];
			for (int j = 0; j < K; j++) {
				if (i == k) break;
				A[i][j] -= A[k][j] * temp;
				I[i][j] -= I[k][j] * temp;
			}
		}
	}
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double** Dmat, double* dvec, double* xvec, int K) { //K is the dimension of xvec
	if (K == 1) {
		xvec[0] = 1;
	}
	else {
		double** Dmat_inv, * x_star, s1, s2;
		double Dmat_inv_sum = 0, x_star_sum = 0, * Dmat_inv_rowsum;
		int num_negatives = 0;
		double* ind_negative;
		Dmat_inv = make2Darray(K, K);
		x_star = (double*)malloc(K * sizeof(double));
		Dmat_inv_rowsum = (double*)malloc(K * sizeof(double));
		ind_negative = (double*)malloc(K * sizeof(double));
		inverse(Dmat, K, Dmat_inv);
		for (int k = 0; k < K; k++) {
			s1 = 0;
			s2 = 0;
			for (int k1 = 0; k1 < K; k1++) {
				s1 += Dmat_inv[k][k1] * dvec[k1];
				Dmat_inv_sum += Dmat_inv[k][k1];
				s2 += Dmat_inv[k][k1];
			}
			x_star[k] = s1;
			Dmat_inv_rowsum[k] = s2;
			x_star_sum += s1;
		}
		for (int k = 0; k < K; k++) {
			xvec[k] = x_star[k] + (1 - x_star_sum) / Dmat_inv_sum * Dmat_inv_rowsum[k];
			if (xvec[k] < 0) {
				num_negatives++;
				ind_negative[k] = 1;
			}
			else {
				ind_negative[k] = 0;
			}
		}
		free(x_star);
		free(Dmat_inv_rowsum);
		delet2Darray(Dmat_inv, K, K);

		if (num_negatives == 0) {
			free(ind_negative);
		}
		else {
			int Knew = K - num_negatives, i, j;
			double** Dmat_new, * dvec_new, * xvec_sub;
			Dmat_new = make2Darray(Knew, Knew);
			dvec_new = (double*)malloc((Knew) * sizeof(double));
			xvec_sub = (double*)malloc((Knew) * sizeof(double));
			i = 0;

			for (int k1 = 0; k1 < K; k1++) {
				if (ind_negative[k1] == 0) {
					dvec_new[i] = dvec[k1];
					j = 0;
					for (int k2 = 0; k2 < K; k2++) {
						if (ind_negative[k2] == 0) {
							Dmat_new[i][j] = Dmat[k1][k2];
							j++;
						}
					}
					i++;
				}
				else {
					xvec[k1] = 0;
				}
			}

			quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
			i = 0;
			for (int k = 0; k < K; k++) {
				if (ind_negative[k] == 0) {
					xvec[k] = xvec_sub[i];
					i++;
				}
			}
			free(dvec_new);
			free(xvec_sub);
			delet2Darray(Dmat_new, Knew, Knew);
		}
	}

}

//the objective function when updating P_t
double val2(double** P_t, double* Eps_M, int num_celltypes, int num_CpGsites, int num_cancer_subtypes,
	double** Ometh, double**** E_sigma, double**** E_mu, double** E_Z, int n) {
	int K, M, Q;
	K = num_celltypes;
	M = num_CpGsites;
	Q = num_cancer_subtypes;

	double s1 = 0, s2 = 0, s3 = 0, sum = 0;
	int k, k2, m, q;

#pragma omp parallel for reduction(+:s1,s2,s3,sum) 
	for (m = 0; m < M; m++) {
		s1 = 0;
		s3 = 0;
		for (k = 0; k < K; k++) {
			for (k2 = 0; k2 < K; k2++) {
				s1 += P_t[k][n] * E_sigma[n][m][k][k2] * P_t[k2][n];
			}
		}

		for (q = 0; q < Q; q++) {
			s2 = 0;
			for (k = 0; k < K; k++) {
				s2 += E_mu[n][m][q][k] * P_t[k][n];
			}
			s2 = E_Z[n][q] * pow(Ometh[m][n] - s2, 2);
			s3 += s2;
		}
		sum += (s1 + s3) / (2 * Eps_M[m]);
	}

	return sum;
}

//========================================================================================================================
//functions in the algorithm
//========================================================================================================================
double normal_density_log(double x, double mean, double sd) {
	return -0.5 * log(2 * M_PI) - log(sd) - pow(x - mean, 2) / (2 * sd * sd);
}

double normal_density(double x, double mean, double var) {
	return exp(-pow(x - mean, 2) / (2 * var)) / sqrt(2 * M_PI * var);
}

// update P_t
void Update_P(double** P_t, double* Eps_M, int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes,
	double** Ometh, double**** E_sigma, double**** E_mu, double** E_Z) {

	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	double s1, s2;
	int k, k2, n, m, q;

	double tar2a, tar2b;
	double* xvec = (double*)malloc(K * sizeof(double));
	double* xvec0 = (double*)malloc(K * sizeof(double));
	double* dvec = (double*)malloc(K * sizeof(double));
	double** Dmat = make2Darray(K, K);
	double** Dmat_inv = make2Darray(K, K);
	double* tmp1 = (double*)malloc(Q * sizeof(double));
	double* tmp2 = (double*)malloc(Q * sizeof(double));

	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			for (k2 = k; k2 < K; k2++) {
				for (q = 0; q < Q; q++) {
					tmp1[q] = 0;
					tmp2[q] = 0;
					for (m = 0; m < M; m++) {
						tmp1[q] += (E_sigma[n][m][k][k2] + E_mu[n][m][q][k] * E_mu[n][m][q][k2]) / Eps_M[m];
						tmp2[q] += Ometh[m][n] * E_mu[n][m][q][k] / Eps_M[m];
					}
					tmp1[q] *= E_Z[n][q];
					tmp2[q] *= E_Z[n][q];
				}
				s1 = 0; s2 = 0;
				for (q = 0; q < Q; q++) {
					s1 += tmp1[q];
					s2 += tmp2[q];
				}

				Dmat[k][k2] = s1;
				if (k2 > k) {
					Dmat[k2][k] = Dmat[k][k2];
				}
			}
			dvec[k] = s2;
		}

		for (k = 0; k < K; k++) {
			xvec0[k] = P_t[k][n];
		}

		tar2a = val2(P_t, Eps_M, K, M, Q, Ometh, E_sigma, E_mu, E_Z, n);
		quadprog(Dmat, dvec, xvec, K);

		for (k = 0; k < K; k++) {
			P_t[k][n] = xvec[k];
		}

		tar2b = val2(P_t, Eps_M, K, M, Q, Ometh, E_sigma, E_mu, E_Z, n);

		if (tar2b > tar2a) {//avoid offset effects when estimates are stable
			for (k = 0; k < K; k++) {
				P_t[k][n] = xvec0[k];
			}
		}
	}

	free(xvec); free(xvec0); free(dvec); free(tmp1); free(tmp2);
	delet2Darray(Dmat, K, K);
	delet2Darray(Dmat_inv, K, K);
}

void Zstep(int n, int num_cancer_subtypes, double** E_Z, double* tmp) {

	int Q;
	Q = num_cancer_subtypes;
	double s1, s2 = 0;

	s1 = tmp[0];
	for (int q = 1; q < Q; q++)
		s1 = MAX(s1, tmp[q]);

	for (int q = 0; q < Q; q++) {
		s2 += exp(tmp[q] - s1);
	}

	for (int q = 0; q < Q; q++) {
		E_Z[n][q] = exp(tmp[q] - s1) / s2;
	}
}

//conduct E step
void EstepI(double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M,
	double** Eps_MK, int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes,
	int num_covariates, double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {

	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	double s1 = 0, s2 = 0;

#pragma omp parallel for reduction(+:s1,s2) schedule(static,1)
	for (int n = 0; n < N; n++) {
		for (int m = 0; m < M; m++) {
			// update E_sigma
			s1 = 0;
			for (int k = 0; k < K; k++) {
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			for (int k1 = 0; k1 < K; k1++) {
				for (int k2 = k1; k2 < K; k2++) {
					E_sigma[n][m][k1][k2] = -1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k2][n] * Eps_MK[m][k1] * Eps_MK[m][k2] / Eps_M[m];
					E_sigma[n][m][k2][k1] = E_sigma[n][m][k1][k2];
				}
				E_sigma[n][m][k1][k1] = Eps_MK[m][k1] - 1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k1][n] * Eps_MK[m][k1] * Eps_MK[m][k1] / Eps_M[m];
			}

			// update E_mu
			for (int q = 0; q < Q; q++) {
				for (int k = 0; k < K; k++) {
					s1 = 0;
					for (int k1 = 0; k1 < K; k1++) {
						s2 = 0;
						for (int p = 0; p < P; p++)
							s2 += Beta[m][k1][p] * Covariates[n][p];
						s1 += (Ometh[m][n] * P_t[k1][n] / Eps_M[m] + (Mu[m][k1] + s2 + Gamma[m][k1][q]) / Eps_MK[m][k1]) * E_sigma[n][m][k1][k];
					}
					E_mu[n][m][q][k] = s1;
				}
			}
		}
	}
}

void EstepII(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M,
	double** Eps_MK, int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes,
	int num_covariates, double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {

	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	double s1 = 0, s2 = 0, s3 = 0;
	double* tmp = (double*)malloc(Q * sizeof(double));

	// #pragma omp parallel for reduction(+:s1,s2) schedule(static,1)
	for (int n = 0; n < N; n++) {
		for (int q = 0; q < Q; q++)
			tmp[q] = Pi_t[q];
		for (int m = 0; m < M; m++) {
			// update E_sigma
			s1 = 0;
			for (int k = 0; k < K; k++) {
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			for (int k1 = 0; k1 < K; k1++) {
				for (int k2 = k1; k2 < K; k2++) {
					E_sigma[n][m][k1][k2] = -1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k2][n] * Eps_MK[m][k1] * Eps_MK[m][k2] / Eps_M[m];
					E_sigma[n][m][k2][k1] = E_sigma[n][m][k1][k2];
				}
				E_sigma[n][m][k1][k1] = Eps_MK[m][k1] - 1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k1][n] * Eps_MK[m][k1] * Eps_MK[m][k1] / Eps_M[m];
			}

			// update E_mu
			for (int q = 0; q < Q; q++) {
				for (int k = 0; k < K; k++) {
					s1 = 0;
					for (int k1 = 0; k1 < K; k1++) {
						s2 = 0;
						for (int p = 0; p < P; p++)
							s2 += Beta[m][k1][p] * Covariates[n][p];
						s1 += (Ometh[m][n] * P_t[k1][n] / Eps_M[m] + (Mu[m][k1] + s2 + Gamma[m][k1][q]) / Eps_MK[m][k1]) * E_sigma[n][m][k1][k];
					}
					E_mu[n][m][q][k] = s1;
				}
			}

			// update E_Z
			s1 = Eps_M[m];
			for (int k = 0; k < K; k++)
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];

			for (int q = 0; q < Q; q++) {
				s2 = 0;
				for (int k = 0; k < K; k++) {
					s3 = 0;
					for (int p = 0; p < P; p++)
						s3 += Beta[m][k][p] * Covariates[n][p];
					s2 += (Mu[m][k] + s3 + Gamma[m][k][q]) * P_t[k][n];
				}
				tmp[q] += normal_density_log(Ometh[m][n], s2, sqrt(s1));
			}
		}
		Zstep(n, num_cancer_subtypes, E_Z, tmp);
	}
	free(tmp);
}

// conduct M step
void MstepI(double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M, double** Eps_MK,
	int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes, int num_covariates,
	double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {
	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double s1 = 0, s2 = 0, s3 = 0;

		//==============================
		// coordinate descent algorithm
		//==============================
#pragma omp parallel for reduction(+:s1,s2,s3)
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - s2 - Gamma[m][k][q]) * E_Z[n][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (int n = 0; n < N; n++) {
				s1 = 0;
				for (int p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (int q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][q][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_Z[n][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k = 0; k < K; k++) {
				for (int k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (int q = 0; q < Q; q++) {
				s2 = 0;
				for (int k = 0; k < K; k++)
					s2 += E_mu[n][m][q][k] * P_t[k][n];
				s3 += (s1 + pow(Ometh[m][n] - s2, 2)) * E_Z[n][q];
			}
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, num_cancer_subtypes, Ometh, E_sigma, E_mu, E_Z);
}

// conduct M step
void MstepII(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M, double** Eps_MK,
	int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes, int num_covariates,
	double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {
	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double s1 = 0, s2 = 0, s3 = 0;

		//==============================
		// update Pi_t
		//==============================
	#pragma omp parallel for reduction(+:s1)
		for (int q = 0; q < Q; q++) {
			s1 = 0;
			for (int n = 0; n < N; n++) {
				s1 += E_Z[n][q];
			}
			Pi_t[q] = s1 / N;
		}

		//==============================
		// coordinate descent algorithm
		//==============================
#pragma omp parallel for reduction(+:s1,s2,s3)
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - s2 - Gamma[m][k][q]) * E_Z[n][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));


			//==============================
			// update Beta
			//==============================
			for (int p = 0; p < P; p++) {
				s1 = 0;
				s3 = 0;
				for (int n = 0; n < N; n++) {
					s2 = 0;
					for (int p2 = 0; p2 < P; p2++) {
						if (p2 == p) {
							continue;
						}
						s2 += Beta[m][k][p2] * Covariates[n][p2];
					}
					for (int q = 0; q < Q; q++) {
						s1 += (E_mu[n][m][q][k] - Mu[m][k] - s2 - Gamma[m][k][q]) * Covariates[n][p] * E_Z[n][q];
					}
					s3 += (Covariates[n][p] * Covariates[n][p]);
				}
				Beta[m][k][p] = s1 / s3;
				//for (int q = 0; q < Q; q++) {
				//	Beta[m][k][p] = SIGN(s1 / s3) * MAXX(ABS(s1 / s3, 0), ABS(Gamma[m][k][q], 0));
				//}
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (int n = 0; n < N; n++) {
				s1 = 0;
				for (int p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (int q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][q][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_Z[n][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k = 0; k < K; k++) {
				for (int k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (int q = 0; q < Q; q++) {
				s2 = 0;
				for (int k = 0; k < K; k++)
					s2 += E_mu[n][m][q][k] * P_t[k][n];
				s3 += (s1 + pow(Ometh[m][n] - s2, 2)) * E_Z[n][q];
			}
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	// Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, num_cancer_subtypes, Ometh, E_sigma, E_mu, E_Z);
}

void MstepIII(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M, double** Eps_MK,
	int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes, int num_covariates,
	double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {
	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double s1 = 0, s2 = 0, s3 = 0;

		//==============================
		// update Pi_t
		//==============================
	#pragma omp parallel for reduction(+:s1)
		for (int q = 0; q < Q; q++) {
			s1 = 0;
			for (int n = 0; n < N; n++) {
				s1 += E_Z[n][q];
			}
			Pi_t[q] = s1 / N;
		}

		//==============================
		// coordinate descent algorithm
		//==============================
#pragma omp parallel for reduction(+:s1,s2,s3)
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - s2 - Gamma[m][k][q]) * E_Z[n][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			// update Gamma
			//==============================
			for (int q = 1; q < Q; q++) {
				s1 = 0;
				s3 = 0;
				for (int n = 0; n < N; n++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - Mu[m][k] - s2) * E_Z[n][q];
					s3 += E_Z[n][q];
				}
				Gamma[m][k][q] = s1 / s3;
				//for (int p = 0; p < P; p++) {
				//	Gamma[m][k][q] = SIGN(s1 / s3) * MAXX(ABS(s1 / s3, 0), ABS(Beta[m][k][p], 0));
				//}
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (int n = 0; n < N; n++) {
				s1 = 0;
				for (int p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (int q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][q][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_Z[n][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k = 0; k < K; k++) {
				for (int k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (int q = 0; q < Q; q++) {
				s2 = 0;
				for (int k = 0; k < K; k++)
					s2 += E_mu[n][m][q][k] * P_t[k][n];
				s3 += (s1 + pow(Ometh[m][n] - s2, 2)) * E_Z[n][q];
			}
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	// Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, num_cancer_subtypes, Ometh, E_sigma, E_mu, E_Z);
}


void MstepIV(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M, double** Eps_MK,
	int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes, int num_covariates,
	double** Ometh, double** Covariates, double**** E_sigma, double**** E_mu, double** E_Z) {
	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double s1 = 0, s2 = 0, s3 = 0;

	//==============================
	// update Pi_t
	//==============================
#pragma omp parallel for reduction(+:s1)
	for (int q = 0; q < Q; q++) {
		s1 = 0;
		for (int n = 0; n < N; n++) {
			s1 += E_Z[n][q];
		}
		Pi_t[q] = s1 / N;
	}

	//==============================
	// coordinate descent algorithm
	//==============================
#pragma omp parallel for reduction(+:s1,s2,s3)
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - s2 - Gamma[m][k][q]) * E_Z[n][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			// update Gamma
			//==============================
			for (int q = 1; q < Q; q++) {
				s1 = 0;
				s3 = 0;
				for (int n = 0; n < N; n++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][q][k] - Mu[m][k] - s2) * E_Z[n][q];
					s3 += E_Z[n][q];
				}
				Gamma[m][k][q] = s1 / s3;
				//for (int p = 0; p < P; p++) {
				//	Gamma[m][k][q] = SIGN(s1 / s3) * MAXX(ABS(s1 / s3, 0), ABS(Beta[m][k][p], 0));
				//}

			}

			//==============================
			// update Beta
			//==============================
			for (int p = 0; p < P; p++) {
				s1 = 0;
				s3 = 0;
				for (int n = 0; n < N; n++) {
					s2 = 0;
					for (int p2 = 0; p2 < P; p2++) {
						if (p2 == p) {
							continue;
						}
						s2 += Beta[m][k][p2] * Covariates[n][p2];
					}
					for (int q = 0; q < Q; q++) {
						s1 += (E_mu[n][m][q][k] - Mu[m][k] - s2 - Gamma[m][k][q]) * Covariates[n][p] * E_Z[n][q];
					}
					s3 += (Covariates[n][p] * Covariates[n][p]);
				}
				Beta[m][k][p] = s1 / s3;
				//for (int q = 0; q < Q; q++) {
				//	Beta[m][k][p] = SIGN(s1 / s3) * MAXX(ABS(s1 / s3, 0), ABS(Gamma[m][k][q], 0));
				//}
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (int n = 0; n < N; n++) {
				s1 = 0;
				for (int p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (int q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][q][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_Z[n][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k = 0; k < K; k++) {
				for (int k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (int q = 0; q < Q; q++) {
				s2 = 0;
				for (int k = 0; k < K; k++)
					s2 += E_mu[n][m][q][k] * P_t[k][n];
				s3 += (s1 + pow(Ometh[m][n] - s2, 2)) * E_Z[n][q];
			}
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	// Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, num_cancer_subtypes, Ometh, E_sigma, E_mu, E_Z);
}

//calculate the observed data likelihood
double observed_log_likelihood(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M,
	double** Eps_MK, int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes,
	int num_covariates, double** Ometh, double** Covariates) {

	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double sum = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	double ploglike = 0;

#pragma omp parallel for reduction(+:s1,s2,s3,s4,sum,ploglike)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			s1 = 0;
			s2 = 0;
			s3 = 0;
			sum = 0;

			for (int k = 0; k < K; k++) {
				s1 += P_t[k][n] * Mu[m][k];
				s3 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
				for (int p = 0; p < P; p++)
					s2 += P_t[k][n] * Beta[m][k][p] * Covariates[n][p];
			}

			for (int q = 0; q < Q; q++) {
				s4 = 0;
				for (int k = 0; k < K; k++) {
					s4 += P_t[k][n] * Gamma[m][k][q];
				}
				sum += Pi_t[q] * normal_density(Ometh[m][n], s1 + s2 + s4, s3 + Eps_M[m]);
			}
			ploglike += log(sum);
		}
	}
	return (-2.0) * ploglike + log(N) * (Q - 1 + (K - 1) * N + M * (1 + 2 * K + (Q - 1) * K + P * K));
}

//calculate the BIC to select an appropriate K
double BIC(double ploglike, int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes,
	int num_covariates) {
	int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	return (-2.0) * ploglike + log(N) * (Q - 1 + (K - 1) * N + M * (1 + 2 * K + (Q - 1) * K + P * K));
}

//carry out the EM algorithm
void clusterHIREAlgo(double* Pi_t, double** P_t, double** Mu, double*** Beta, double*** Gamma, double* Eps_M, double** Eps_MK,
	int num_celltypes, int num_samples, int num_CpGsites, int num_cancer_subtypes, int num_covariates, double** Ometh,
	double** Covariates, double** E_Z, double* BICpointer, double tol, int num_iteration) {

	int K, N, M, Q, P, iter = 0, num_neg = 0;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	double**** E_sigma = make4Darray(N, M, K, K);
	double**** E_mu = make4Darray(N, M, Q, K);
	double ploglike = -100000, ploglike1 = 0;

	while (iter < num_iteration * .1) {
		EstepI(P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		MstepI(P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
			K, N, M, Q, P, Ometh, Covariates);
		iter++;
		Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
		//Rprintf("Iteration: %d\n", iter);
	}

	//while (iter < num_iteration * .25) {
	//	EstepII(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
	//		Ometh, Covariates, E_sigma, E_mu, E_Z);
	//	MstepIII(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
	//		Ometh, Covariates, E_sigma, E_mu, E_Z);
	//	ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
	//		K, N, M, Q, P, Ometh, Covariates);
	//	iter++;
	//	Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
	//	//Rprintf("Iteration: %d\n", iter);
	//}

	while (iter < num_iteration * .4) {
		EstepII(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		MstepIV(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
			K, N, M, Q, P, Ometh, Covariates);
		iter++;
		Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
	}

	while (iter < num_iteration) {
		EstepI(P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		MstepIV(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z);
		ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
			K, N, M, Q, P, Ometh, Covariates);
		iter++;
		Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
	}

	//Rprintf("Iteration: %d\n", iter);
	//while (iter < num_iteration && (ABS(ploglike, ploglike1) / ABS(ploglike1, 0) > tol))
	// double ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh);
	*BICpointer = ploglike;
	delet4Darray(E_mu, N, M, Q, K);
	delet4Darray(E_sigma, N, M, K, K);
}