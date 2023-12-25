#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define AMAX(a, aa, b, ab)  (((a) > (b)) ? (aa) : (ab)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define ABS(a, b)  (((a) > (b)) ? (a-b) : (b-a)) 

//========================================================================================================================
//Define multiple dimensional arrays (Memory is not continuous here)
//========================================================================================================================

float**** make4Darray(short int a1, short int a2, short int a3, short int a4) {
	float**** tmp;
	tmp = (float****)malloc(a1 * sizeof(float***));
	for (short int n = 0; n < a1; n++) {
		tmp[n] = (float***)malloc(a2 * sizeof(float**));
		for (short int m = 0; m < a2; m++) {
			tmp[n][m] = (float**)malloc(a3 * sizeof(float*));
			for (short int k = 0; k < a3; k++) {
				tmp[n][m][k] = (float*)malloc(a4 * sizeof(float));
			}
		}
	}
	return tmp;
}

void delet4Darray(float**** tmp, short int a1, short int a2, short int a3, short int a4) {
	for (short int n = 0; n < a1; n++) {
		for (short int m = 0; m < a2; m++) {
			for (short int k = 0; k < a3; k++) {
				free(tmp[n][m][k]);
			}
			free(tmp[n][m]);
		}
		free(tmp[n]);
	}
	free(tmp);
}

float*** make3Darray(short int a1, short int a2, short int a3) {
	float*** tmp;
	tmp = (float***)malloc(a1 * sizeof(float**));
	for (short int n = 0; n < a1; n++) {
		tmp[n] = (float**)malloc(a2 * sizeof(float*));
		for (short int m = 0; m < a2; m++) {
			tmp[n][m] = (float*)malloc(a3 * sizeof(float));
		}
	}
	return tmp;
}

void delet3Darray(float*** tmp, short int a1, short int a2, short int a3) {
	for (short int n = 0; n < a1; n++) {
		for (short int m = 0; m < a2; m++) {
			free(tmp[n][m]);
		}
		free(tmp[n]);
	}
	free(tmp);
}

float** make2Darray(short int a1, short int a2) {
	float** tmp;
	tmp = (float**)malloc(a1 * sizeof(float*));
	for (short int n = 0; n < a1; n++) {
		tmp[n] = (float*)malloc(a2 * sizeof(float));
	}
	return tmp;
}

void delet2Darray(float** tmp, short int a1, short int a2) {
	for (short int n = 0; n < a1; n++) {
		free(tmp[n]);
	}
	free(tmp);
}

//========================================================================================================================
//calculate the inverse of a matrix (the code is from http://www.sourcecodesworld.com/source/show.asp?ScriptID=1086)
//========================================================================================================================
void inverse(float** A, short int K, float** I) {
	short int i, j, k;

	for (i = 0; i < K; i++) {
		for (j = 0; j < K; j++) {
			if (i == j)
				I[i][j] = 1;
			else
				I[i][j] = 0;
		}
	}

	/*---------------LoGiC starts here------------------*/		//procedure to make the matrix A to unit matrix
	for (k = 0; k < K; k++) {
		float temp = A[k][k];
		for (j = 0; j < K; j++) {
			A[k][j] /= temp;
			I[k][j] /= temp;
		}
		for (i = 0; i < K; i++) {
			temp = A[i][k];
			for (j = 0; j < K; j++) {
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
void quadprog(float** Dmat, float* dvec, float* xvec, short int K) { //K is the dimension of xvec
	if (K == 1) {
		xvec[0] = 1;
	}
	else {
		float** Dmat_inv, * x_star, s1, s2;
		float Dmat_inv_sum = 0, x_star_sum = 0, * Dmat_inv_rowsum;
		short int num_negatives = 0;
		float* ind_negative;
		Dmat_inv = make2Darray(K, K);
		x_star = (float*)malloc(K * sizeof(float));
		Dmat_inv_rowsum = (float*)malloc(K * sizeof(float));
		ind_negative = (float*)malloc(K * sizeof(float));
		inverse(Dmat, K, Dmat_inv);
		short int k, k1, k2;

		for ( k = 0; k < K; k++) {
			s1 = 0;
			s2 = 0;
			for ( k1 = 0; k1 < K; k1++) {
				s1 += Dmat_inv[k][k1] * dvec[k1];
				Dmat_inv_sum += Dmat_inv[k][k1];
				s2 += Dmat_inv[k][k1];
			}
			x_star[k] = s1;
			Dmat_inv_rowsum[k] = s2;
			x_star_sum += s1;
		}
		for ( k = 0; k < K; k++) {
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
			short int Knew = K - num_negatives, i, j;
			float** Dmat_new, * dvec_new, * xvec_sub;
			Dmat_new = make2Darray(Knew, Knew);
			dvec_new = (float*)malloc((Knew) * sizeof(float));
			xvec_sub = (float*)malloc((Knew) * sizeof(float));
			i = 0;

			for ( k1 = 0; k1 < K; k1++) {
				if (ind_negative[k1] == 0) {
					dvec_new[i] = dvec[k1];
					j = 0;
					for ( k2 = 0; k2 < K; k2++) {
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
			for ( k = 0; k < K; k++) {
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
float val2(float** P_t, float* Eps_M, short int num_celltypes, short int num_CpGsites, short int num_cancer_subtypes,
	float** Ometh, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z, short int n) {
	short int K, M, Q;
	K = num_celltypes;
	M = num_CpGsites;
	Q = num_cancer_subtypes;

	float s1 = 0, s2 = 0, s3 = 0, sum = 0;
	short int k, k2, m, q;

#pragma omp parallel for reduction(+:s1,s2,s3,sum) 
	for (m = 0; m < M; m++) {
		s1 = 0;
		for (k = 0; k < K; k++) {
			for (k2 = 0; k2 < K; k2++) {
				s1 += P_t[k][n] * E_sigma[n][m][k][k2] * P_t[k2][n];
			}
		}

		s2 = 0;
		for (k = 0; k < K; k++) {
			s2 += E_mu[n][m][k] * P_t[k][n];
		}
		s2 = E_z[n][m][q] * pow(Ometh[m][n] - s2, 2);
		sum += (s1 + s2) / (2 * Eps_M[m]);
	}

	return sum;
}

//========================================================================================================================
//functions in the algorithm
//========================================================================================================================
float normal_density_log(float x, float mean, float sd) {
	return -0.5 * log(2 * M_PI) - log(sd) - pow(x - mean, 2) / (2 * sd * sd);
}

float normal_density(float x, float mean, float var) {
	return exp(-pow(x - mean, 2) / (2 * var)) / sqrt(2 * M_PI * var);
}

// update P_t
void Update_P(float** P_t, float* Eps_M, short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes,
	float** Ometh, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {

	short int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	float s1, s2;
	short int k, k2, n, m, q;

	float tar2a, tar2b;
	float* xvec = (float*)malloc(K * sizeof(float));
	float* xvec0 = (float*)malloc(K * sizeof(float));
	float* dvec = (float*)malloc(K * sizeof(float));
	float** Dmat = make2Darray(K, K);
	float** Dmat_inv = make2Darray(K, K);

	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			for (k2 = k; k2 < K; k2++) {
				s1 = 0;
				s2 = 0;
				for (m = 0; m < M; m++) {
					s1 += (E_sigma[n][m][k][k2] + E_mu[n][m][k] * E_mu[n][m][k2]) / Eps_M[m];
					s2 += Ometh[m][n] * E_mu[n][m][k] / Eps_M[m];
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

		tar2a = val2(P_t, Eps_M, K, M, Q, Ometh, E_sigma, E_mu, E_Z, E_z, n);
		quadprog(Dmat, dvec, xvec, K);

		for (k = 0; k < K; k++) {
			P_t[k][n] = xvec[k];
		}

		tar2b = val2(P_t, Eps_M, K, M, Q, Ometh, E_sigma, E_mu, E_Z, E_z, n);

		if (tar2b > tar2a) {//avoid offset effects when estimates are stable
			for (k = 0; k < K; k++) {
				P_t[k][n] = xvec0[k];
			}
		}
	}

	free(xvec); free(xvec0); free(dvec);
	delet2Darray(Dmat, K, K);
	delet2Darray(Dmat_inv, K, K);
}

void Zstep(short int n, short int num_cancer_subtypes, float** E_Z, float* tmp) {

	short int Q;
	Q = num_cancer_subtypes;
	float s1, s2 = 0;
	short int q;

	s1 = tmp[0];
	for ( q = 1; q < Q; q++)
		s1 = MAX(s1, tmp[q]);

	for ( q = 0; q < Q; q++) {
		s2 += exp(tmp[q] - s1);
	}

	for ( q = 0; q < Q; q++) {
		E_Z[n][q] = exp(tmp[q] - s1) / s2;
	}
}

//conduct E step
void EstepI(short int* Labels, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M,
	float** Eps_MK, short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes,
	short int num_covariates, float** Ometh, float** Covariates, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {

	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	float s1 = 0, s2 = 0;
	short int n, m, k, k1, k2, q, p;

#pragma omp parallel for reduction(+:s1,s2) schedule(static,1)
	for ( n = 0; n < N; n++) {
		Labels[n] = 0;
		for (q = 1; q < Q; q++) {
			Labels[n] = AMAX(E_Z[n][Labels[n]], Labels[n], E_Z[n][q], q);
		}
		
		for ( m = 0; m < M; m++) {
			// update E_sigma
			s1 = 0;
			for ( k = 0; k < K; k++) {
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			for ( k1 = 0; k1 < K; k1++) {
				for ( k2 = k1; k2 < K; k2++) {
					E_sigma[n][m][k1][k2] = -1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k2][n] * Eps_MK[m][k1] * Eps_MK[m][k2] / Eps_M[m];
					E_sigma[n][m][k2][k1] = E_sigma[n][m][k1][k2];
				}
				E_sigma[n][m][k1][k1] = Eps_MK[m][k1] - 1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k1][n] * Eps_MK[m][k1] * Eps_MK[m][k1] / Eps_M[m];
			}

			// update E_mu
			for (k = 0; k < K; k++) {
				s1 = 0;
				for (k1 = 0; k1 < K; k1++) {
					s2 = 0;
					for (p = 0; p < P; p++)
						s2 += Beta[m][k1][p] * Covariates[n][p];
					s1 += (Ometh[m][n] * P_t[k1][n] / Eps_M[m] + (Mu[m][k1] + s2 + Gamma[m][k1][Labels[n]]) / Eps_MK[m][k1]) * E_sigma[n][m][k1][k];
				}
				E_mu[n][m][k] = s1;
			}
		}
	}
}

void EstepII(short int* Labels, float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M,
	float** Eps_MK, int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes,
	short int num_covariates, float** Ometh, float** Covariates, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {

	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	float s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	float* tmp = (float*)malloc(Q * sizeof(float));
	float* tp2 = (float*)malloc(Q * sizeof(float));
	short int n, m, k, k1, k2, q, p;

	// #pragma omp parallel for reduction(+:s1,s2) schedule(static,1)
	for ( n = 0; n < N; n++) {
		for (q = 0; q < Q; q++) {
			tmp[q] = Pi_t[q];
		}
		Labels[n] = 0;
		for (q = 1; q < Q; q++) {
			Labels[n] = AMAX(E_Z[n][Labels[n]], Labels[n], E_Z[n][q], q);
		}

		for ( m = 0; m < M; m++) {
			// update E_sigma
			s1 = 0;
			for (k = 0; k < K; k++) {
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			for (k1 = 0; k1 < K; k1++) {
				for (k2 = k1; k2 < K; k2++) {
					E_sigma[n][m][k1][k2] = -1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k2][n] * Eps_MK[m][k1] * Eps_MK[m][k2] / Eps_M[m];
					E_sigma[n][m][k2][k1] = E_sigma[n][m][k1][k2];
				}
				E_sigma[n][m][k1][k1] = Eps_MK[m][k1] - 1.0 / (1.0 + s1 / Eps_M[m]) * P_t[k1][n] * P_t[k1][n] * Eps_MK[m][k1] * Eps_MK[m][k1] / Eps_M[m];
			}

			// update E_mu
			for (k = 0; k < K; k++) {
				s1 = 0;
				for (k1 = 0; k1 < K; k1++) {
					s2 = 0;
					for (p = 0; p < P; p++)
						s2 += Beta[m][k1][p] * Covariates[n][p];
					s1 += (Ometh[m][n] * P_t[k1][n] / Eps_M[m] + (Mu[m][k1] + s2 + Gamma[m][k1][Labels[n]]) / Eps_MK[m][k1]) * E_sigma[n][m][k1][k];
				}
				E_mu[n][m][k] = s1;
			}

			// update E_Z
			s1 = Eps_M[m];
			for ( k = 0; k < K; k++)
				s1 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];

			s4 = 0;
			for ( q = 0; q < Q; q++) {
				s2 = 0;
				for ( k = 0; k < K; k++) {
					s3 = 0;
					for ( p = 0; p < P; p++)
						s3 += Beta[m][k][p] * Covariates[n][p];
					s2 += (Mu[m][k] + s3 + Gamma[m][k][q]) * P_t[k][n];
				}
				tmp[q] += normal_density_log(Ometh[m][n], s2, sqrt(s1));
				tp2[q] = Pi_t[q] * normal_density(Ometh[m][n], s2, s1);
				s4 += tp2[q];
			}
			for (q = 0; q < Q; q++) {
				E_z[n][m][q] = tp2[q] / s4;
			}
		}
		Zstep(n, num_cancer_subtypes, E_Z, tmp);
	}
	free(tmp);
	free(tp2);
}


void MstepI(float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M, float** Eps_MK,
	short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes, short int num_covariates,
	float** Ometh, float** Covariates, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {
	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	float s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	short int n, m, k, k1, q, p, p2, q2;

//	//==============================
//	// update Pi_t
//	//==============================
//#pragma omp parallel for reduction(+:s1)
//	for (int q = 0; q < Q; q++) {
//		s1 = 0;
//		for (int n = 0; n < N; n++) {
//			s1 += E_Z[n][q];
//		}
//		Pi_t[q] = s1 / N;
//	}

	//==============================
	// coordinate descent algorithm
	//==============================
#pragma omp parallel for reduction(+:s1,s2,s3,s4)
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][k] - s2 - Gamma[m][k][q]) * E_z[n][m][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (n = 0; n < N; n++) {
				s1 = 0;
				for (p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_z[n][m][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (n = 0; n < N; n++) {
			s1 = 0;
			for (k = 0; k < K; k++) {
				for (k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (q = 0; q < Q; q++) {
				s2 = 0;
				for (k = 0; k < K; k++)
					s2 += E_mu[n][m][k] * P_t[k][n];
				s3 += (s1 + pow(Ometh[m][n] - s2, 2)) * E_Z[n][q];
			}
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, num_cancer_subtypes, Ometh, E_sigma, E_mu, E_Z, E_z);
}

void MstepII(float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M, float** Eps_MK,
	short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes, short int num_covariates,
	float** Ometh, float** Covariates, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {
	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	float s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	short int n, m, k, k1, q, p, p2, q2;

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
#pragma omp parallel for reduction(+:s1,s2,s3,s4)
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][k] - s2 - Gamma[m][k][q]) * E_z[n][m][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			// update Gamma
			//==============================
			for (q = 1; q < Q; q++) {
				s1 = 0;
				s3 = 0;
				for (n = 0; n < N; n++) {
					s2 = 0;
					s4 = 0;
					for (p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					//for (q2 = 1; q2 < Q; q2++) {
					//	s4 += Gamma[m][k][q2];
					//}
					//s4 -= Gamma[m][k][q];
					s1 += (E_mu[n][m][k] - Mu[m][k] - s2 - s4) * E_z[n][m][q];
					s3 += E_z[n][m][q];
				}
				Gamma[m][k][q] = s1 / s3;
			}

			//==============================
			// update Beta
			//==============================
			for (p = 0; p < P; p++) {
				s1 = 0;
				s3 = 0;
				for (n = 0; n < N; n++) {
					s2 = 0;
					for (p2 = 0; p2 < P; p2++) {
						s2 += Beta[m][k][p2] * Covariates[n][p2];
					}
					s2 -= Beta[m][k][p] * Covariates[n][p];
					for (q = 0; q < Q; q++) {
						s1 += (E_mu[n][m][k] - Mu[m][k] - s2 - Gamma[m][k][q]) * Covariates[n][p] * E_z[n][m][q];
					}
					s3 += (Covariates[n][p] * Covariates[n][p]);
				}
				Beta[m][k][p] = s1 / s3;
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (n = 0; n < N; n++) {
				s1 = 0;
				for (p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_z[n][m][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (n = 0; n < N; n++) {
			s1 = 0;
			for (k = 0; k < K; k++) {
				for (k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (q = 0; q < Q; q++) {
				s2 = 0;
				for (k = 0; k < K; k++)
					s2 += E_mu[n][m][k] * P_t[k][n];
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

void MstepIII(float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M, float** Eps_MK,
	short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes, short int num_covariates,
	float** Ometh, float** Covariates, float**** E_sigma, float*** E_mu, float** E_Z, float*** E_z) {
	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	float s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	short int n, m, k, k1, q, p, p2, q2;

	//==============================
	// coordinate descent algorithm
	//==============================
#pragma omp parallel for reduction(+:s1,s2,s3,s4)
	for (m = 0; m < M; m++) {
		for (k = 0; k < K; k++) {
			//==============================
			// update Mu
			//==============================
			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s2 = 0;
					for (int p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					s1 += (E_mu[n][m][k] - s2 - Gamma[m][k][q]) * E_z[n][m][q];
				}
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			//==============================
			// update Gamma
			//==============================
			for (q = 1; q < Q; q++) {
				s1 = 0;
				s3 = 0;
				for (n = 0; n < N; n++) {
					s2 = 0;
					s4 = 0;
					for (p = 0; p < P; p++)
						s2 += Beta[m][k][p] * Covariates[n][p];
					//for (q2 = 1; q2 < Q; q2++) {
					//	s4 += Gamma[m][k][q2];
					//}
					//s4 -= Gamma[m][k][q];
					s1 += (E_mu[n][m][k] - Mu[m][k] - s2 - s4) * E_z[n][m][q];
					s3 += E_z[n][m][q];
				}
				Gamma[m][k][q] = s1 / s3;
			}

			//==============================
			// update Beta
			//==============================
			for (p = 0; p < P; p++) {
				s1 = 0;
				s3 = 0;
				for (n = 0; n < N; n++) {
					s2 = 0;
					for (p2 = 0; p2 < P; p2++) {
						s2 += Beta[m][k][p2] * Covariates[n][p2];
					}
					s2 -= Beta[m][k][p] * Covariates[n][p];
					for (q = 0; q < Q; q++) {
						s1 += (E_mu[n][m][k] - Mu[m][k] - s2 - Gamma[m][k][q]) * Covariates[n][p] * E_z[n][m][q];
					}
					s3 += (Covariates[n][p] * Covariates[n][p]);
				}
				Beta[m][k][p] = s1 / s3;
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (n = 0; n < N; n++) {
				s1 = 0;
				for (p = 0; p < P; p++)
					s1 += Beta[m][k][p] * Covariates[n][p];
				for (q = 0; q < Q; q++) {
					s3 += (pow(E_mu[n][m][k] - Mu[m][k] - s1 - Gamma[m][k][q], 2) + E_sigma[n][m][k][k]) * E_z[n][m][q];
				}
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (n = 0; n < N; n++) {
			s1 = 0;
			for (k = 0; k < K; k++) {
				for (k1 = 0; k1 < K; k1++)
					s1 += P_t[k][n] * E_sigma[n][m][k][k1] * P_t[k1][n];
			}

			for (q = 0; q < Q; q++) {
				s2 = 0;
				for (k = 0; k < K; k++)
					s2 += E_mu[n][m][k] * P_t[k][n];
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
double observed_log_likelihood(float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M,
	float** Eps_MK, short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes,
	short int num_covariates, float** Ometh, float** Covariates) {

	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;
	double sum = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
	double ploglike = 0;
	short int n, m, k, q, p;

#pragma omp parallel for reduction(+:s1,s2,s3,s4,sum,ploglike)
	for ( m = 0; m < M; m++) {
		for ( n = 0; n < N; n++) {
			s1 = 0;
			s2 = 0;
			s3 = 0;
			sum = 0;

			for ( k = 0; k < K; k++) {
				s1 += P_t[k][n] * Mu[m][k];
				s3 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
				for ( p = 0; p < P; p++)
					s2 += P_t[k][n] * Beta[m][k][p] * Covariates[n][p];
			}

			for ( q = 0; q < Q; q++) {
				s4 = 0;
				for ( k = 0; k < K; k++) {
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
double BIC(double ploglike, short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes,
	short int num_covariates) {
	short int K, N, M, Q, P;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	return (-2.0) * ploglike + log(N) * (Q - 1 + (K - 1) * N + M * (1 + 2 * K + (Q - 1) * K + P * K));
}

void ZInitial(float** E_Z, float*** E_z, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes) {
	short int N, M, Q, n, m, q;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;

	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++) {
			for (q = 0; q < Q; q++) {
				E_z[n][m][q] = E_Z[n][q];
			}
		}
	}
}

//carry out the EM algorithm
void clusterHIREAlgo(float* Pi_t, float** P_t, float** Mu, float*** Beta, float*** Gamma, float* Eps_M, float** Eps_MK,
	short int num_celltypes, short int num_samples, short int num_CpGsites, short int num_cancer_subtypes, short int num_covariates, float** Ometh,
	float** Covariates, float** E_Z, double* BICpointer, float tol, short int num_iteration) {

	short int K, N, M, Q, P, iter = 0, i;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_cancer_subtypes;
	P = num_covariates;

	
	short int* Labels = (short int*)malloc(N * sizeof(short int));
	float**** E_sigma = make4Darray(N, M, K, K);
	float*** E_mu = make3Darray(N, M, K);
	float*** E_z = make3Darray(N, M, Q);
	double ploglike = -100000, ploglike1 = 0;

	ZInitial(E_Z, E_z, N, M, Q);

	while (iter < num_iteration * .1) {
		EstepI(Labels, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		MstepI(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		//ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
		//	K, N, M, Q, P, Ometh, Covariates);
		iter++;
		//Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
		Rprintf("Iteration: %d\n", iter);
	}

	while (iter < num_iteration * .3) {
		EstepII(Labels, Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		MstepII(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		//ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
		//	K, N, M, Q, P, Ometh, Covariates);
		iter++;
		//Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
		Rprintf("Disease subtype: ");
		for (i = 0; i < 10; i++) {
			Rprintf("%d ", Labels[i]);
		}
		Rprintf(", with iteration and gamma: %d %f\n", iter, Gamma[0][0][1]);
	}

	while (iter < num_iteration) {
		EstepI(Labels, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		MstepIII(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK, K, N, M, Q, P,
			Ometh, Covariates, E_sigma, E_mu, E_Z, E_z);
		//ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Beta, Gamma, Eps_M, Eps_MK,
		//	K, N, M, Q, P, Ometh, Covariates);
		iter++;
		//Rprintf("Iteration: %d with BIC: %f\n", iter, ploglike);
		Rprintf("Iteration: %d\n", iter);
	}

	//Rprintf("Iteration: %d\n", iter);
	//while (iter < num_iteration && (ABS(ploglike, ploglike1) / ABS(ploglike1, 0) > tol))
	// double ploglike = observed_log_likelihood(Pi_t, P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh);
	*BICpointer = ploglike;
	delet3Darray(E_mu, N, M, K);
	delet3Darray(E_z, N, M, Q);
	delet4Darray(E_sigma, N, M, K, K);
}