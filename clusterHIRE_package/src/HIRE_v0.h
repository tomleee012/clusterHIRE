#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b)) 
#define ABS(a, b)  (((a) > (b)) ? (a-b) : (b-a)) 

//========================================================================================================================
//Define multiple dimensional arrays (Memory is not continuous here)
//========================================================================================================================

double**** make4Darray(int a1, int a2, int a3, int a4){
	double ****tmp;
	tmp = (double ****)malloc(a1*sizeof(double ***));
	for(int i=0; i<a1; i++){
		tmp[i] = (double ***) malloc(a2*sizeof(double **));
		for(int j=0; j<a2; j++){
			tmp[i][j] = (double **) malloc(a3*sizeof(double*));
			for(int k=0; k<a3; k++){
				tmp[i][j][k] = (double *)malloc(a4*sizeof(double));
			}
		}
	}
	return tmp;
}

void delet4Darray(double ****tmp, int a1, int a2, int a3, int a4){
	for(int i=0; i<a1; i++){
		for(int j=0; j<a2; j++){
			for(int k=0; k<a3; k++){
				free(tmp[i][j][k]);
			}
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double*** make3Darray(int a1, int a2, int a3){
	double ***tmp;
	tmp = (double ***)malloc(a1*sizeof(double **));
	for(int i=0; i<a1; i++){
		tmp[i] = (double **) malloc(a2*sizeof(double *));
		for(int j=0; j<a2; j++){
			tmp[i][j] = (double *) malloc(a3*sizeof(double));
		}
	}
	return tmp;
}

void delet3Darray(double ***tmp, int a1, int a2, int a3){
	for(int i=0; i<a1; i++){
		for(int j=0; j<a2; j++){
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double** make2Darray(int a1, int a2){
	double **tmp;
	tmp = (double **)malloc(a1*sizeof(double *));
	for(int i=0; i<a1; i++){
		tmp[i] = (double *) malloc(a2*sizeof(double));
	}
	return tmp;
}

void delet2Darray(double **tmp, int a1, int a2){
	for(int i=0; i<a1; i++){
		free(tmp[i]);
	}
	free(tmp);
}


//========================================================================================================================
//calculate the inverse of a matrix (the code is from http://www.cs.rochester.edu/u/brown/Crypto/assts/projects/adj.html)
//========================================================================================================================

/*Recursive definition of determinate using expansion by minors*/
double Determinant(double **a,int n){
   int i,j,j1,j2;
   double det = 0;
   double **m=NULL;

   if (n < 1) { /* Error */
	  
   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (double **) malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = (double *) malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*Find the cofactor matrix of a square matrix*/

void CoFactor(double **a, int n, double **b){
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = (double **) malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = (double *) malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

/*Transpose of a square matrix, do it in place*/
void Transpose(double **a,int n){
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}

/*calculate the inverse*/
void inverse(double **a, int n, double **a_inv){
	double det;
	double  **cofac_a;
	cofac_a = (double **)malloc(n*sizeof(double*));
	for(int i =0; i <n;i++){
		cofac_a[i] = (double *)malloc(n*sizeof(double));
	}
	CoFactor(a, n, cofac_a);
	Transpose(cofac_a, n); //turn the cofacotor matrix into the adjoint matrix
	det = Determinant(a, n); 
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			a_inv[i][j] = cofac_a[i][j] / det;
		}
	}

	for(int i=0; i<n; i++){
		free(cofac_a[i]);
	}
	free(cofac_a);
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double **Dmat, double *dvec, double *xvec, int K){ //K is the dimension of xvec
	if(K==1){
		xvec[0] = 1;
	}else{
		double **Dmat_inv, *x_star, s1, s2;
		double Dmat_inv_sum=0, x_star_sum=0, *Dmat_inv_rowsum;
		int num_negatives=0;
		double *ind_negative;
		Dmat_inv = (double **)malloc(K*sizeof(double*));
		x_star = (double *)malloc(K*sizeof(double));
		Dmat_inv_rowsum = (double *)malloc(K*sizeof(double));
		ind_negative = (double *)malloc(K*sizeof(double));
		for(int k=0; k<K; k++){
			Dmat_inv[k] = (double *)malloc(K*sizeof(double));
		}
		inverse(Dmat, K, Dmat_inv);
		for(int k=0; k<K; k++){
			s1 = 0;
			s2 = 0;
			for(int k1=0; k1<K; k1++){
				s1 += Dmat_inv[k][k1]*dvec[k1];
				Dmat_inv_sum += Dmat_inv[k][k1];
				s2 += Dmat_inv[k][k1];
			}
			x_star[k] = s1;
			Dmat_inv_rowsum[k] = s2;
			x_star_sum += s1;
		}
		for(int k=0; k<K; k++){
			xvec[k] = x_star[k] + (1-x_star_sum)/Dmat_inv_sum*Dmat_inv_rowsum[k];
			if(xvec[k]<0){
				num_negatives++;
				ind_negative[k] = 1; 
			}else{
				ind_negative[k] = 0;
			}
		}
		free(x_star);
		free(Dmat_inv_rowsum);
		for(int k=0; k<K; k++){
			free(Dmat_inv[k]);
		}
		free(Dmat_inv);
		
		if(num_negatives == 0){
			free(ind_negative);
		}else{
			int Knew = K-num_negatives, i, j;
			double ** Dmat_new, *dvec_new, *xvec_sub;
			Dmat_new = (double **)malloc((Knew)*sizeof(double*));
			for(int k=0; k<Knew; k++){
				Dmat_new[k] = (double *) malloc((Knew)*sizeof(double));
			}
			dvec_new = (double *) malloc((Knew)*sizeof(double));
			xvec_sub = (double *) malloc((Knew)*sizeof(double));
			i = 0;

			for(int k1=0; k1<K; k1++){
				if(ind_negative[k1]==0){
					dvec_new[i] = dvec[k1];
					j = 0;
					for(int k2=0; k2<K; k2++){
						if(ind_negative[k2]==0){
							Dmat_new[i][j] = Dmat[k1][k2];
							j++;
						}
					}
					i++;
				}else{
					xvec[k1] = 0;
				}
			}

			quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
			i=0;
			for(int k=0; k<K; k++){
				if(ind_negative[k]==0){
					xvec[k] = xvec_sub[i];
					i++;
				}
			}
			free(dvec_new);
			free(xvec_sub);
			for(int k=0; k<Knew;k++){
				free(Dmat_new[k]);
			}
			free(Dmat_new);
		}
		
	}
	
} 


//the objective function when updating P_t
double val2(double **P_t, double *Eps_M,  int num_celltypes,   
					int num_CpGsites,  double **Ometh, double ****E_sigma, double ***E_mu, int n){
	int K,M;
	K = num_celltypes;
	M = num_CpGsites;

	double s1,s2,sum=0;

	for(int m=0; m<M; m++){
			s1=0;
			s2=0;
			for(int k=0; k<K; k++){
				for(int k2=0; k2<K; k2++){
					s1+= P_t[k][n]*E_sigma[n][m][k][k2]*P_t[k2][n]; 
				}
				s2+= E_mu[n][m][k]*P_t[k][n];
			}
			s2 = Ometh[m][n] - s2;
			s2 = s2*s2;
			sum += (s1+s2)/(2*Eps_M[m]);
	}

	return sum;												
}

//========================================================================================================================
//functions in the algorithm
//========================================================================================================================

double normal_density_log(double x, double mean, double sd) {
	return -0.5 * log(2 * M_PI) - log(sd) - pow(x - mean, 2) / (2 * sd * sd);
}

// no use without lasso penalty
double positive_part(double x){
	if(x >= 0){
		return x;
	}else{
		return 0;
	}
}

double absolute(double x){
	if(x >= 0){
		return x;
	}else{
		return -x;
	}
}

double sign(double x){
	if(x > 0){
		return 1.0;
	}else if(x == 0){
		return 0.0;
	}else{
		return -1.0;
	}
}

// update P_t
void Update_P(double** P_t, double* Eps_M, int num_celltypes, int num_samples, int num_CpGsites,
	double** Ometh, double**** E_sigma, double*** E_mu) {

	int K, N, M;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	double s1, s2;

	double* dvec, ** Dmat, ** Dmat_inv, * xvec, * xvec0, tar2a, tar2b;
	xvec = (double*)malloc(K * sizeof(double));
	xvec0 = (double*)malloc(K * sizeof(double));
	dvec = (double*)malloc(K * sizeof(double));
	Dmat = (double**)malloc(K * sizeof(double*));
	Dmat_inv = (double**)malloc(K * sizeof(double*));
	for (int k = 0; k < K; k++) {
		Dmat[k] = (double*)malloc(K * sizeof(double));
		Dmat_inv[k] = (double*)malloc(K * sizeof(double));
	}
	for (int n = 0; n < N; n++) {

		for (int k = 0; k < K; k++) {
			for (int k2 = k; k2 < K; k2++) {
				s1 = 0;
				for (int m = 0; m < M; m++) {
					s1 += (E_sigma[n][m][k][k2] + E_mu[n][m][k] * E_mu[n][m][k2]) / (Eps_M[m]);
				}
				Dmat[k][k2] = s1;
				if (k2 > k) {
					Dmat[k2][k] = Dmat[k][k2];
				}
			}

			s2 = 0;
			for (int m = 0; m < M; m++) {
				s2 += Ometh[m][n] * E_mu[n][m][k] / (Eps_M[m]);
			}
			dvec[k] = s2;
		}

		for (int k = 0; k < K; k++) {
			xvec0[k] = P_t[k][n];
		}

		tar2a = val2(P_t, Eps_M, K, M, Ometh, E_sigma, E_mu, n);

		quadprog(Dmat, dvec, xvec, K);

		for (int k = 0; k < K; k++) {
			P_t[k][n] = xvec[k];
		}
		tar2b = val2(P_t, Eps_M, K, M, Ometh, E_sigma, E_mu, n);

		if (tar2b > tar2a) {//avoid offset effects when estimates are stable
			for (int k = 0; k < K; k++) {
				P_t[k][n] = xvec0[k];
			}
		}
	}
	free(xvec);
	free(xvec0);
	free(dvec);
	for (int k = 0; k < K; k++) {
		free(Dmat[k]);
		free(Dmat_inv[k]);
	}
	free(Dmat);
	free(Dmat_inv);
}

//conduct E step
void Estep(double **P_t, double **Mu, double ***Gamma, double *Eps_M, 
					double **Eps_MK, int num_celltypes, int num_samples, 
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
					double ****E_sigma, double ***E_mu ){
	int K,N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_covariates;
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

			for (int k1 = 0; k1 < K; k1++) {
				s2 = 0;
				for (int k2 = 0; k2 < K; k2++) {
					s1 = 0;
					for (int q = 0; q < Q; q++)
						s1 += Gamma[m][k2][q] * X_scaled[q][n];
					s2 += (Ometh[m][n] * P_t[k2][n] / Eps_M[m] + (Mu[m][k2] + s1) / Eps_MK[m][k2]) * E_sigma[n][m][k2][k1];
				}
				E_mu[n][m][k1] = s2;
			}
		}
	}
					
}

// conduct M step
void Mstep(double **P_t, double **Mu, double ***Gamma, double *Eps_M, 
					double **Eps_MK, int num_celltypes, int num_samples, 
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
					double ****E_sigma, double ***E_mu, int start){
	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_covariates;
	double s1 = 0, s2 = 0, s3 = 0;
	//=================================
	// use coordinate descent
	//=================================

	//use coordinate descent algorithm to update Mu and Gamma
#pragma omp parallel for reduction(+:s1,s2,s3)
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			s1 = 0;
			for (int n = 0; n < N; n++) {
				s2 = 0;
				for (int q = 0; q < Q; q++)
					s2 += Gamma[m][k][q] * X_scaled[q][n];
				s1 += E_mu[n][m][k] - s2;
			}
			Mu[m][k] = MAX(0, MIN(1, s1 / N));

			for (int q = start; q < Q; q++) {
				s1 = 0;
				s3 = 0;
				for (int n = 0; n < N; n++) {
					s2 = 0;
					for (int q2 = 0; q2 < Q; q2++) {
						if (q2 == q) continue;
						s2 += Gamma[m][k][q2] * X_scaled[q2][n];
					}
					s1 += (E_mu[n][m][k] - Mu[m][k] - s2) * X_scaled[q][n];
					s3 += (X_scaled[q][n] * X_scaled[q][n]);
				}
				Gamma[m][k][q] = s1 / s3;
			}

			//==============================
			//update Eps_MK
			//==============================
			s3 = 0;
			for (int n = 0; n < N; n++) {
				s2 = 0;
				for (int q = 0; q < Q; q++) {
					s2 += Gamma[m][k][q] * X_scaled[q][n];
				}
				s3 += pow(E_mu[n][m][k] - Mu[m][k] - s2, 2) + E_sigma[n][m][k][k];
			}
			Eps_MK[m][k] = s3 / N;
		}
		//==============================
		//update Eps_M
		//==============================
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k1 = 0; k1 < K; k1++) {
				for (int k2 = 0; k2 < K; k2++)
					s1 += P_t[k1][n] * E_sigma[n][m][k1][k2] * P_t[k2][n];
			}
			s2 = 0;
			for (int k = 0; k < K; k++)
				s2 += E_mu[n][m][k] * P_t[k][n];
			s3 += (s1 + pow(Ometh[m][n] - s2, 2));
		}
		Eps_M[m] = s3 / N;
	}

	//==============================
	//update P_t
	//==============================
	Update_P(P_t, Eps_M, num_celltypes, num_samples, num_CpGsites, Ometh, E_sigma, E_mu);
}


//calculate the observed data likelihood
double observed_log_likelihood(double** P_t, double** Mu, double*** Gamma, double* Eps_M,
	double** Eps_MK, int num_celltypes, int num_samples,
	int num_CpGsites, int num_covariates, double** Ometh, double** X_scaled) {
	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_covariates;
	double sum = 0, s1 = 0, s2 = 0, s3 = 0;

#pragma omp parallel for reduction(+:s1,s2,s3,sum)
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			s1 = s2 = s3 = 0;
			for (int k = 0; k < K; k++) {
				s1 += P_t[k][n] * Mu[m][k];
			}

			for (int k = 0; k < K; k++) {
				for (int q = 0; q < Q; q++) {
					s2 += P_t[k][n] * Gamma[m][k][q] * X_scaled[q][n];
				}
			}

			for (int k = 0; k < K; k++) {
				s3 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			sum += normal_density_log(Ometh[m][n], s1 + s2, sqrt(s3 + Eps_M[m]));
		}
	}
	return sum;
}

//calculate the BIC to select an appropriate K
double BIC(double** P_t, double** Mu, double*** Gamma, double* Eps_M,
	double** Eps_MK, int num_celltypes, int num_samples,
	int num_CpGsites, int num_covariates, double** Ometh, double** X_scaled) {
	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_covariates;
	double sum = 0, s1, s2, s3;
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			s1 = s2 = s3 = 0;
			for (int k = 0; k < K; k++) {
				s1 += P_t[k][n] * Mu[m][k];
			}

			for (int k = 0; k < K; k++) {
				for (int q = 0; q < Q; q++) {
					s2 += P_t[k][n] * Gamma[m][k][q] * X_scaled[q][n];
				}
			}

			for (int k = 0; k < K; k++) {
				s3 += P_t[k][n] * P_t[k][n] * Eps_MK[m][k];
			}
			sum += normal_density_log(Ometh[m][n], s1 + s2, sqrt(s3 + Eps_M[m]));
		}
	}

	return (-2.0) * sum + log(N) * ((K - 1) * N + M * (1 + 2 * K + Q * K));
}

//carry out the EM algorithm
void EmEwas(double **P_t, double **Mu, double ***Gamma, double *Eps_M, 
					double **Eps_MK, int num_celltypes, int num_samples, 
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled, 
					double *BICpointer, double tol, int num_iteration, int start){ //tol is the tolerance deciding when the algorithm stops 
						
	int K, N, M, Q, iter=0;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_covariates;
			
	double*** E_mu = make3Darray(N, M, K);
	double**** E_sigma = make4Darray(N, M, K, K);

	double ploglike0 = -100000, ploglike1 = 1;
	while ((ABS(ploglike1, ploglike0) / ABS(ploglike0, 0) > tol) && (iter < num_iteration)) { //relative tolerance
		ploglike0 = ploglike1;

		Estep(P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh, X_scaled, E_sigma, E_mu);
		Mstep(P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh, X_scaled, E_sigma, E_mu, start);
		
		iter++;
		ploglike1 = observed_log_likelihood(P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh, X_scaled);
		Rprintf("Iteration: %d\t observed-data log likelihood: %lf\n", iter, ploglike1);

	}

	*BICpointer = BIC(P_t, Mu, Gamma, Eps_M, Eps_MK, K, N, M, Q, Ometh, X_scaled);
	delet3Darray(E_mu, N, M, K);
	delet4Darray(E_sigma, N, M, K, K);
}

