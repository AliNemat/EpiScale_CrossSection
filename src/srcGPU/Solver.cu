#include "Solver.h"

using namespace std ; 

/*
vector<double>  Solver::solve3Diag(const vector <double> & lDiag, const vector <double> & diag, const vector <double> & uDiag,
				                   const vector <double> & rHS) {

   // --- Initialize cuSPARSE
    cusparseHandle_t handle;    cusparseCreate(&handle);

    const int N     =diag.size();        // --- Size of the linear system

    // --- Lower diagonal, diagonal and upper diagonal of the system matrix
    double *h_ld = (double*)malloc(N * sizeof(double));
    double *h_d  = (double*)malloc(N * sizeof(double));
    double *h_ud = (double*)malloc(N * sizeof(double));
    double *h_x = (double *)malloc(N * sizeof(double)); 

    for (int k = 0; k < N ; k++) {
       h_ld[k] =lDiag.at(k);
	   h_d[k]  =diag.at(k);
       h_ud[k] =uDiag.at(k);
	   h_x[k]  =rHS.at(k); 
    }
   // for (int k = 0; k < N; k++) 
    
	double *d_ld;   cudaMalloc(&d_ld, N * sizeof(double));
    double *d_d;    cudaMalloc(&d_d,  N * sizeof(double));
    double *d_ud;   cudaMalloc(&d_ud, N * sizeof(double));
    double *d_x;       cudaMalloc(&d_x, N * sizeof(double));   
    cout << "lower diagonal elements" << endl ; 
    for (int k=0; k<N; k++) printf("%f\n", h_ld[k]);
    cout << "diagonal elements" << endl ; 
    for (int k=0; k<N; k++) {
		cout << k <<", " <<h_d[k] << endl ;  
	}
    cout << "upper diagonal elements" << endl ; 
    for (int k=0; k<N; k++) printf("%f\n", h_ud[k]);
    cout << "RHS elements" << endl ; 
    for (int k=0; k<N; k++) printf("%f\n", h_x[k]);
    cudaMemcpy(d_ld, h_ld, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d,  h_d,  N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ud, h_ud, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, h_x, N * sizeof(double), cudaMemcpyHostToDevice);

    // --- Allocating and defining dense host and device data vectors
 //   h_x[0] = 100.0;  h_x[1] = 200.0; h_x[2] = 400.0; h_x[3] = 500.0; h_x[4] = 300.0;

    // --- Allocating the host and device side result vector
    //double *h_y = (double *)malloc(N * sizeof(double)); 
    //double *d_y;        cudaMalloc(&d_y, N * sizeof(double)); 
    cusparseDgtsv(handle, N, 1, d_ld, d_d, d_ud, d_x, N);

    cudaMemcpy(h_x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);
    cout << "results" << endl ; 
    for (int k=0; k<N; k++) printf("%f\n", h_x[k]);
    cout << "finished results" << endl ; 
	vector < double> ans ; 
	for (int k=0; k<N; k++) {
	   ans.push_back(h_x[k]); 

	}

	return ans ; 

}
 */
 vector<double>  Solver::SOR3DiagPeriodic(const vector <bool> & nodeIsActive, const vector <double> & lDiag, const vector <double> & diag, const vector <double> & uDiag,
	                                      const vector <double> & rHS, 
										  const vector <int> & prevIndex,
										  const vector <int> & nextIndex,
										  vector<double>  & firstGuess) {
				
	const int maxIteration=2500 ; 
	const double beta=1.2 ; 
	
	vector <double> ans =firstGuess ;
	vector<double> ansOld=firstGuess ; 
	double maxError=10000 ;  // Just a big number to go inside while loop for the fist time. 
    const int N     =diag.size();        // --- Size of the linear system
	
	//Simple iterative SOR solver
	int numIterator=0 ; 
	while (maxError>1.0E-6 && numIterator<maxIteration) {
		maxError=0 ; 
		numIterator ++ ; 
		for (int i=0; i<N ; i++) {
			if (!nodeIsActive.at(i)) {
				continue ; 
			}
			ansOld.at(i)=ans.at(i) ; 
			ans.at(i)=beta*(rHS.at(i)-lDiag.at(i)*ans[prevIndex.at(i)]-uDiag.at(i)*ans[nextIndex.at(i)])/diag.at(i)+ (1-beta)*ansOld.at(i); 
			if ( abs (ansOld.at(i)-ans.at(i)) >maxError) {
				maxError=abs (ansOld.at(i)-ans.at(i));  
			}
		}
	}
	cout << "In SOR solver after " << numIterator <<" iteration, maximum difference in two successuve iterations is "<< maxError << endl ;  ; 
	return ans; 
 }
