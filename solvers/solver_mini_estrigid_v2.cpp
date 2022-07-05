#include <Eigen/Dense>
#include "mex.h"




using namespace Eigen;




MatrixXcd solver_mini_estrigid_v2(const VectorXd& data)
{
	// Compute coefficients
    const double* d = data.data();
    VectorXd coeffs(9);
    coeffs[0] = d[0];
    coeffs[1] = d[1];
    coeffs[2] = d[2];
    coeffs[3] = d[3];
    coeffs[4] = d[4];
    coeffs[5] = d[5];
    coeffs[6] = d[6];
    coeffs[7] = 1;
    coeffs[8] = -1;



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,7,1,0,1,7,7,2,0,7,3,2,1,7,4,2,7,8,3,7,7 };
	static const int coeffs1_ind[] = { 6,8,6,4,8,5,4,3,6,5,8,5,8,7,7 };
	static const int C0_ind[] = { 0,7,8,9,17,18,23,24,27,30,32,33,35,36,40,43,45,47,49,54,58 } ;
	static const int C1_ind[] = { 3,5,8,11,14,16,17,19,25,27,28,33,34,37,44 };

	Matrix<double,8,8> C0; C0.setZero();
	Matrix<double,8,6> C1; C1.setZero();
	for (int i = 0; i < 21; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 15; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,8,6> C12 = C0.partialPivLu().solve(C1);



	// Setup action matrix
	// Matrix<double,8, 6> RR;
    MatrixXd RR(8, 6);	
	RR << -C12.bottomRows(2), Matrix<double,6,6>::Identity(6, 6);

	static const int AM_ind[] = { 5,4,0,6,7,1 };
	// Matrix<double, 6, 6> AM;
    MatrixXd AM(6, 6);
	for (int i = 0; i < 6; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 2, 6> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 6, 6> > es(AM);
	ArrayXcd D = es.eigenvalues();	
	ArrayXXcd V = es.eigenvectors();

V = (V / V.row(0).array().replicate(6, 1)).eval();


    sols.row(0) = V.row(1).array();
    sols.row(1) = D.transpose().array();









	return sols;
}
// Action =  y
// Quotient ring basis (V) = 1,x,x*y,y,y^2,y^3,
// Available monomials (RR*V) = x*y^2,y^4,1,x,x*y,y,y^2,y^3,





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:mini_estrigid_v2:nrhs", "One input required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:mini_estrigid_v2:nlhs", "One output required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:mini_estrigid_v2:notDouble", "Input data must be type double.");
	}
	if(mxGetNumberOfElements(prhs[0]) % 7 != 0) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:mini_estrigid_v2:incorrectSize", "Input size must be multiple of 7.");
	}
	int n_instances = mxGetNumberOfElements(prhs[0]) / 7;
	double *input = mxGetPr(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix(2,6*n_instances,mxCOMPLEX);
	double* zr = mxGetPr(plhs[0]);
	double* zi = mxGetPi(plhs[0]);
	for(int k = 0; k < n_instances; k++) {
		const VectorXd data = Map<const VectorXd>(input + k*7, 7);
		MatrixXcd sols = solver_mini_estrigid_v2(data);
		Index offset = k*sols.size();
		for (Index i = 0; i < sols.size(); i++) {
			zr[i+offset] = sols(i).real();
			zi[i+offset] = sols(i).imag();
		}
	}
}


