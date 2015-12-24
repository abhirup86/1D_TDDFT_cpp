#include "stdafx.h"

#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

//Use lapack for eigen value/vector solving
int EigenSolverLAPACK(float *matrix, const int &nrow, float *evalue) {
	char jobz = 'V', uplo = 'U';
	int n = nrow, lda = nrow, info = 0;
	int lwork = 2*n*n + 6*n + 1, liwork = 5*n + 3, iwork = liwork;
	float *work = new float[lwork]();
	info = LAPACKE_ssyevd(LAPACK_COL_MAJOR, jobz, uplo, n, matrix, lda, evalue);
	//ssyevd(&jobz, &uplo, &n, matrix, &lda, evalue, work, &lwork, &iwork, &liwork, &info);
	delete [] work;
	if(info == 0) return 0;
	else return -1;
}

int EigenSolverLAPACK(double *matrix, const int &nrow, double *evalue) {
	char jobz = 'V', uplo = 'U';
	int n = nrow, lda = nrow, info = 0;
	int lwork = 2*n*n + 6*n + 1, liwork = 5*n + 3, iwork = liwork;
	double *work = new double[lwork]();
	info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, jobz, uplo, n, matrix, lda, evalue);
	//dsyevd(&jobz, &uplo, &n, matrix, &lda, evalue, work, &lwork, &iwork, &liwork, &info);
	delete [] work;
	if(info == 0) return 0;
	else return -1;
}

int EigenSolverLAPACK(complex<float> *matrix, const int &nrow, float *evalue) {
	char jobz = 'V', uplo = 'U';
	int n = nrow, lda = nrow, info = 0;
	int lwork = 2*n*n + 6*n + 1, liwork = 5*n + 3, iwork = liwork;
	int lrwork = 2*n*n + 65*n + 1;
	complex<float> *work = new complex<float>[lwork]();
	float *rwork = new float[lrwork]();
	cheevd(&jobz, &uplo, &n, matrix, &lda, evalue, work, &lwork, rwork, &lrwork, &iwork, &liwork, &info);
	delete [] work;
	delete [] rwork;
	if(info == 0) return 0;
	else return -1;
}

int EigenSolverLAPACK(complex<double> *matrix, const int &nrow, double *evalue) {
	char jobz = 'V', uplo = 'U';
	int n = nrow, lda = nrow, info = 0;
	int lwork = 2*n*n + 6*n + 1, liwork = 5*n + 3, iwork = liwork;
	int lrwork = 2*n*n + 65*n + 1;
	complex<double> *work = new complex<double>[lwork]();
	double *rwork = new double[lrwork]();
	zheevd(&jobz, &uplo, &n, matrix, &lda, evalue, work, &lwork, rwork, &lrwork, &iwork, &liwork, &info);
	delete [] work;
	delete [] rwork;
	if(info == 0) return 0;
	else return -1;
}

int LinearSolverLAPACK(float *matrix, const int &nrow, float *vec) {
	int n = nrow, nrhs = 1, lda = nrow, ldb = nrow, info = 0;
	int *ipiv = new int[n]();
	info = LAPACKE_sgesv(LAPACK_COL_MAJOR, n, nrhs, matrix, lda, ipiv, vec, ldb);
	delete [] ipiv;
	if(info == 0) return 0;
	else return -1;
}

int LinearSolverLAPACK(double *matrix, const int &nrow, double *vec) {
	int n = nrow, nrhs = 1, lda = nrow, ldb = nrow, info = 0;
	int *ipiv = new int[n]();
	info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, matrix, lda, ipiv, vec, ldb);
	delete [] ipiv;
	if(info == 0) return 0;
	else return -1;
}

int LinearSolverLAPACK(complex<float> *matrix, const int &nrow, complex<float> *vec) {
	int n = nrow, nrhs = 1, lda = nrow, ldb = nrow, info = 0;
	int *ipiv = new int[n]();
	info = LAPACKE_cgesv(LAPACK_COL_MAJOR, n, nrhs, matrix, lda, ipiv, vec, ldb);
	delete [] ipiv;
	if(info == 0) return 0;
	else return -1;
}

int LinearSolverLAPACK(complex<double> *matrix, const int &nrow, complex<double> *vec) {
	int n = nrow, nrhs = 1, lda = nrow, ldb = nrow, info = 0;
	int *ipiv = new int[n]();
	info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, matrix, lda, ipiv, vec, ldb);
	delete [] ipiv;
	if(info == 0) return 0;
	else return -1;
}

// result = matrix*vec
void LinearAlgebraMul(const double *matrix, const double *vec, const int &dim, double *result) {
	for(int i = 0; i < dim; ++i) result[i] = 0.0;
	for(int i = 0; i < dim; ++i)
		for(int j = 0; j < dim; ++j)
			result[i] += matrix[i*dim+j]*vec[j];
	//dzgemv can be used as well
}

// result = matrix*vec
void LinearAlgebraMul(const double *matrix, const complex<double> *vec, const int &dim, complex<double> *result) {
	for(int i = 0; i < dim; ++i) result[i] = complex<double>(0.0, 0.0);
	for(int i = 0; i < dim; ++i)
		for(int j = 0; j < dim; ++j)
			result[i] += matrix[i*dim+j]*vec[j];
	//dzgemv can be used as well
}

// vec = matrix*vec (in place)
void LinearAlgebraMul(const double *matrix, double *vec, const int &dim) {
	double *tmp = new double[dim]();
	for(int i = 0; i < dim; ++i)
		for(int j = 0; j < dim; ++j)
			tmp[i] += matrix[i*dim+j]*vec[j];
	// switch pointers.
	for(int i = 0; i < dim; ++i) vec[i] = tmp[i];
	delete [] tmp;
	//dzgemv can be used as well
}

// vec = matrix*vec (in place)
void LinearAlgebraMul(const double *matrix, complex<double> *vec, const int &dim) {
	complex<double> *tmp = new complex<double>[dim]();
	for(int i = 0; i < dim; ++i)
		for(int j = 0; j < dim; ++j)
			tmp[i] += matrix[i*dim+j]*vec[j];
	// switch pointers.
	for(int i = 0; i < dim; ++i) vec[i] = tmp[i];
	delete [] tmp;
	//dzgemv can be used as well
}

// vec = matrix*vec (in place)
void LinearAlgebraMul(const complex<double> *matrix, complex<double> *vec, const int &dim) {
	complex<double> *tmp = new complex<double>[dim]();
	for(int i = 0; i < dim; ++i)
		for(int j = 0; j < dim; ++j)
			tmp[i] += matrix[i*dim+j]*vec[j];
	// switch pointers.
	for(int i = 0; i < dim; ++i) vec[i] = tmp[i];
	delete [] tmp;
	//dzgemv can be used as well
}

// result = alpha*vec
void LinearAlgebraMul(const double alpha, const double *vec, const int &dim, double *result) {
	for(int i = 0; i < dim; ++i) result[i] = 0.0;
	for(int i = 0; i < dim; ++i) result[i] = alpha*vec[i];
}

// result = alpha*vec
void LinearAlgebraMul(const complex<double> alpha, const complex<double> *vec, const int &dim, complex<double> *result) {
	for(int i = 0; i < dim; ++i) result[i] = complex<double>(0.0, 0.0);
	for(int i = 0; i < dim; ++i) result[i] = alpha*vec[i];
}

// vec = alpha*vec(in place)
void LinearAlgebraMul(const double alpha, double *vec, const int &dim) {
	double *tmp = new double[dim]();
	for(int i = 0; i < dim; ++i) tmp[i] = alpha*vec[i];
	delete [] vec;
	vec = tmp;
}

// vec = alpha*vec(in place)
void LinearAlgebraMul(const complex<double> alpha, complex<double> *vec, const int &dim) {
	complex<double> *tmp = new complex<double>[dim]();
	for(int i = 0; i < dim; ++i) tmp[i] = alpha*vec[i];
	delete [] vec;
	vec = tmp;
}

// vec = a*vec
//void LinearAlgebraMul(const complex<double> alpha, const complex<double> *vec, const int &dim, complex<double> *result) {
//	//add alpha*vec to result
//	cblas_zaxpy (dim, &alpha, vec, 1, result, 1);
//}
//

#endif //LINEARALGEBRA_H
