// Диплом.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int nx = 10;
const int ny = 10;
const int unx = nx + 1;
const int uny = ny;
const int vnx = nx;
const int vny = ny + 1;
const double t = 0.01;
const double h = 0.1;
const double Re = 30.;

void iteration(double** u, double** v);
void Pard(double** pMatrix, double* rightPart, int N, double* p);
double convectionDiffusion(double** w, double u, double v, int i, int j);
double** presureSolver(double** u, double** v);
void pBorder(double** pMatrix, double* rightPart, double** u, double** v, int& matrixLineNumber);
void pInnerPart(double** pMatrix, double* rightPart, double** u, double** v, int& matrixLineNumber);
void border(double** u, double** v);
void addRegularConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addLeftMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addRightMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addLeftAndTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addLeftAndBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addRightAndTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void addRightAndBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v);
void printVel(double** u, double** v);

int main()
{
	double** u0 = NULL;
	double** v0 = NULL;

	u0 = new double* [uny];
	for (int i = 0; i < uny; i++)
	{
		u0[i] = new double[unx];
		for (int j = 0; j < unx; j++) {
			u0[i][j] = 0;
		}
	}

	v0 = new double* [vny];
	for (int i = 0; i < vny; i++)
	{
		v0[i] = new double[vnx];
		for (int j = 0; j < vnx; j++) {
			v0[i][j] = 0;
		}
	}

	for (int i = 0; i < 10000; i++) {
		iteration(u0, v0);
	}
	
	printVel(u0, v0);

	ofstream out("Kartina.dat"); //ќткрываем выходной файл
	if (!out)
	{
		cout << "Error";
		return 1;
	}
	for (int i = 0; i < ny; i++) {
		for (int j = 0; j < nx; j++) {
			out << i << " " << j << " " << u0[i][j] << " " << v0[i][j] << "\n";
		}
	}

	out.close();

	return 0;
};

void printVel(double** u, double** v) {
	printf("\n");
	for (int i = 0; i < uny; i++) {
		for (int j = 0; j < unx; j++) {
			printf("u(%d, %d)=%lf\n", i, j, u[i][j]);
		}
	}
	for (int i = 0; i < vny; i++) {
		for (int j = 0; j < vnx; j++) {
			printf("v(%d, %d)=%lf\n", i, j, v[i][j]);
		}
	}
}

void iteration(double** u, double** v) {
	double** u12 = NULL;
	double** v12 = NULL;

	border(u, v);

	u12 = new double* [uny];
	for (int i = 0; i < uny; i++)
	{
		u12[i] = new double[unx];
		for (int j = 0; j < unx; j++) {
			u12[i][j] = 0;
		}
	}

	v12 = new double* [vny];
	for (int i = 0; i < vny; i++)
	{
		v12[i] = new double[vnx];
		for (int j = 0; j < vnx; j++) {
			v12[i][j] = 0;
		}
	}

	for (int i = 1; i < uny - 1; i++) {
		for (int j = 1; j < unx - 1; j++) {
			u12[i][j] = u[i][j] - t * convectionDiffusion(u, u[i][j], (v[i][j] + v[i][j + 1] + v[i + 1][j] + v[i + 1][j + 1])/4., i, j);
		}
	}

	for (int i = 1; i < vny - 1; i++) {
		for (int j = 1; j < vnx - 1; j++) {
			v12[i][j] = v[i][j] - t * convectionDiffusion(v, (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1])/4., v[i][j], i, j);
		}
	}

	border(u12, v12);

	double** p = presureSolver(u12, v12);

	for (int i = 1; i < uny - 1; i++) {
		for (int j = 1; j < unx - 1; j++) {
			u[i][j] = u12[i][j] - (t / h) * (p[i][j] - p[i][j - 1]);
		}
	}

	for (int i = 1; i < vny - 1; i++) {
		for (int j = 1; j < vnx - 1; j++) {
			v[i][j] = v12[i][j] - (t / h) * (p[i][j] - p[i - 1][j]);
		}
	}

	border(u, v);
}

double** presureSolver(double** u, double** v) {
	int matrixLineNumber = 0;
	double* p = new double[nx * ny];
	double** pMatrix = new double*[nx * ny];
	for (int i = 0; i < nx * ny; i++) {
		pMatrix[i] = new double [nx * ny];
	}
	double* rightPart = new double[nx * ny];


	pBorder(pMatrix, rightPart, u, v, matrixLineNumber);

	Pard(pMatrix, rightPart, nx * ny, p);

	double** result = new double* [ny];
	for (int i = 0; i < ny; i++) {
		result[i] = new double[nx];
		for (int j = 0; j < nx; j++) {
			result[i][j] = p[i * nx + j];
		}
	}

	return result;
}

void pBorder(double** pMatrix, double* rightPart, double** u, double** v, int &matrixLineNumber) {
	for (int i = 1; i < nx - 1; i++) {
		addBottomMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, i, 0, u, v);
		addTopMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, i, ny - 1, u, v);
	}

	for (int i = 1; i < ny - 1; i++) {
		addRightMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, nx - 1, i, u, v);
		addLeftMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, 0, i, u, v);
	}

	addLeftAndBottomMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, 0, 0, u, v);
	addLeftAndTopMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, 0, ny - 1, u, v);
	addRightAndBottomMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, nx - 1, 0, u, v);
	addRightAndTopMissConditionPressure(pMatrix, rightPart, matrixLineNumber++, nx - 1, ny - 1, u, v);
}

void pInnerPart(double** pMatrix, double* rightPart, double** u, double** v, int& matrixLineNumber) {
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			addRegularConditionPressure(pMatrix, rightPart, matrixLineNumber++, i, j, u, v);
		}
	}
}

void border(double** u, double** v) {
	for (int i = 0; i < uny; i++) {
		u[i][unx - 1] = 0;
		u[i][0] = 0;
	}

	for (int j = 0; j < unx; j++) {
		u[0][j] = 0;
		u[uny - 1][j] = 1;
	}

	for (int i = 0; i < vny; i++) {
		v[i][vnx - 1] = 0;
		v[i][0] = 0;
	}

	for (int j = 0; j < vnx; j++) {
		v[0][j] = 0;
		v[vny - 1][j] = 0;
	}
}

double convectionDiffusion(double** w, double u, double v, int i, int j) {
	double result = 0.;

	result += u * (w[i][j + 1] - w[i][j - 1]) / (2 * h) + v * (w[i][j + 1] - w[i][j - 1]) / (2 * h);
	result -= (1 / Re) * ((w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (h * h) + (w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (h * h));

	return result;
};

void addRegularConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -4;
	presureLine[lineNumber * nx + columnNumber + 1] = 1;
	presureLine[lineNumber * nx + columnNumber - 1] = 1;
	presureLine[(lineNumber - 1) * nx + columnNumber] = 1;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx*ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber*nx + columnNumber] = -3;
	presureLine[lineNumber*nx + columnNumber + 1] = 1;
	presureLine[lineNumber*nx + columnNumber - 1] = 1;
	presureLine[(lineNumber - 1)*nx + columnNumber] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -3;
	presureLine[lineNumber * nx + columnNumber + 1] = 1;
	presureLine[lineNumber * nx + columnNumber - 1] = 1;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addLeftMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -3;
	presureLine[lineNumber * nx + columnNumber + 1] = 1;
	presureLine[(lineNumber - 1) * nx + columnNumber] = 1;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addRightMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -3;
	presureLine[lineNumber * nx + columnNumber - 1] = 1;
	presureLine[(lineNumber - 1) * nx + columnNumber] = 1;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}


void addLeftAndTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -2;
	presureLine[(lineNumber - 1) * nx + columnNumber] = 1;
	presureLine[lineNumber * nx + columnNumber + 1] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addLeftAndBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -2;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;
	presureLine[lineNumber * nx + columnNumber + 1] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addRightAndTopMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -2;
	presureLine[(lineNumber - 1) * nx + columnNumber] = 1;
	presureLine[lineNumber * nx + columnNumber - 1] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void addRightAndBottomMissConditionPressure(double** matrix, double* rightPart, int n, int columnNumber, int lineNumber, double** u, double** v) {
	double* presureLine = new double[nx * ny];

	for (int k = 0; k < nx * ny; k++)
	{
		presureLine[k] = 0;
	}

	presureLine[lineNumber * nx + columnNumber] = -2;
	presureLine[(lineNumber + 1) * nx + columnNumber] = 1;
	presureLine[lineNumber * nx + columnNumber - 1] = 1;

	matrix[n] = presureLine;

	rightPart[n] = (h / t) * (u[lineNumber][columnNumber + 1] - u[lineNumber][columnNumber] + v[lineNumber + 1][columnNumber] - v[lineNumber][columnNumber]);
}

void Pard(double** pMatrix, double* rightPart, int N, double* p)
{

	int* ia, * ja;
	double* A, * b, * x;
	ia = new int[N + 1];
	ja = new int[N * N];
	A = new double[N * N];
	b = new double[N];
	x = new double[N];

	double xk, yk, xj, yj, nyk, nxk;

	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	MKL_INT nrhs = 1;
	void* pt[64];
	double ddum;
	MKL_INT idum;
	MKL_INT mtype;

	//------------------------------------------------------
	for (int i = 0; i < N + 1; i++)
		ia[i] = i * N;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			ja[i * N + j] = j;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A[i * N + j] = pMatrix[i][j];
		}
	}

	for (int j = 0; j < N; j++)
		b[j] = rightPart[j];

	//------------------------------------------------------

	mtype = 11;       /* Real unsymmetric matrix */

/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (int i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;		  /* Zero-based numeration */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;         /* Which factorization to use. */
	msglvl = 0;           /* Print statistical information  */
	error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N, A, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}

	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N, A, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}

	phase = 33;

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N, A, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&N, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	for (int i = 0; i < N; i++)
		p[i] = x[i];

	delete[] x;
	delete[] b;
	delete[] ia;
	delete[] ja;
	delete[] A;
}