#include "Clustering.h"
#include <math.h>




namespace CLUSTER {

	Clustering::Clustering()
	{
	}


	Clustering::~Clustering()
	{
	}

	int Clustering::DoCluster(int n, int m, int k, const double *arrData, int* arrLabels) {
		_n = n;
		_m = m;
		_k = k;
		_arrData = arrData;
		_arrLabels = arrLabels;
		doCluster();
		return _k;
	}

	// add vector v to vector vSum
	void vectorAdd(double* vSum, const double* v, int nLen) {
		for (size_t i = 0; i < nLen; i++)
		{
			vSum[i] += v[i];
		}
	}

	// divide a vector
	void vectorDiv(double* v, double dbScale, int nLen) {
		for (size_t i = 0; i < nLen; i++)
		{
			v[i] /= dbScale;
		}
	}

	// reset a vector
	void vectorReset(double*v, int nLen) {
		for (size_t i = 0; i < nLen; i++)
		{
			v[i] = 0;
		}
	}

	// calculate distance between the two vectors
	double vectorDis(const double *v1, const double* v2, int nLen) {
		double sum = 0;

		for (size_t i = 0; i < nLen; i++)
		{
			double bias = v1[i] - v2[i];
			sum += bias*bias;
		}
		return sqrt(sum);
	}


	// copy vector v1 to v0
	void vectorCopy(double* v0, const double* v1, int nLen) {
		for (size_t i = 0; i < nLen; i++)
		{
			v0[i] = v1[i];
		}
	}
}