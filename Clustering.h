// ============================================================
// Interface for clustering
// Written by Mingdong
// 2017/05/12
// ============================================================
#pragma once
namespace CLUSTER {
	class Clustering
	{
	public:
		Clustering();
		virtual ~Clustering();
	protected:
		int _n;						// number of items
		int _m;						// dimensions of items
		int _k;						// number of clusters
		const double *_arrData;		// data items
		int* _arrLabels;			// labels
	protected:
		/*
		implementation of clustering
		*/
		virtual void doCluster() = 0;
	public:
		/*
			set input data
			n:			number of vectors
			m:			length of vectors
			k:			number of clusters
			arrData:	data buffer
			arrLabels:	return the labels. [length:n][0~r-1]

			return number of clusters
		*/
		int DoCluster(int n, int m, int k, const double *arrData, int* arrLabels);

		// set parameters for DBSCAN
		virtual void SetDBSCANParams(int minPts, double eps) {};
	};

	// add vector v to vector vSum
	void vectorAdd(double* vSum, const double* v, int nLen);

	// divide a vector
	void vectorDiv(double* v, double dbScale, int nLen);

	// reset a vector
	void vectorReset(double*v, int nLen);

	// calculate distance between the two vectors
	double vectorDis(const double *v1, const double* v2, int nLen);


	// copy vector v1 to v0
	void vectorCopy(double* v0, const double* v1, int nLen);
}