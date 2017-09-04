// ============================================================
// DBSCAN clustering
// Written by Mingdong
// do not use clustering interface
// 2017/05/26
// ============================================================
#pragma once
#include "Clustering.h"
class DBSCANClustering// : public Clustering
{
public:
	DBSCANClustering();
	virtual ~DBSCANClustering();

public:
	/*
		do the clustering
		h,w: the height and width of the clustering field
		nMinPts, dbEps: clustering paramter
		arrData: input data, 0/1, 1 means this point is to cluster
		return the number of clusters
	*/
	int DoCluster(int h, int w, int nMinPts, double dbEps, const int *arrData, int* arrLabels);

private:
	// clustering parameters
	int _nH;
	int _nW;
	int _nMinPts;
	double _dbEps;
	const int* _arrData;
	int* _arrLabels;

	// indicate the current label
	int _nCurrentlabel;

private:
	int doCluster();
};

