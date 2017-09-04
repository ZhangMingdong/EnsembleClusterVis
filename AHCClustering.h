// ============================================================
// agglomerative hierarchical clustering (AHC)
// Written by Mingdong
// 2017/05/12
// ============================================================
#pragma once
#include "Clustering.h"
class AHCClustering : public Clustering
{
public:
	AHCClustering();
	virtual ~AHCClustering();
protected:
	virtual void doCluster();
	// build the distance matrix
	void buildDistanceMatrix(double** disMatrix);
private:
	// use the max or min distance between clusters
	void doClusterImpMaxOrMin(bool bMax = true);
	// use the distance of the cluster center
	void doClusterImp();
	// clean the resource
	void clean();
private:
	// distance matrix
	double** _disMatrix;
	// array of states
	int* _arrState;
	// centroid of each cluster[x0,y0,x1,y1,...]
	double* _arrCentroid;
};

