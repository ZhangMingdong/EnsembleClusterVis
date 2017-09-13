// ============================================================
// DBSCAN clustering
// Written by Mingdong
// do not use clustering interface
// 2017/05/26
// ============================================================
#pragma once
#include "Clustering.h"
namespace CLUSTER {
	class DBSCANClustering : public Clustering
	{
	public:
		DBSCANClustering();
		virtual ~DBSCANClustering();
	private:
		// clustering parameters
		int _nMinPts=100;
		double _dbEps=.5;
	protected:
		virtual void doCluster();
	public:		
		virtual void SetDBSCANParams(int minPts, double eps) { _nMinPts = minPts; _dbEps = eps; };
	};

}