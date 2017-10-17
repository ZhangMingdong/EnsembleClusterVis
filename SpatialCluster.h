#pragma once
#include <vector>
#include <MathTypes.hpp>

class UncertaintyRegion {
public:
	// clustered point
	std::vector<IPoint2> _vecPoints;
	// cluster area, number of grid points
	int _nArea = 0;
};

/*
	class representing a spatial cluster
*/
class SpatialCluster
{
public:
	SpatialCluster();
	~SpatialCluster();
private:
	// width of the grid
	int _nW = 0;
	// height of the grid
	int _nH = 0;
	// clutering threshold
	double _dbThreshold = 0;
	// the grid to cluster
	const double* _arrGrid = NULL;
	// store state of the grid points
	int* _arrStates = NULL;
public:
	// do the cluster algorithm
	void DoCluster(const double* arrGrid, int nW, int nH, double dbThreshold, std::vector<UncertaintyRegion>& result);
private:
	// spread for a grid point
	void spread(int nI, int nJ, int nCluster, UncertaintyRegion& cluster);
};



