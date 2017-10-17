#include "SpatialCluster.h"
#include <QDebug>


SpatialCluster::SpatialCluster()
{
}


SpatialCluster::~SpatialCluster()
{
}

#define UNPROCESSED -2
#define THRESHOLD	-1

void SpatialCluster::DoCluster(const double* arrGrid, int nW, int nH, double dbThreshold, std::vector<UncertaintyRegion>& result) {
	// 0.set parameter
	_arrGrid = arrGrid;
	_nW = nW;
	_nH = nH;
	_dbThreshold = dbThreshold;
	// 1.initialize states
	_arrStates = new int[_nW*_nH];
	for (size_t i = 0,length= _nW*_nH; i < length; i++)
	{
		_arrStates[i] = UNPROCESSED;
	}
	// process
	int nCurrentCluster = 0;	// current cluster sequence
	for (size_t i = 0; i < _nH; i++)
	{
		for (size_t j = 0; j < _nW; j++) {
			int nIndex = i*_nW + j;
			if (_arrStates[nIndex]== UNPROCESSED)
			{
				if (_arrGrid[nIndex]<_dbThreshold)
				{
					if (false)
					{
						qDebug() << _arrGrid[nIndex];
					}
					_arrStates[nIndex] = THRESHOLD;
				}
				else {
					result.push_back(UncertaintyRegion());
					spread(i, j, nCurrentCluster, result[nCurrentCluster]);
					nCurrentCluster++;
				}
			}
		}
	}

	delete []_arrStates;
}


void  SpatialCluster::spread(int nI, int nJ, int nCluster, UncertaintyRegion& cluster) {
	_arrStates[nI*_nW + nJ] = nCluster;
	cluster._vecPoints.push_back(IPoint2(nI, nJ));
	cluster._nArea++;
	if (nI > 0 && _arrStates[(nI - 1)*_nW + nJ] == UNPROCESSED && _arrGrid[(nI - 1)*_nW + nJ] >_dbThreshold) spread(nI - 1, nJ, nCluster, cluster);
	if (nI <_nH - 1 && _arrStates[(nI + 1)*_nW + nJ] == UNPROCESSED && _arrGrid[(nI + 1)*_nW + nJ] >_dbThreshold) spread(nI + 1, nJ, nCluster, cluster);
	if (nJ > 0 && _arrStates[nI*_nW + nJ - 1] == UNPROCESSED && _arrGrid[nI*_nW + nJ - 1] >_dbThreshold) spread(nI, nJ - 1, nCluster, cluster);
	if (nJ <_nW - 1 && _arrStates[nI*_nW + nJ + 1] == UNPROCESSED && _arrGrid[nI*_nW + nJ + 1] >_dbThreshold) spread(nI, nJ + 1, nCluster, cluster);
}