#include "AHCClustering.h"
#include <math.h>

#include <QDebug>


namespace CLUSTER {


	AHCClustering::AHCClustering()
	{
	}


	AHCClustering::~AHCClustering()
	{
	}

	void AHCClustering::doCluster() {
		// 0.initialize the distance matrix
		_disMatrix = new double*[_n];
		buildDistanceMatrix(_disMatrix);
		// 1.initialized the centroid
		_arrCentroid = new double[_n * _m];
		for (size_t i = 0; i < _n*_m; i++)
		{
			_arrCentroid[i] = _arrData[i];
		}

		// 2.clustering
		for (int i = 0; i < _n; i++)
		{
			_arrLabels[i] = i;
		}
		// state of each grid point
		_arrState = new int[_n];
		for (int i = 0; i < _n; i++) _arrState[i] = 1;	// 0 means this grid has been merged to other cluster

		bool bNewImp = true;
		if (bNewImp)
		{
			doClusterImp();
		}
		else {
			doClusterImpMaxOrMin();
		}


		// 3.release the resource
		clean();

		//	realign the label
		int nBias = 0;
		for (int i = 0; i < _n; i++)
		{
			// if i is a valid label
			bool bValidLabel = false;
			for (int j = 0; j < _n; j++)
			{
				if (_arrLabels[j] == i) {
					_arrLabels[j] -= nBias;
					bValidLabel = true;
				}
			}
			if (!bValidLabel) nBias++;
		}
	}

	void AHCClustering::doClusterImpMaxOrMin(bool bMax) {
		// start from every grid point as a cluster
		int nClusters = _n;

		while (nClusters > _k)
		{
			// 2.1.Find the nearest points pair nI,nJ
			int nI, nJ;
			double fMinDis = 10000.0;
			for (int i = 0; i < _n; i++)
			{
				if (_arrState[i] == 0) continue;

				for (int j = 0; j < i; j++)
				{
					if (_arrState[j] == 0) continue;
					if (_arrLabels[i] == _arrLabels[j]) continue;
					if (_disMatrix[i][j] < fMinDis)
					{
						nI = i;
						nJ = j;
						fMinDis = _disMatrix[i][j];
					}
				}
			}

			int nSourceLabel = _arrLabels[nI];
			int nTargetLabel = _arrLabels[nJ];


			// 2.2.search all the point with source label, reset distance and state
			for (int i = 0; i < _n; i++)
			{
				if (_arrLabels[i] == nSourceLabel)
				{
					// 2.2.1.reset the label
					_arrLabels[i] = nTargetLabel;
					if (_arrState[i] == 1)
					{
						// 2.2.2.reset distance
						for (int j = 0; j < _n; j++)
						{
							if (j == i || j == nJ) continue;
							// as we merge i to nJ
							// for each j, if dis(i,j)>dis(nJ,j), we use dis(i,j) to replace dis(nJ,j)
							double fDisIJ;
							if (i < j) fDisIJ = _disMatrix[j][i];
							else fDisIJ = _disMatrix[i][j];
							if (j < nJ) {
								_disMatrix[nJ][j] = (bMax ? fmax(_disMatrix[nJ][j], fDisIJ) : fmin(_disMatrix[nJ][j], fDisIJ));
							}
							else {
								_disMatrix[j][nJ] = (bMax ? fmax(_disMatrix[j][nJ], fDisIJ) : fmin(_disMatrix[j][nJ], fDisIJ));
							}
						}
						_arrState[i] = 0;
					}
				}
			}
			nClusters--;
		}
	}

	void AHCClustering::doClusterImp() {
		// start from every grid point as a cluster
		int nClusters = _n;

		while (nClusters > _k)
		{
			// 1.Find the nearest points pair nI,nJ
			int nI, nJ;
			double fMinDis = 1000000.0;
			for (int i = 0; i < _n; i++)
			{
				if (_arrState[i] == 0) continue;

				for (int j = 0; j < i; j++)
				{
					if (_arrState[j] == 0) continue;
					if (_arrLabels[i] == _arrLabels[j]) continue;
					if (_disMatrix[i][j] < fMinDis)
					{
						nI = i;
						nJ = j;
						fMinDis = _disMatrix[i][j];
					}
				}
			}

			int nSourceLabel = _arrLabels[nI];
			int nTargetLabel = _arrLabels[nJ];


			// 2.search all the point with source label, reset the label
			for (int i = 0; i < _n; i++)
			{
				if (_arrLabels[i] == nSourceLabel)
				{
					// 2.2.1.reset the label
					_arrLabels[i] = nTargetLabel;
					_arrState[i] = 0;
				}
			}
			// 2.reset the distance
			// calculate centroid
			vectorReset(_arrCentroid + nJ*_m, _m);
			int num = 0;
			for (size_t i = 0; i < _n; i++)
			{
				if (_arrLabels[i] == nTargetLabel)
				{
					vectorAdd(_arrCentroid + i*_m, _arrData + i*_m, _m);
					num++;
				}
			}
			vectorDiv(_arrCentroid, num, _m);
			// calculate distance
			for (size_t i = 0; i < _n; i++)
			{
				if (_arrState[i] == 0) continue;
				for (size_t j = 0; j < i; j++)
				{
					if (_arrState[j] == 0) continue;
					if (_arrLabels[i] == nTargetLabel || _arrLabels[j] == nTargetLabel&&_arrLabels[i] != _arrLabels[j])
					{
						int index1 = i;
						int index2 = j;
						if (i < j)
						{
							index1 = j;
							index2 = i;
						}
						_disMatrix[index1][index2] = vectorDis(_arrCentroid + index1*_m, _arrCentroid + index2*_m, _m);
					}
				}
			}
			nClusters--;
		}


	}

	void AHCClustering::buildDistanceMatrix(double** disMatrix) {
		// 0.allocate resource
		for (int i = 0; i < _n; i++)
		{
			disMatrix[i] = new double[_n];
		}
		// 1.calculate the distance
		for (int i = 0; i < _n; i++)
		{
			for (int j = 0; j < i; j++)
			{

				disMatrix[i][j] = vectorDis(_arrData + i*_m, _arrData + j*_m, _m);
			}
		}
	}

	void AHCClustering::clean() {
		for (int i = 0; i < _n; i++)
		{
			delete[]_disMatrix[i];
		}
		delete[] _disMatrix;
		delete[] _arrState;
		delete[] _arrCentroid;
	}

}