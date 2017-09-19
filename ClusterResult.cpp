#include "ClusterResult.h"

#include <iostream>

using namespace std;

/*
allign these result
nK: number of clusters
arrMap: map of the ensemble members to its new sequence
*/
void ClusterResult::Align(ClusterResult* arrResult, int nK) {
	int arrMap[g_nEnsembles];
	// 0.initialize the map to natral sequence
	for (size_t i = 0; i < g_nEnsembles; i++) arrMap[i] = i;

	// 1.swap sorting
	for (size_t i = 0; i < g_nEnsembles - 1; i++)
	{
		for (size_t j = i + 1; j < g_nEnsembles; j++) {
			bool bSwap = false;
			for (size_t k = 0; k < nK; k++)
			{
				if (arrResult[k]._arrLabels[arrMap[i]]>arrResult[k]._arrLabels[arrMap[j]]) {
					bSwap = true;
					break;
				}
				else if (arrResult[k]._arrLabels[arrMap[i]] < arrResult[k]._arrLabels[arrMap[j]]) {
					break;
				}
			}
			if (bSwap)
			{
				int nTemp = arrMap[i];
				arrMap[i] = arrMap[j];
				arrMap[j] = nTemp;
			}
		}
	}

	// 2.reset the ordering
	int arrMap2[g_nEnsembles];
	for (size_t i = 0; i < g_nEnsembles; i++)
	{
		arrMap2[arrMap[i]] = i;
	}

	for (size_t i = 0; i < nK; i++)
	{
		arrResult[i].AlighWith(arrMap2);
	}
}

void ClusterResult::PushLabel(int nIndex, int nLabel) {
	_vecItems[nLabel].push_back(nIndex);
	_arrLabels[nIndex] = nLabel;
}


void ClusterResult::Match(ClusterResult& mc) {
	int mtxWeight[5][5];					// weight matrix
	int arrMatch[5];						// the result of the perfect match
	int nClusters = 5;						// number of clusters
											// 0.initialize matrix
	for (size_t i = 0; i < nClusters; i++)
	{
		for (size_t j = 0; j < nClusters; j++) {
			mtxWeight[i][j] = 0;
		}
	}
	// 1.calculate weight
	for (size_t i = 0; i < _nM; i++)
	{
		mtxWeight[_arrLabels[i]][mc._arrLabels[i]]++;
	}

	for (size_t i = 0; i < nClusters; i++)
	{
		for (size_t j = 0; j < nClusters; j++) {
			cout << mtxWeight[i][j] << "\t";
		}
		cout << endl;
	}

	// 2.calculate max weight perfect matching using Hungarian algorithm
	// use approximate algorithm, always searching for the biggest value
	for (size_t l = 0; l < nClusters; l++) {
		// find biggest value
		int iMax = 0;
		int jMax = 0;
		int vMax = -1;
		for (size_t i = 0; i < nClusters; i++)
		{
			for (size_t j = 0; j < nClusters; j++) {
				if (mtxWeight[i][j]>vMax) {
					vMax = mtxWeight[i][j];
					iMax = i;
					jMax = j;
				}
			}
		}
		arrMatch[iMax] = jMax;
		for (size_t i = 0; i < nClusters; i++)
		{
			mtxWeight[i][jMax] = mtxWeight[iMax][i] = -1;
		}
	}
	/*
	qDebug() << "match:";
	for (size_t i = 0; i < nClusters; i++)
	{
	qDebug() << arrMatch[i];
	}
	*/
	// 3.reset labels
	for (size_t i = 0; i < _nM; i++)
	{
		_arrLabels[i] = arrMatch[_arrLabels[i]];
	}
	generateItemsByLabels();

}


void ClusterResult::AlighWith(int* arrMap) {
	for (size_t i = 0; i < _nK; i++)
	{
		int nLen = _vecItems[i].size();
		for (size_t j = 0; j < nLen - 1; j++)
		{
			for (size_t k = j + 1; k < nLen; k++) {
				if (arrMap[_vecItems[i][j]]>arrMap[_vecItems[i][k]])
				{
					int nTemp = _vecItems[i][j];
					_vecItems[i][j] = _vecItems[i][k];
					_vecItems[i][k] = nTemp;
				}
			}
		}
	}
}

void ClusterResult::generateItemsByLabels() {
	for (size_t i = 0; i < 5; i++)
	{
		_vecItems[i].clear();
	}
	for (size_t i = 0; i < _nM; i++)
	{
		_vecItems[_arrLabels[i]].push_back(i);
	}
}

void ClusterResult::Sort() {
	// array recording count of each cluster
	int arrCount[g_nClusterMax];
	for (size_t i = 0; i < g_nClusterMax; i++)
	{
		arrCount[i] = _vecItems[i].size();
	}
	// get the projection
	int arrMap[g_nClusterMax];
	for (size_t i = 0; i < 5; i++)
	{
		int nMax = -1;
		int nMaxIndex = -1;
		for (size_t j = 0; j < g_nClusterMax; j++)
		{
			if (arrCount[j] > nMax) {
				nMax = arrCount[j];
				nMaxIndex = j;
			}
		}
		arrMap[nMaxIndex] = i;
		arrCount[nMaxIndex] = -1;
	}
	// reset label
	for (size_t i = 0; i < _nM; i++)
	{
		_arrLabels[i] = arrMap[_arrLabels[i]];
	}
	generateItemsByLabels();
}

void ClusterResult::Reset(int nM, int nK) {
	_nM = nM;
	_nK = nK;
	for (int i = 0; i < 5; i++) {
		_vecItems[i].clear();
	}
}