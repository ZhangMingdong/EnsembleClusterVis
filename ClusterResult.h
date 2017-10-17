#pragma once
#include <vector>
#include "def.h"

// record the cluster result
struct ClusterResult {
	int _nM = 0;								// dimension of each item. also length of _arrLabels;
	int _nK = 0;								// number of clusters 1~5;
	std::vector<int> _vecItems[g_nClusterMax];	// items of each cluster, only the first _nK are used
	int _arrLabels[g_nEnsembles];				// labels of each ensemble member
												// push a label
	void PushLabel(int nIndex, int nLabel);
	// match this cluster index with c, minimized the items that changing cluster index
	void Match(ClusterResult& c);
	// align with a map
	void AlighWith(int* arrMap);
	// reset the result
	void Reset(int nM, int nK);
	// sort clusters by number of items in them
	void Sort();
	// generate items by the array of labels
	void generateItemsByLabels();

	/*
		allign these result
		nK: number of clusters
	*/
	static void Align(ClusterResult* arrResult, int nK);

	/*
	allign these result
	nK: number of clusters
	nBase: alignment base
	*/
	static void Align(ClusterResult* arrResult, int nK,int nBase);
};