#include "Sequence2DJeung.h"

#include <algorithm>
#include <iterator>
#include <iostream>

using namespace std;

Sequence2DJeung::Sequence2DJeung()
{
}

Sequence2DJeung::~Sequence2DJeung()
{
}

// detect the trend
void Sequence2DJeung::TrendDetectDB() {
	TrendDetect();
}

// detect the trend
void Sequence2DJeung::TrendDetect() {

	sortIndices();

	groupTimeSteps();

	trendDetectBasedOnGroups();

	// print the result
	printGroup();
}

// get the event list
std::vector<DBEpsilonEvent> Sequence2DJeung::GetEvents() {
	vector<DBEpsilonEvent> vecEvents;
	return vecEvents;
}

void Sequence2DJeung::groupTimeSteps() {
	for each (vector<IndexAndValue> vecIV in _vectOrderedIndex)
	{
		vector<set<int>> vecC;
		// for each time steps
		for (size_t i = 0; i < _nEns; i++)
		{
			if (i==0|| vecIV[i]._dbValue-vecIV[i-1]._dbValue>_dbEpsilon)
			{
				// add a new group
				set<int> setC;
				setC.insert(vecIV[i]._nIndex);
				vecC.push_back(setC);
			}
			else {
				vecC[vecC.size() - 1].insert(vecIV[i]._nIndex);
			}
		}
		_vecClusters.push_back(vecC);
	}
}

// add a group gNew to list vecCandidate
// remove duplicate
void AddGroupToList(Group gNew, vector<Group> &vecCandidate) {
	for (Group& g : vecCandidate) {
		if (g._member==gNew._member)
		{
			if (gNew._nS < g._nS) g._nS = gNew._nS;
			return;
		}
	}
	vecCandidate.push_back(gNew);
}

void Sequence2DJeung::trendDetectBasedOnGroups() {
	// store the candidates during the process of the calculation
	vector<Group> vecCandidates;

	// first time step -- add the clusters with element more than M to the vector of candidates
	for (set<int> setC : _vecClusters[0]) {
		if (setC.size()>_nM)
		{
			Group g;
			g._nS = g._nE = 0;
			g._member = setC;
			vecCandidates.push_back(g);
		}
	}

	// subsequent time steps
	for (size_t i = 1; i < _nLength; i++)
	{
		vector<Group> vecCandidatesNew;
		int nClusters = _vecClusters[i].size();

		// set assigned of all the clusters false
		vector<bool> vecAssigned(nClusters);
		for (size_t j = 0; j < nClusters; j++)
		{
			vecAssigned[j] = false;
		}

		// for each candidate -- check the intersection with all the clusters
		for (Group& candidate : vecCandidates) {		
			candidate._bAssigned = false;
			// for each cluster
			for (int iC = 0; iC < nClusters; iC++) {
				set<int> setC = _vecClusters[i][iC];
				if (!setC.size() > _nM) continue;
				// calculate intersection
				std::set<int> common_data;
				set_intersection(setC.begin(), setC.end(), candidate._member.begin(), candidate._member.end(),
					inserter(common_data, common_data.begin()));
				if (common_data.size()>_nM)
				{
					if (common_data == candidate._member)
					{
						// prolong the candidate
						candidate._nE = i;
						AddGroupToList(candidate, vecCandidatesNew);
						//vecCandidatesNew.push_back(candidate);
						candidate._bAssigned = true;
					}
					else {
						// add new candidate for the intersection
						Group gNew;
						gNew._nS = gNew._nE = candidate._nS;
						gNew._member = common_data;
						AddGroupToList(gNew, vecCandidatesNew);
						//vecCandidatesNew.push_back(gNew);
					}
					if (common_data==setC)
					{
						vecAssigned[iC] = true;
					}
				}
			}
		}
		// add all the unassigend candidate to the group list
		for (Group candidate : vecCandidates) {
			if (!candidate._bAssigned)
			{
				candidate._nE = i-1;
				addGroup(candidate);
			}
		}
		// add all the unassigned clusters to the candidate list
		for (int iC = 0; iC < nClusters; iC++) {
			if (!vecAssigned[iC]&& _vecClusters[i][iC].size()>_nM)
			{
				Group gNew;
				gNew._nE = gNew._nS = i;
				gNew._member = _vecClusters[i][iC];
				AddGroupToList(gNew, vecCandidatesNew);
				//vecCandidatesNew.push_back(gNew);
			}
		}

		vecCandidates = vecCandidatesNew;
	}


	// add all the unassigend candidate to the group list
	for (Group candidate : vecCandidates) {
		if (!candidate._bAssigned)
		{
			addGroup(candidate);
		}
	}
}