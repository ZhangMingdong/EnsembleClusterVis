#include "Sequence2DBuchin.h"
#include <ctime>
#include <iostream>
#include <algorithm>

using namespace std;

Sequence2DBuchin::Sequence2DBuchin()
{
}

Sequence2DBuchin::~Sequence2DBuchin()
{
}

void Sequence2DBuchin::TrendDetect() {

}

void Sequence2DBuchin::TrendDetectDB() {
	// sort the indices
	sortIndices();

	// find all the dbEvents();
	findDBEvents();
	
	// generate groups
	detectGroupsFromDBEvents();

	// print the calculated result
	printDBGroup();
}

void Sequence2DBuchin::findDBEvents() {
	// search each time step, from the 2nd one
	for (size_t iStep = 1; iStep < _nLength; iStep++)
	{
		// check each index
		for (size_t iSeq1 = 0; iSeq1 < _nEns; iSeq1++)
		{
			// index and value of seq1
			IndexAndValue iv11 = _vectOrderedIndex[iStep][iSeq1];
			// index and value of seq1 in last time step
			IndexAndValue iv10 = _vectOrderedIndex[iStep-1][_vecRanks[iStep - 1][iv11._nIndex]];
			for (size_t iSeq2 = 0; iSeq2 < iSeq1; iSeq2++)
			{
				if (iSeq1 == iSeq2) continue;
				IndexAndValue iv21 = _vectOrderedIndex[iStep][iSeq2];
				IndexAndValue iv20 = _vectOrderedIndex[iStep - 1][_vecRanks[iStep - 1][iv21._nIndex]];
				bool bNear1 = abs(iv11._dbValue - iv21._dbValue) < _dbEpsilon;
				bool bNear0 = abs(iv10._dbValue - iv20._dbValue) < _dbEpsilon;
				if (bNear0&&!bNear1)
				{
					// split event
					double dbTime = getEventTime(iStep-1, iStep, iv11._nIndex, iv21._nIndex);
					if (isDBEventFree(dbTime, iv11._nIndex, iv21._nIndex) == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(1, dbTime, iv11._nIndex, iv21._nIndex));
					}
				}
				else if (bNear1&&!bNear0) {
					// merge event
					double dbTime = getEventTime(iStep - 1, iStep, iv11._nIndex, iv21._nIndex);
					if (isDBEventFree(dbTime, iv11._nIndex, iv21._nIndex) == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(0, dbTime, iv11._nIndex, iv21._nIndex));
					}
				}
				else if (!bNear0&&!bNear1 && (iv11._dbValue - iv21._dbValue)*(iv10._dbValue - iv20._dbValue) < 0) {
					// merge and split
					double dbTime1;
					double dbTime2;
					getCrossEventTimes(dbTime1, dbTime2, iStep - 1, iStep, iv11._nIndex, iv21._nIndex);
					// merge event
					if (isDBEventFree(dbTime1, iv11._nIndex, iv21._nIndex) == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(0, dbTime1, iv11._nIndex, iv21._nIndex));
					}
					// split event
					if (isDBEventFree(dbTime2, iv11._nIndex, iv21._nIndex) == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(1, dbTime2, iv11._nIndex, iv21._nIndex));
					}
				}
			}
		}
	}
	// sort them
	sort(_vecFreeDBEvents.begin(), _vecFreeDBEvents.end(), DBEventCompare);
	
	
	cout << "number of events: " << _vecFreeDBEvents.size() << endl;
	for (DBEpsilonEvent evt : _vecFreeDBEvents) {
		cout << "(" << evt._dbTime << "," << evt._nType << "," << evt._nIndex1 << "," << evt._nIndex2 << "," << ")" << endl;
	}
}

void Sequence2DBuchin::findDBEvents_1() {
	// search each time step, from the 2nd one
	for (size_t iStep = 0; iStep < _nLength; iStep++)
	{
		// check each index
		for (size_t iSeq1 = 0; iSeq1 < _nEns - 1; iSeq1++)
		{
			size_t iSeq2 = iSeq1 + 1;
			IndexAndValue iv1 = _vectOrderedIndex[iStep][iSeq1];
			IndexAndValue iv2 = _vectOrderedIndex[iStep][iSeq2];
			IndexAndValue ivBefore2 = _vectOrderedIndex[iStep][iSeq1];
			while (iv2._dbValue - ivBefore2._dbValue<_dbEpsilon)
			{
				// while the two indices are directly epsilon-connected
				// check for start event
				if (iStep>0 && abs(_vecSequences[iv1._nIndex][iStep - 1] - _vecSequences[iv2._nIndex][iStep - 1])>_dbEpsilon)	// apart in last step
				{
					double dbTime = getEventTime(iStep - 1, iStep, iv1._nIndex, iv2._nIndex);
					int nEventFreeState = isDBEventFree(dbTime, iv1._nIndex, iv2._nIndex);
					if (nEventFreeState == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(0, dbTime, iv1._nIndex, iv2._nIndex));
					}
					else if (nEventFreeState == -1) {
						_vecFreeDBEvents.push_back(DBEpsilonEvent(0, dbTime, iv2._nIndex, iv1._nIndex));

					}

				}
				// check for end event
				if (iStep<_nLength - 1 && abs(_vecSequences[iv1._nIndex][iStep + 1] - _vecSequences[iv2._nIndex][iStep + 1])>_dbEpsilon)	// apart in next step
				{
					double dbTime = getEventTime(iStep, iStep + 1, iv1._nIndex, iv2._nIndex);
					int nEventFreeState = isDBEventFree(dbTime, iv1._nIndex, iv2._nIndex);
					if (nEventFreeState == 1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(1, dbTime, iv1._nIndex, iv2._nIndex));
					}
					else if (nEventFreeState == -1)
					{
						_vecFreeDBEvents.push_back(DBEpsilonEvent(1, dbTime, iv2._nIndex, iv1._nIndex));
					}
				}
				// loop
				iSeq2++;
				if (iSeq2 == _nEns)
				{
					break;
				}
				else {
					ivBefore2 = iv2;
					iv2 = _vectOrderedIndex[iStep][iSeq2];
				}
			}
		}
	}
	// sort them
	sort(_vecFreeDBEvents.begin(), _vecFreeDBEvents.end(), DBEventCompare);


	cout << "number of events: " << _vecFreeDBEvents.size() << endl;
	for (DBEpsilonEvent evt : _vecFreeDBEvents) {
		cout << "(" << evt._dbTime << "," << evt._nType << "," << evt._nIndex1 << "," << evt._nIndex2 << "," << ")" << endl;
	}
}

int Sequence2DBuchin::isDBEventFree(double dbTime, int nIndex1, int nIndex2) {
	double dbV1 = getValue(nIndex1, dbTime);
	double dbV2 = getValue(nIndex2, dbTime);
	if (dbV1 > dbV2) {
		double dbTemp = dbV1;
		dbV1 = dbV2;
		dbV2 = dbTemp;
	}
	for (size_t i = 0; i < _nEns; i++)
	{
		if (i != nIndex1 && i != nIndex2) {
			double dbValue=getValue(i, dbTime);
			if (dbValue > dbV1&&dbValue < dbV2) return 0;
		}
	}	
	return dbV1<dbV2? 1:-1;
}

void Sequence2DBuchin::detectGroupsFromDBEvents() {
	// build groups at time 0;
	vector<DBGroup> vecCandidates;
	int nCurrentGroup = 0;
	vecCandidates.push_back(DBGroup());
	vecCandidates[0]._dbS = 0;
	vecCandidates[0]._member.insert(_vectOrderedIndex[0][0]._nIndex);
	for (size_t i = 1; i < _nEns; i++)
	{
		if (_vectOrderedIndex[0][i]._dbValue - _vectOrderedIndex[0][i - 1]._dbValue>_dbEpsilon) {
			vecCandidates.push_back(DBGroup());
			nCurrentGroup++;
			vecCandidates[nCurrentGroup]._dbS = 0;
		}
		vecCandidates[nCurrentGroup]._member.insert(_vectOrderedIndex[0][i]._nIndex);
	}
	// treat every events
	for each (DBEpsilonEvent evt in _vecFreeDBEvents)
	{
		if (evt._nType==0){
			// merge event - a new candidate will be added
			DBGroup gNew;
			gNew._dbS = evt._dbTime;
			for each (DBGroup g in vecCandidates)
			{
				// if a group contain one of the members, add all its members into the new group
				if (g._member.find(evt._nIndex1) != g._member.end() || g._member.find(evt._nIndex2) != g._member.end()) {
					for each (int nIndex in g._member)
					{
						gNew._member.insert(nIndex);
					}
				}
			}
			vecCandidates.push_back(gNew);
		}
		else {	// split event		
			// get the order of items at this time
			vector<IndexAndValue> vecIV;
			for (size_t i = 0; i < _nEns; i++)
			{
				vecIV.push_back(IndexAndValue(i, getValue(i, evt._dbTime)));
			}
			sort(vecIV.begin(), vecIV.end(), IndexAndValueCompare);

			// split the candidates
			vector<DBGroup> vecNewCandidates;
			for (int i = vecCandidates.size() - 1; i >= 0;i--)
			{
				DBGroup g = vecCandidates[i];
				// if a group contain at least of the two members
				if (g._member.find(evt._nIndex1) != g._member.end() || g._member.find(evt._nIndex2) != g._member.end()) {
					DBGroup gNew[2];
					gNew[0]._dbS = g._dbS;
					gNew[1]._dbS = g._dbS;
					int nCurrentGroup = 0;
					for each (IndexAndValue iv in vecIV)
					{
						if (g._member.find(iv._nIndex)!=g._member.end())
						{
							gNew[nCurrentGroup]._member.insert(iv._nIndex);
						}
						if (nCurrentGroup == 0 && (iv._nIndex == evt._nIndex1 || iv._nIndex == evt._nIndex2)) {
							nCurrentGroup++;
						}
					}
					if (gNew[0]._member.size()>0&&gNew[1]._member.size()>0)
					{
						vecNewCandidates.push_back(gNew[0]);
						vecNewCandidates.push_back(gNew[1]);
						g._dbE = evt._dbTime;
						addDBGroup(g);
						vecCandidates.erase(vecCandidates.begin() + i);
					}				
				}

			}
			// add new candidate
			for each (DBGroup gNew in vecNewCandidates)
			{
				bool bDuplicate = false;
				for(DBGroup& g : vecCandidates)
				{
					if (g._member == gNew._member) {
						if (gNew._dbS<g._dbS)
						{
							g._dbS = gNew._dbS;
						}
						bDuplicate = true;
						break;
					}
				}
				if (!bDuplicate)
				{
					vecCandidates.push_back(gNew);
				}
			}
		}
	}

	// add all the candidates when finished all the events
	for each (DBGroup g in vecCandidates)
	{
		g._dbE = _nLength - 1;
		addDBGroup(g);
	}
}