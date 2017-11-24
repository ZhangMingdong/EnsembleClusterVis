#include "Sequence2D.h"
#include <algorithm>
#include <iostream>
#include <ctime>

using namespace std;


// map a group to a unified hash number
unsigned long hashFunc(const Group& g) {
	unsigned long ulNum = 0;
	unsigned long ulBase = 1;
	for each (int index in g._member)
	{
		ulNum += (ulBase << index);
	}
	return ulNum;
}

// map a group to a unified hash number
unsigned long hashFunc(const vector<IndexAndValue> vecIV) {
	unsigned long ulNum = 0;
	unsigned long ulBase = 1;
	for each (IndexAndValue iv in vecIV)
	{
		ulNum += (ulBase << iv._nIndex);
	}
	return ulNum;
}

Sequence2D::Sequence2D()
{
}

Sequence2D::~Sequence2D()
{
}

void Sequence2D::AddSequence(std::vector<double> seq) {
	_vecSequences.push_back(seq);
}

void Sequence2D::Init(int nLength, double dbMin, double dbMax, double dbEpsilon) {
	_nLength = nLength;
	_dbMin = dbMin;
	_dbMax = dbMax;
	_dbEpsilon = dbEpsilon;
	_nEns = _vecSequences.size();
}

void Sequence2D::TrendDetectRootOnly() {
	// candidates of groups
	Group group0;
	group0._nS = 0;
	group0._nE = 0;
	for (size_t i = 0; i < _nEns; i++)
	{
		group0._member.push_back(i);
	}

	// start from the first step
	trendDetectStep_SplitOnly(group0, 0);

	// print the result
	printGroup();

}

void Sequence2D::printGroup() {
	for each (Group g in _vecGroups)
	{
#ifdef REVERSE
		cout << _nWidth - 1 - g._nE << "\t" << _nWidth - 1 - g._nS << "\t";
#else
		cout <<"(" << g._nS << "," << g._nE << "):";
#endif		
		for each (int nMember in g._member)
		{
			cout << nMember << ",";
		}
		cout << endl;
	}
}

void Sequence2D::TrendDetect_2() {
	// start from each position
	for (size_t i = 0; i < _nLength - _nDelta; i++)
	{
		Group group0;
		group0._nS = i;
		for (size_t i = 0; i < _nEns; i++)
		{
			group0._member.push_back(i);
		}

		// start from the first step

		// 1.generate segments
		vector<vector<IndexAndValue>> vecSegments = generateSegments(group0, i);
		for each (vector<IndexAndValue> seg in vecSegments)
		{
			startTrendDetect(seg, i);
		}
	}
	// remove duplicated
	for (int i = _vecGroups.size() - 1; i >= 0; i--)
	{
		for (int j = 0, length = _vecGroups.size(); j < length; j++) {
			if (i != j&& _vecGroups[j].Contain(_vecGroups[i]))
			{
				_vecGroups.erase(_vecGroups.begin() + i);
				break;
			}
		}
	}
	// print the result
	printGroup();
}

void Sequence2D::TrendDetect_1() {
	// start from each position
	for (size_t i = 0; i < _nLength - _nDelta; i++)
	{
		Group group0;
		group0._nS = i;
		for (size_t i = 0; i < _nEns; i++)
		{
			group0._member.push_back(i);
		}

		// start from the first step
		trendDetectStep(group0, i);
	}
	// remove duplicated
	for (int i = _vecGroups.size() - 1; i >= 0; i--)
	{
		for (int j = 0, length = _vecGroups.size(); j < length; j++) {
			if (i != j&& _vecGroups[j].Contain(_vecGroups[i]))
			{
				_vecGroups.erase(_vecGroups.begin() + i);
				break;
			}
		}
	}
	// print the result
	for each (Group g in _vecGroups)
	{
#ifdef REVERSE
		cout << _nWidth - 1 - g._nE << "\t" << _nWidth - 1 - g._nS << "\t";
#else
		cout << g._nS << "\t" << g._nE << "\t";
#endif		
		for each (int nMember in g._member)
		{
			cout << nMember << "\t";
		}
		cout << endl;
	}
}

// cluster comparison function, used for sort
bool IndexAndValueCompare(IndexAndValue iv1, IndexAndValue iv2) {
	return iv1._dbValue < iv2._dbValue;
}

bool DBEventCompare(DBEpsilonEvent e1, DBEpsilonEvent e2) {
	return e1._dbTime < e2._dbTime;
}

void Sequence2D::trendDetectStep(Group candidate, int nStep) {
	// 0.finished if nStep is the end
	if (nStep == _nLength)
	{
		candidate._nE = nStep - 1;
		addGroup(candidate);
		return;
	}

	// 1.generate segments
	vector<vector<IndexAndValue>> vecNewPointArray = generateSegments(candidate, nStep);

	// 4.handle the splited points
	if (vecNewPointArray.size() == 1)
	{
		// no splitting
		trendDetectStep(candidate, nStep + 1);
	}
	else {
		// generate a new group of the candidate
		candidate._nE = nStep - 1;
		addGroup(candidate);
		// generate new candidates
		for (size_t i = 0, length = vecNewPointArray.size(); i < length; i++)
		{
			int length2 = vecNewPointArray[i].size();
			if (length2 > _nM) {
				Group newCandidate;
				newCandidate._nS = candidate._nS;
				for (size_t j = 0; j < length2; j++)
				{
					newCandidate._member.push_back(vecNewPointArray[i][j]._nIndex);
				}
				// reset the start time of the new candidate
				newCandidate._nS = calculateStart(newCandidate, nStep);
				trendDetectStep(newCandidate, nStep + 1);
			}
		}
	}
}

void Sequence2D::trendDetectStep_SplitOnly(Group candidate, int nStep) {
	// 0.finished if nStep is the end
	if (nStep == _nLength)
	{
		candidate._nE = nStep - 1;
		addGroup(candidate);
		return;
	}

	// 1.generate segments
	vector<vector<IndexAndValue>> vecNewPointArray = generateSegments(candidate, nStep);

	// 4.handle the splited points
	if (vecNewPointArray.size() == 1)
	{
		// no splitting
		trendDetectStep_SplitOnly(candidate, nStep + 1);
	}
	else {
		// generate a new group of the candidate
		candidate._nE = nStep - 1;
		addGroup(candidate);
		// generate new candidates
		for (size_t i = 0, length = vecNewPointArray.size(); i < length; i++)
		{
			int length2 = vecNewPointArray[i].size();
			if (length2 > _nM) {
				Group newCandidate;
				newCandidate._nS = candidate._nS;
				for (size_t j = 0; j < length2; j++)
				{
					newCandidate._member.push_back(vecNewPointArray[i][j]._nIndex);
				}
				// reset the start time of the new candidate
				newCandidate._nS = calculateStart(newCandidate, nStep);

				if (newCandidate._nS==0)
					trendDetectStep_SplitOnly(newCandidate, nStep + 1);
			}
		}
	}
}

void Sequence2D::addGroup(Group candidate) {
	if (candidate._nE - candidate._nS>_nDelta && candidate._member.size()>_nM)
	{
		sort(candidate._member.begin(), candidate._member.end());
		_vecGroups.push_back(candidate);
	}
}

int Sequence2D::calculateStart(Group g, int nCurrent) {
	int nS = nCurrent - 1;
	while (nS >= g._nS)
	{
		// 1.get the point in the new timestep
		vector<IndexAndValue> vecNewPoints;
		for (size_t i = 0, length = g._member.size(); i < length; i++)
		{
			vecNewPoints.push_back(IndexAndValue(g._member[i], _vecSequences[g._member[i]][nS]));
		}

		// 2.sort them according to the value
		sort(vecNewPoints.begin(), vecNewPoints.end(), IndexAndValueCompare);


		int nCurrentStartPos = 0;	// start position of current group;
		for (size_t i = 0; i < vecNewPoints.size() - 1; i++)
		{
			if (vecNewPoints[i + 1]._dbValue - vecNewPoints[i]._dbValue>_dbEpsilon)
				return nS + 1;
		}

		nS--;
	}

	return g._nS;
}

vector<vector<IndexAndValue>> Sequence2D::generateSegments(const Group& candidate, int nStep) {
	// 1.get the point in the new timestep
	vector<IndexAndValue> vecNewPoints;
	for (size_t i = 0, length = candidate._member.size(); i < length; i++)
	{
		vecNewPoints.push_back(IndexAndValue(candidate._member[i], _vecSequences[candidate._member[i]][nStep]));
	}

	// 2.sort them according to the value
	sort(vecNewPoints.begin(), vecNewPoints.end(), IndexAndValueCompare);

	// 3.split the points in this timestep
	vector<vector<IndexAndValue>> vecNewPointArray;
	int nCurrentStartPos = 0;	// start position of current group;
	for (size_t i = 0; i < vecNewPoints.size() - 1; i++)
	{
		if (vecNewPoints[i + 1]._dbValue - vecNewPoints[i]._dbValue>_dbEpsilon) {
			// space larger than epsilon, split
			vector<IndexAndValue> newArray;
			for (size_t j = nCurrentStartPos; j <= i; j++)
			{
				newArray.push_back(vecNewPoints[j]);
			}
			vecNewPointArray.push_back(newArray);
			nCurrentStartPos = i + 1;
		}
	}
	// add the last segment
	vector<IndexAndValue> newArray;
	for (size_t j = nCurrentStartPos; j <vecNewPoints.size(); j++)
	{
		newArray.push_back(vecNewPoints[j]);
	}
	vecNewPointArray.push_back(newArray);

	return vecNewPointArray;
}

void Sequence2D::startTrendDetect(vector<IndexAndValue> segment, int nStep) {
	// check if this segment have been calculated
	unsigned long uiKey = hashFunc(segment);

	/*
	cout << nStep << " : ";
	for each (IndexAndValue iv in segment)
	{
	cout << iv._nIndex << "\t";
	}
	cout << endl;
	*/
	if (_hashTrends.find(uiKey) != _hashTrends.end())
	{
		vector<int> vecSteps = _hashTrends[uiKey];
		int nPos = 0;
		while (nPos<vecSteps.size() && vecSteps[nPos] <= nStep)
		{
			nPos++;
		}
		if (nPos % 2 == 1) {
			//	cout << "skip" << endl;
			return;
		}
	}
	else {
		vector<int> vecSteps;
		vecSteps.push_back(nStep);
		_hashTrends[uiKey] = vecSteps;
	}
	// 1.search this segment as a group
	Group g;
	g._nS = nStep;
	for each (IndexAndValue iv in segment)
		g._member.push_back(iv._nIndex);
	trendDetectStep(g, nStep + 1);
	int nLen = segment.size();
	if (nLen < _nM + 1) return;
	// 2.check if there's any small groups
	for (size_t i = 1; i < nLen - 1; i++)
	{
		if (segment[i + 1]._dbValue - segment[i - 1]._dbValue < _dbEpsilon) {
			// this points can be removed
			vector<IndexAndValue> newSeg;
			for (size_t j = 0; j < nLen; j++)
			{
				if (j != i) newSeg.push_back(segment[j]);
			}
			startTrendDetect(newSeg, nStep);
		}
	}
	// 3.remove the first or the last

	vector<IndexAndValue> newSeg1;
	vector<IndexAndValue> newSeg2;
	for (size_t i = 1; i < nLen; i++)
	{
		newSeg1.push_back(segment[i - 1]);
		newSeg2.push_back(segment[i]);
	}
	startTrendDetect(newSeg1, nStep);
	startTrendDetect(newSeg2, nStep);
}

int Sequence2D::GetGroupSize() {
	return _vecGroups.size();
}

Group Sequence2D::GetGroup(int nIndex) {
	return _vecGroups[nIndex];
}

double Sequence2D::GetValue(int nSeq, int nIndex) {
	return _vecSequences[nSeq][nIndex];
}

//===============epsilon event version
void Sequence2D::TrendDetect() {
	SortIndices();

	FindEvents();

	TrendDetectBasedOnEvents();

	// print the result
	printGroup();
}

void Sequence2D::SortIndices() {
	// for each time steps
	for (size_t i = 0; i < _nLength; i++)
	{
		// add each index and its value
		vector<IndexAndValue> vecIndices;
		for (size_t j = 0; j < _nEns; j++)
		{
			vecIndices.push_back(IndexAndValue(j, _vecSequences[j][i]));
		}

		// sort them
		sort(vecIndices.begin(), vecIndices.end(), IndexAndValueCompare);

		// get rank
		vector<int> vecRank(_nEns);
		for (size_t j = 0; j < _nEns; j++)
		{
			vecRank[vecIndices[j]._nIndex] = j;
		}

		// add to the vector
		_vectOrderedIndex.push_back(vecIndices);
		_vecRanks.push_back(vecRank);
	}
}

void Sequence2D::FindEvents() {
	// search each time step
	for (size_t iStep = 0; iStep < _nLength; iStep++)
	{
		// check each index
		for (size_t iSeq1 = 0; iSeq1 < _nEns-1; iSeq1++)
		{
			size_t iSeq2 = iSeq1 + 1;
			IndexAndValue iv1 = _vectOrderedIndex[iStep][iSeq1];
			IndexAndValue iv2 = _vectOrderedIndex[iStep][iSeq2];
			while (iv2._dbValue-iv1._dbValue<_dbEpsilon)
			{
				// while the two indices are directly epsilon-connected
				// check for start event
				if (iStep == 0 ||																				// first time step
					abs(_vecSequences[iv1._nIndex][iStep-1] - _vecSequences[iv2._nIndex][iStep-1])>_dbEpsilon)	// apart in last step
				{
					_vecStartEvents.push_back(EpsilonEvent(0, iStep, iv1._nIndex, iv2._nIndex));

				}
				// check for end event
				if (iStep == _nLength - 1 ||																		// last time step
					abs(_vecSequences[iv1._nIndex][iStep + 1] - _vecSequences[iv2._nIndex][iStep + 1])>_dbEpsilon)	// apart in next step
				{
					_vecEndEvents.push_back(EpsilonEvent(1, iStep, iv1._nIndex, iv2._nIndex));

				}
				// loop
				iSeq2++; 
				if (iSeq2 == _nEns)
				{
					break;
				}
				else {
					iv2 = _vectOrderedIndex[iStep][iSeq2];
				}
			}
		}
	}
}

void Sequence2D::TrendDetectBasedOnEvents() {
	for each (EpsilonEvent startE in _vecStartEvents)
	{
		int nTimeS = startE._nTime;
		for each(EpsilonEvent endE in _vecEndEvents) {
			// check each event pair
			int nTimeE = endE._nTime;
			if (nTimeE - nTimeS < _nDelta) continue;

			vector<vector<IndexAndValue>> vecOrderedIndex;
			for (size_t i = nTimeS; i <= nTimeE; i++)
			{
				vecOrderedIndex.push_back(_vectOrderedIndex[i]);
			}
			generateGroupFromEventPair(startE, endE,vecOrderedIndex);
		}
	}
}

void Sequence2D::generateGroupFromEventPair(EpsilonEvent eS, EpsilonEvent eE, vector<vector<IndexAndValue>> vecOrderedIndex) {
	int nTimeS = eS._nTime;									// time of start event
	int nTimeE = eE._nTime;									// time of end event
	int nIndexS1 = eS._nIndex1;
	int nIndexE1 = eE._nIndex1;
	int nIndexS2 = eS._nIndex2;
	int nIndexE2 = eE._nIndex2;

	for (size_t i = 0,length=nTimeE-nTimeS+1; i < length; i++)
	{
		// foreach time step
		for (size_t j = 0, length2 = vecOrderedIndex[i].size() - 1; j < length2; j++) {
			// check each space of two values
			if (vecOrderedIndex[i][j+1]._dbValue- vecOrderedIndex[i][j]._dbValue>_dbEpsilon)
			{
				// a face taller than epsilon has been found(i,j->j+1)

				int nIndex1 = vecOrderedIndex[i][j]._nIndex;
				int nIndex2 = vecOrderedIndex[i][j + 1]._nIndex;
				int nRank1 = _vecRanks[nTimeS + i][nIndex1];
				int nRank2 = _vecRanks[nTimeS + i][nIndex2];

				if (_vecRanks[nTimeS + i][nIndexS2] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexE2] <= nRank1)
				{
					// both event on top of the face
					for (size_t ii = 0; ii < length; ii++)
					{
						for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
						{
							// remove the members bellow
							if (_vecRanks[nTimeS + i][vecOrderedIndex[ii][jj]._nIndex]>nRank1)
								vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
						}
					}
					generateGroupFromEventPair(eS, eE, vecOrderedIndex);

				}
				else if (_vecRanks[nTimeS + i][nIndexS1] >= nRank2 &&
						 _vecRanks[nTimeS + i][nIndexE1] >= nRank2)
				{
					// both event on bottom of the face
					for (size_t ii = 0; ii < length; ii++)
					{
						for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
						{
							// remove the members above
							if (_vecRanks[nTimeS + i][vecOrderedIndex[ii][jj]._nIndex]<nRank2)
								vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
						}
					}
					generateGroupFromEventPair(eS, eE, vecOrderedIndex);
				}
				// the two events are not in the same side
				return;
			}
		}
	}
	// no face taller than epsilon
	// check if both the event is free
	int nPosS1=-1;
	int nPosE1=-1;
	int nPosS2=-1;
	int nPosE2=-1;
	int nLen = vecOrderedIndex[0].size();
	int nTimeRange = nTimeE - nTimeS;
	for (size_t i = 0; i < nLen; i++)
	{
		if (vecOrderedIndex[0][i]._nIndex == nIndexS1) nPosS1 = i;
		else if (vecOrderedIndex[0][i]._nIndex == nIndexS2) nPosS2 = i;
		if (vecOrderedIndex[nTimeRange][i]._nIndex == nIndexE1) nPosE1 = i;
		else if (vecOrderedIndex[nTimeRange][i]._nIndex == nIndexE2) nPosE2 = i;
	}
	if (nPosS2==nPosS1+1 && nPosE2==nPosE1+1)
	{
		Group newGroup;
		newGroup._nS = nTimeS;
		newGroup._nE = nTimeE;
		for each (IndexAndValue iv in vecOrderedIndex[0])
		{
			newGroup._member.push_back(iv._nIndex);
		}
		addGroup(newGroup);

	}

}


// ===========double time step
void Sequence2D::TrendDetectDB() {

	clock_t begin = clock();

	SortIndices();

	FindDBEvents();

	TrendDetectBasedOnDBEvents();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	// print the result
	printDBGroup();


	cout << "Trend Detection using time: " << elapsed_secs << endl;
}

void Sequence2D::FindDBEvents() {
	// search each time step
	for (size_t iStep = 0; iStep < _nLength; iStep++)
	{
		// check each index
		for (size_t iSeq1 = 0; iSeq1 < _nEns - 1; iSeq1++)
		{
			size_t iSeq2 = iSeq1 + 1;
			IndexAndValue iv1 = _vectOrderedIndex[iStep][iSeq1];
			IndexAndValue iv2 = _vectOrderedIndex[iStep][iSeq2];
			while (iv2._dbValue - iv1._dbValue<_dbEpsilon)
			{
				// while the two indices are directly epsilon-connected
				// check for start event
				if (iStep == 0)
				{
					// first time step
					_vecStartDBEvents.push_back(DBEpsilonEvent(0, 0, iv1._nIndex, iv2._nIndex));
				} 																			
				else if(abs(_vecSequences[iv1._nIndex][iStep - 1] - _vecSequences[iv2._nIndex][iStep - 1])>_dbEpsilon)	// apart in last step
				{
					_vecStartDBEvents.push_back(DBEpsilonEvent(0, getEventTime(iStep-1,iStep,iv1._nIndex,iv2._nIndex), iv1._nIndex, iv2._nIndex));

				}
				// check for end event
				if (iStep == _nLength - 1)
				{
					// last time step
					_vecEndDBEvents.push_back(DBEpsilonEvent(1, _nLength - 1, iv1._nIndex, iv2._nIndex));
				}
				else if(abs(_vecSequences[iv1._nIndex][iStep + 1] - _vecSequences[iv2._nIndex][iStep + 1])>_dbEpsilon)	// apart in next step
				{
					double dbTime = getEventTime(iStep, iStep + 1, iv1._nIndex, iv2._nIndex);
					_vecEndDBEvents.push_back(DBEpsilonEvent(1, dbTime, iv1._nIndex, iv2._nIndex));

				}
				// loop
				iSeq2++;
				if (iSeq2 == _nEns)
				{
					break;
				}
				else {
					iv2 = _vectOrderedIndex[iStep][iSeq2];
				}
			}
		}
	}
}

void Sequence2D::TrendDetectBasedOnDBEvents() {
	for each (DBEpsilonEvent startE in _vecStartDBEvents)
	{
		double dbTimeS = startE._dbTime;
		for each(DBEpsilonEvent endE in _vecEndDBEvents) {
			// check each event pair
			double dbTimeE = endE._dbTime;
			if (dbTimeE - dbTimeS < _nDelta) continue;
			int nTimeS = (int)(dbTimeS + .999);
			int nTimeE = (int)(dbTimeE);

			vector<vector<IndexAndValue>> vecOrderedIndex;
			for (size_t i = nTimeS; i <= nTimeE; i++)
			{
				vecOrderedIndex.push_back(_vectOrderedIndex[i]);
			}
			generateGroupFromDBEventPair(startE, endE, vecOrderedIndex);
		}
	}
}

void Sequence2D::generateGroupFromDBEventPair(DBEpsilonEvent eS, DBEpsilonEvent eE, vector<vector<IndexAndValue>> vecOrderedIndex) {
	double dbTimeS = eS._dbTime;									// time of start event
	double dbTimeE = eE._dbTime;									// time of end event
	int nTimeS = (int)(dbTimeS + .999);
	int nTimeE = (int)(dbTimeE);
	int nIndexS1 = eS._nIndex1;
	int nIndexE1 = eE._nIndex1;
	int nIndexS2 = eS._nIndex2;
	int nIndexE2 = eE._nIndex2;
	
	// 1.check there is face taller than epsilon
	for (size_t i = 0, length = nTimeE - nTimeS + 1; i < length; i++)
	{
		// foreach time step
		for (size_t j = 0, length2 = vecOrderedIndex[i].size() - 1; j < length2; j++) {
			// check each space of two values
			if (vecOrderedIndex[i][j + 1]._dbValue - vecOrderedIndex[i][j]._dbValue>_dbEpsilon)
			{
				// a face taller than epsilon has been found(i,j->j+1)

				int nIndex1 = vecOrderedIndex[i][j]._nIndex;
				int nIndex2 = vecOrderedIndex[i][j + 1]._nIndex;
				int nRank1 = _vecRanks[nTimeS + i][nIndex1];
				int nRank2 = _vecRanks[nTimeS + i][nIndex2];

				if (_vecRanks[nTimeS + i][nIndexS1] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexS2] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexE1] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexE2] <= nRank1)
				{
					// both event on top of the face
					for (size_t ii = 0; ii < length; ii++)
					{
						for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
						{
							// remove the members bellow
							if (_vecRanks[nTimeS + i][vecOrderedIndex[ii][jj]._nIndex]>nRank1)
								vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
						}
					}
					generateGroupFromDBEventPair(eS, eE, vecOrderedIndex);

				}
				else if (_vecRanks[nTimeS + i][nIndexS1] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexS2] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexE1] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexE2] >= nRank2)
				{
					// both event on bottom of the face
					for (size_t ii = 0; ii < length; ii++)
					{
						for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
						{
							// remove the members above
							if (_vecRanks[nTimeS + i][vecOrderedIndex[ii][jj]._nIndex]<nRank2)
								vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
						}
					}
					generateGroupFromDBEventPair(eS, eE, vecOrderedIndex);
				}
				// the two events are not in the same side
				return;
			}
		}
	}
	// no face taller than epsilon

	// 2.check if both the event is free
	double dbVS1 = getValue(nIndexS1, dbTimeS);
	double dbVS2 = getValue(nIndexS2, dbTimeS);
	double dbVE1 = getValue(nIndexE1, dbTimeE);
	double dbVE2 = getValue(nIndexE2, dbTimeE);

	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		if (iv._nIndex != nIndexS1 && iv._nIndex != nIndexS2)
		{
			double dbV = getValue(iv._nIndex, dbTimeS);
			if ((dbV - dbVS1)*(dbV - dbVS2) < 0) return;
		}
		if (iv._nIndex != nIndexE1 && iv._nIndex != nIndexE2)
		{
			double dbV = getValue(iv._nIndex, dbTimeE);
			if ((dbV - dbVE1)*(dbV - dbVE2) < 0) return;
		}
	}
	// both the event is free

	// 3.add this group
	DBGroup newGroup;
	newGroup._dbS = dbTimeS;
	newGroup._dbE = dbTimeE;
	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		newGroup._member.push_back(iv._nIndex);
	}
	addDBGroup(newGroup);
}

double Sequence2D::getEventTime(int nTime1, int nTime2, int nIndex1, int nIndex2) {
	double dbV11 = _vecSequences[nIndex1][nTime1];
	double dbV12 = _vecSequences[nIndex1][nTime2];
	double dbV21 = _vecSequences[nIndex2][nTime1];
	double dbV22 = _vecSequences[nIndex2][nTime2];

	// at least one end point of seg1 is less than seg2
	if (dbV21 < dbV11 || dbV22 < dbV12){
		return nTime1 + (_dbEpsilon - dbV11 + dbV21) / (dbV12 - dbV11 - dbV22 + dbV21);
	}
	else{
		return nTime1 + (-_dbEpsilon - dbV11 + dbV21) / (dbV12 - dbV11 - dbV22 + dbV21);
	}
	
}

double Sequence2D::getValue(int nIndex, double dbTime) {
	int nTime1 = (int)dbTime;
	int nTime2 = nTime1 + 1;
	double dbV1 = _vecSequences[nIndex][nTime1];
	if (nTime1 == _nLength - 1) return dbV1;
	double dbV2 = _vecSequences[nIndex][nTime2];
	return dbV1 + (dbTime - nTime1)*(dbV2 - dbV1);
}

void Sequence2D::addDBGroup(DBGroup candidate) {
	sort(candidate._member.begin(), candidate._member.end());
	// check if this candidates has been added
	bool bDuplicate = false;
	if (candidate._dbS<.001)
	{
		for each (DBGroup g in _vecDBGroups)
		{
			if (abs(candidate._dbE - g._dbE) < .001 && candidate._member.size()==g._member.size()) {
				bDuplicate = true;
				for (size_t i = 0,length=g._member.size(); i < length; i++)
				{
					if (candidate._member[i]!=g._member[i])
					{
						bDuplicate = false;
						break;
					}
				}
			}
		}
	}
	if (bDuplicate) return;
	if (candidate._dbE - candidate._dbS>_nDelta && candidate._member.size()>_nM)
	{
		_vecDBGroups.push_back(candidate);
	}
}

int Sequence2D::GetDBGroupSize() {
	return _vecDBGroups.size();
}

DBGroup Sequence2D::GetDBGroup(int nIndex) {
	return _vecDBGroups[nIndex];
}

void Sequence2D::printDBGroup() {
	for each (DBGroup g in _vecDBGroups)
	{
#ifdef REVERSE
		cout << _nWidth - 1 - g._nE << "\t" << _nWidth - 1 - g._nS << "\t";
#else
		cout << "(" << g._dbS << "," << g._dbE << "):";
#endif		
		for each (int nMember in g._member)
		{
			cout << nMember << ",";
		}
		cout << endl;
	}
	cout << "number of groups: " << _vecDBGroups.size() << endl;
}

double Sequence2D::GetDBValue(double dbTime, int nIndex) {
	return getValue(nIndex, dbTime);
}

// ===========improved algorithm
void Sequence2D::TrendDetectDBImproved() {

	clock_t begin = clock();


	SortIndices();

	FindDBEvents();

	TrendDetectBasedOnDBEventsImproved();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	// print the result
	printDBGroup();
	cout << "Trend Detection using time: " << elapsed_secs << endl;
}

void Sequence2D::TrendDetectBasedOnDBEventsImproved() {
	// sort end events
	sort(_vecEndDBEvents.begin(), _vecEndDBEvents.end(), DBEventCompare);
	sort(_vecStartDBEvents.begin(), _vecStartDBEvents.end(), DBEventCompare);

	// search for groups
	for each (DBEpsilonEvent startE in _vecStartDBEvents)
	{
		double dbTimeS = startE._dbTime;
		int nTimeS = (int)(dbTimeS + .999);		

		// build the arrangement from the time of the start event
		vector<vector<IndexAndValue>> vecOrderedIndex;
		for (size_t i = nTimeS; i <_nLength; i++)
		{
			vecOrderedIndex.push_back(_vectOrderedIndex[i]);
		}

		for each(DBEpsilonEvent endE in _vecEndDBEvents) {
			// check each event pair
			double dbTimeE = endE._dbTime;
			if (dbTimeE - dbTimeS < _nDelta) continue;
			int nTimeE = (int)(dbTimeE);


			generateGroupFromDBEventPairImproved(startE, endE, vecOrderedIndex);
			if (vecOrderedIndex[0].size() < _nM) break;
		}
	}
}

void Sequence2D::generateGroupFromDBEventPairImproved(DBEpsilonEvent eS, DBEpsilonEvent eE, vector<vector<IndexAndValue>>& vecOrderedIndex) {
//	if (abs(eS._dbTime - 84.6999) < .001&&abs(eE._dbTime - 90) < .001) {
//		cout << "hehe";
//	}
	double dbTimeS = eS._dbTime;									// time of start event
	double dbTimeE = eE._dbTime;									// time of end event
	int nTimeS = (int)(dbTimeS + .999);
	int nTimeE = (int)(dbTimeE);
	int nIndexS1 = eS._nIndex1;
	int nIndexE1 = eE._nIndex1;
	int nIndexS2 = eS._nIndex2;
	int nIndexE2 = eE._nIndex2;
	int nLenArrangement = vecOrderedIndex.size();
	// 0.check if the two members of the end event are still in vecOrderedIndex
	int nCount = 0;
	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		if (iv._nIndex == nIndexE1) nCount++;
		else if (iv._nIndex == nIndexE2) nCount++;
	}
	if (nCount < 2) return;

	// check if the members are epsilon collected at the time of the start event
	vector<IndexAndValue> vecStartTimeArrange;
	vector<int> vecStartTimeRank(_nEns);
	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		vecStartTimeArrange.push_back(IndexAndValue(iv._nIndex, getValue(iv._nIndex, dbTimeS)));
	}
	sort(vecStartTimeArrange.begin(), vecStartTimeArrange.end(), IndexAndValueCompare);
	for (size_t i = 0,length=vecStartTimeArrange.size(); i < length; i++)
	{
		vecStartTimeRank[vecStartTimeArrange[i]._nIndex] = i;
	}

	for (size_t j = 0, length2 = vecStartTimeArrange.size() - 1; j < length2; j++) {
		// check each space of two values
		if (vecStartTimeArrange[j + 1]._dbValue - vecStartTimeArrange[j]._dbValue>_dbEpsilon+.001)
		{
			// a face taller than epsilon has been found(i,j->j+1)

			int nIndex1 = vecStartTimeArrange[j]._nIndex;
			int nIndex2 = vecStartTimeArrange[j + 1]._nIndex;
			int nRank1 = vecStartTimeRank[nIndex1];
			int nRank2 = vecStartTimeRank[nIndex2];

			if (vecStartTimeRank[nIndexS1] <= nRank1 &&
				vecStartTimeRank[nIndexS2] <= nRank1 &&
				vecStartTimeRank[nIndexE1] <= nRank1 &&
				vecStartTimeRank[nIndexE2] <= nRank1)
			{
				// both event on top of the face
				splitArrangement(vecStartTimeRank, nRank1, vecOrderedIndex, true);
				generateGroupFromDBEventPairImproved(eS, eE, vecOrderedIndex);

			}
			else if (vecStartTimeRank[nIndexS1] >= nRank2 &&
				vecStartTimeRank[nIndexS2] >= nRank2 &&
				vecStartTimeRank[nIndexE1] >= nRank2 &&
				vecStartTimeRank[nIndexE2] >= nRank2)
			{
				// both event on bottom of the face
				splitArrangement(vecStartTimeRank, nRank2, vecOrderedIndex, false);
				generateGroupFromDBEventPairImproved(eS, eE, vecOrderedIndex);
			}
			// the two events are not in the same side
			return;
		}
	}

	// 1.check there is face taller than epsilon
	for (size_t i = 0, length = nTimeE - nTimeS + 1; i < length; i++)
	{
		// foreach time step
		for (size_t j = 0, length2 = vecOrderedIndex[i].size() - 1; j < length2; j++) {
			// check each space of two values
			if (vecOrderedIndex[i][j + 1]._dbValue - vecOrderedIndex[i][j]._dbValue>_dbEpsilon)
			{
				// a face taller than epsilon has been found(i,j->j+1)

				int nIndex1 = vecOrderedIndex[i][j]._nIndex;
				int nIndex2 = vecOrderedIndex[i][j + 1]._nIndex;
				int nRank1 = _vecRanks[nTimeS + i][nIndex1];
				int nRank2 = _vecRanks[nTimeS + i][nIndex2];

				if (_vecRanks[nTimeS + i][nIndexS1] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexS2] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexE1] <= nRank1 &&
					_vecRanks[nTimeS + i][nIndexE2] <= nRank1)
				{
					// both event on top of the face
					splitArrangement(nTimeS + i, nRank1, vecOrderedIndex, true);
					generateGroupFromDBEventPairImproved(eS, eE, vecOrderedIndex);

				}
				else if (_vecRanks[nTimeS + i][nIndexS1] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexS2] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexE1] >= nRank2 &&
					_vecRanks[nTimeS + i][nIndexE2] >= nRank2)
				{
					// both event on bottom of the face
					splitArrangement(nTimeS + i, nRank2, vecOrderedIndex, false);
					generateGroupFromDBEventPairImproved(eS, eE, vecOrderedIndex);
				}
				// the two events are not in the same side
				return;
			}
		}
	}
	// no face taller than epsilon

	// 2.check if both the event is free
	double dbVS1 = getValue(nIndexS1, dbTimeS);
	double dbVS2 = getValue(nIndexS2, dbTimeS);
	double dbVE1 = getValue(nIndexE1, dbTimeE);
	double dbVE2 = getValue(nIndexE2, dbTimeE);

	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		if (iv._nIndex != nIndexS1 && iv._nIndex != nIndexS2)
		{
			double dbV = getValue(iv._nIndex, dbTimeS);
			if ((dbV - dbVS1)*(dbV - dbVS2) < 0) return;
		}
		if (iv._nIndex != nIndexE1 && iv._nIndex != nIndexE2)
		{
			double dbV = getValue(iv._nIndex, dbTimeE);
			if ((dbV - dbVE1)*(dbV - dbVE2) < 0) return;
		}
	}
	// both the event is free

	// 3.add this group
	DBGroup newGroup;
	newGroup._dbS = dbTimeS;
	newGroup._dbE = dbTimeE;
	for each (IndexAndValue iv in vecOrderedIndex[0])
	{
		newGroup._member.push_back(iv._nIndex);
	}
	addDBGroup(newGroup);

	// 4 split by this event
	if (nTimeE < _nLength - 1) 
	{
		nTimeE++;		// if nTimeE is not the last timestep, use next timestep
		int nRank1 = _vecRanks[nTimeE][nIndexE1];
		int nRank2 = _vecRanks[nTimeE][nIndexE2];
		if (nRank1>nRank2)
		{
			int nTemp = nRank1;
			nRank1 = nRank2;
			nRank2 = nTemp;
		}

		if (_vecRanks[nTimeE][nIndexS1] <= nRank1 &&
			_vecRanks[nTimeE][nIndexS2] <= nRank1)
		{
			// both event on top of the face
			splitArrangement(nTimeE, nRank1, vecOrderedIndex, true);
		}
		else if (_vecRanks[nTimeE][nIndexS1] >= nRank2 &&
			_vecRanks[nTimeE][nIndexS2] >= nRank2)
		{
			// both event on bottom of the face
			splitArrangement(nTimeE, nRank2, vecOrderedIndex, false);
		}
		else {
			for (size_t i = 0, length = vecOrderedIndex.size(); i < length; i++)
			{
				vecOrderedIndex[i].clear();
			}
		}
	}
	else {
		for (size_t i = 0, length = vecOrderedIndex.size(); i < length; i++)
		{
			vecOrderedIndex[i].clear();
		}
	}
}

void Sequence2D::splitArrangement(int nTime,int nRank, std::vector<std::vector<IndexAndValue>>& vecOrderedIndex, bool bTop) {

	for (size_t ii = 0,length= vecOrderedIndex.size(); ii < length; ii++)
	{
		for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
		{
			// remove the members bellow
			if ((_vecRanks[nTime][vecOrderedIndex[ii][jj]._nIndex]>nRank&&bTop)
				|| (_vecRanks[nTime][vecOrderedIndex[ii][jj]._nIndex]<nRank&&!bTop))
				vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
		}
	}
}

void Sequence2D::splitArrangement(vector<int> vecRank, int nRank, std::vector<std::vector<IndexAndValue>>& vecOrderedIndex, bool bTop) {

	for (size_t ii = 0, length = vecOrderedIndex.size(); ii < length; ii++)
	{
		for (int jj = vecOrderedIndex[ii].size() - 1; jj >= 0; jj--)
		{
			// remove the members bellow
			if ((vecRank[vecOrderedIndex[ii][jj]._nIndex]>nRank&&bTop)
				|| (vecRank[vecOrderedIndex[ii][jj]._nIndex]<nRank&&!bTop))
				vecOrderedIndex[ii].erase(vecOrderedIndex[ii].begin() + jj);
		}
	}
}