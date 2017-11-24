#include "Sequence2DVan.h"

#include <ctime>
#include <algorithm>
#include <iostream>
using namespace std;

Sequence2DVan::Sequence2DVan()
{
}


Sequence2DVan::~Sequence2DVan()
{
}

void Sequence2DVan::TrendDetectDB() {

	clock_t begin = clock();

	sortIndices();

	findDBEvents();

//	trendDetectBasedOnDBEvents();
	trendDetectBasedOnDBEventsImproved();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	// print the result
	printDBGroup();


	cout << "Trend Detection using time: " << elapsed_secs << endl;

}

void Sequence2DVan::trendDetectBasedOnDBEvents() {
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

void Sequence2DVan::trendDetectBasedOnDBEventsImproved() {
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

void Sequence2DVan::findDBEvents() {
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
				else if (abs(_vecSequences[iv1._nIndex][iStep - 1] - _vecSequences[iv2._nIndex][iStep - 1])>_dbEpsilon)	// apart in last step
				{
					_vecStartDBEvents.push_back(DBEpsilonEvent(0, getEventTime(iStep - 1, iStep, iv1._nIndex, iv2._nIndex), iv1._nIndex, iv2._nIndex));

				}
				// check for end event
				if (iStep == _nLength - 1)
				{
					// last time step
					_vecEndDBEvents.push_back(DBEpsilonEvent(1, _nLength - 1, iv1._nIndex, iv2._nIndex));
				}
				else if (abs(_vecSequences[iv1._nIndex][iStep + 1] - _vecSequences[iv2._nIndex][iStep + 1])>_dbEpsilon)	// apart in next step
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