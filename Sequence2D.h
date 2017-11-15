#pragma once
#include <unordered_map>
#include <vector>

//#define REVERSE

// an epsilon event
struct EpsilonEvent {
	int _nType = 0;		// 0:merge,1:split
	int _nTime = 0;
	int _nIndex1 = 0;
	int _nIndex2 = 0;
	EpsilonEvent(int nType, int nTime, int nIndex1, int nIndex2):_nType(nType),_nTime(nTime),_nIndex1(nIndex1),_nIndex2(nIndex2) {

	}
};


struct Group {
	std::vector<int> _member;	// members
	int _nS;					// start time
	int _nE;					// end time
								// check whether containing g
	bool Contain(const Group& g) {
		int nLen = _member.size();
		if (g._member.size() == nLen)
		{
			for (size_t i = 0; i < nLen; i++)
				if (g._member[i] != _member[i]) return false;
			return _nS <= g._nS&&_nE >= g._nE;
		}
		return false;
	}
};
// pair of index and value
struct IndexAndValue {
	int _nIndex = 0;
	double _dbValue = 0;
	IndexAndValue(int nIndex, double dbValue) :_nIndex(nIndex), _dbValue(dbValue) {};
};

/*
	class used to represent the 2d sequences and detect trends
	mingdong
	2017/11/14
*/
class Sequence2D
{
public:
	Sequence2D();
	~Sequence2D();
protected:
	// the sequences, same length
	std::vector<std::vector<double>> _vecSequences;
	int _nLength = 0;
	int _nEns = 0;
	double _dbMin = 0;
	double _dbMax = 0;
	// groups
	std::vector<Group> _vecGroups;
	// trend detecting parameters
	int _nM = 2;
	double _dbEpsilon = 0;
	int _nDelta = 5;

	// hashtable of calculated trends
	std::unordered_map<unsigned long, std::vector<int> > _hashTrends;
public:
	// add sequence
	void AddSequence(std::vector<double> seq);
	// initialization
	void Init(int nLength, double dbMin, double dbMax, double dbEpsilon);
	int GetGroupSize();
	Group GetGroup(int nIndex);
	double GetValue(int nSeq, int nIndex);
	// detect the trend from the root
	void TrendDetectRootOnly();
	// detect the trend
	void TrendDetect();
	// detect the trend
	// old version, too slow
	// 2017/11/13
	void TrendDetect_2();
	// detect the trend
	// old version
	// 2017/11/12
	void TrendDetect_1();
	double GetMin() { return _dbMin; }
	double GetMax() { return _dbMax; }
	int GetEns() { return _nEns; }
	double GetEpsilon() { return _dbEpsilon; }
	int GetLength() { return _nLength; }
private:
	// detect the trend at a single step
	// recursive version, only find the trend start from beginning
	// only split, no merge
	void trendDetectStep_SplitOnly(Group candidate, int nStep);

	// detect the trend at a single step
	// recursive version, only find the trend start from beginning
	void trendDetectStep(Group candidate, int nStep);

	// start detecting trend from a segment at a given time step
	void startTrendDetect(std::vector<IndexAndValue> segment, int nStep);

	// generate segments given a candidate and a timestep
	std::vector<std::vector<IndexAndValue>> generateSegments(const Group& candidate, int nStep);

	// add a new group
	void addGroup(Group g);

	// calculate the start time of g
	// which whould between g.start and nCurrent
	int calculateStart(Group g, int nCurrent);

	// print the calculated groups
	void printGroup();

//===============epsilon event version
private:
	std::vector<EpsilonEvent> _vecStartEvents;
	std::vector<EpsilonEvent> _vecEndEvents;
	// ordered index sequence of each timestep
	std::vector<std::vector<IndexAndValue>> _vectOrderedIndex;	
	// in each time step, the ranks of each index in the ordered list
	std::vector<std::vector<int>> _vecRanks;
private:
	// sort the indices
	void SortIndices();
	// find all the events
	void FindEvents();
	// detect trend based on the found events
	void TrendDetectBasedOnEvents();
	// generate group from a event pair
	void generateGroupFromEventPair(EpsilonEvent eS, EpsilonEvent eE,std::vector<std::vector<IndexAndValue>> vecOrderedIndex);

};

