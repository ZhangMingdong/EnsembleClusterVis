#pragma once
#include <unordered_map>
#include <vector>
#include <set>

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
	std::set<int> _member;	// members
	int _nS;					// start time
	int _nE;					// end time
	bool _bAssigned = false;
								// check whether containing g
	bool Contain(const Group& g) {
		return _member == g._member;
	}
};

struct DBEpsilonEvent {
	int _nType = 0;		// 0:merge,1:split
	double _dbTime = 0;
	int _nIndex1 = 0;
	int _nIndex2 = 0;
	DBEpsilonEvent(int nType, double dbTime, int nIndex1, int nIndex2) :_nType(nType), _dbTime(dbTime), _nIndex1(nIndex1), _nIndex2(nIndex2) {

	}
};

struct DBGroup {
	std::set<int> _member;	// members
	double _dbS;				// start time
	double _dbE;				// end time
	double _dbR;				// radius of the group
};

// pair of index and value
struct IndexAndValue {
	int _nIndex = 0;
	double _dbValue = 0;
	IndexAndValue(int nIndex, double dbValue) :_nIndex(nIndex), _dbValue(dbValue) {};
};
bool IndexAndValueCompare(IndexAndValue iv1, IndexAndValue iv2);

bool DBEventCompare(DBEpsilonEvent e1, DBEpsilonEvent e2);
/*
	class used to represent the 2d sequences and detect trends
	mingdong
	2017/11/14
*/
class Sequence2D
{
public:
	Sequence2D();
	virtual ~Sequence2D();
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
	double _dbEpsilon = 0;
//	int _nM = 2;
//	int _nDelta = 5;

//	int _nM = 1;
//	int _nDelta = 1;

	int _nM = 3;
	int _nDelta = 3;

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
	double GetDBValue(double dbTime, int nIndex);
	// detect the trend from the root
	void TrendDetectRootOnly();
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
protected:
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
protected:
	std::vector<EpsilonEvent> _vecStartEvents;
	std::vector<EpsilonEvent> _vecEndEvents;
	// ordered index sequence of each timestep
	std::vector<std::vector<IndexAndValue>> _vectOrderedIndex;	
	// in each time step, the ranks of each index in the ordered list
	std::vector<std::vector<int>> _vecRanks;
protected:
	// find all the events
	void FindEvents();
	// detect trend based on the found events
	void TrendDetectBasedOnEvents();
	// generate group from a event pair
	void generateGroupFromEventPair(EpsilonEvent eS, EpsilonEvent eE,std::vector<std::vector<IndexAndValue>> vecOrderedIndex);

// ===========double time step

protected:
	std::vector<DBGroup> _vecDBGroups;
	// generate group from a event pair
	void generateGroupFromDBEventPair(DBEpsilonEvent eS, DBEpsilonEvent eE, std::vector<std::vector<IndexAndValue>> vecOrderedIndex);

	// get the time of epsilon
	double getEventTime(int nTime1, int nTime2, int nIndex1, int nIndex2);

	// get the two time for cross event, ordered
	void getCrossEventTimes(double& dbTime1,double& dbTime2,int nTime1, int nTime2, int nIndex1, int nIndex2);

	// add a new group
	void addDBGroup(DBGroup g);
	// print the calculated groups
	void printDBGroup();

// ===========improved algorithm
protected:

	// generate group from a event pair
	void generateGroupFromDBEventPairImproved(DBEpsilonEvent eS, DBEpsilonEvent eE, std::vector<std::vector<IndexAndValue>>& vecOrderedIndex);

	// split the arrange ment
	// bTop==true: keep the members on above
	// bTop==false: keep the members on bellow
	void splitArrangement(int nTime,int nRank, std::vector<std::vector<IndexAndValue>>& vecOrderedIndex, bool bTop);
	void splitArrangement(std::vector<int> vecRank, int nRank, std::vector<std::vector<IndexAndValue>>& vecOrderedIndex, bool bTop);
public:
	enum SequenceType
	{
		ST_Buchin,
		ST_Van,
		ST_Jeung
	};
	static Sequence2D* GenerateInstance(SequenceType type);
	//=================used in both class=================
protected:
	// sort the indices
	void sortIndices();
	// get value of nIndex at time dbTime
	double getValue(int nIndex, double dbTime);
public:
	// detect the trend
	virtual void TrendDetectDB()=0;
	// detect the trend
	virtual void TrendDetect();
	// get the event list
	virtual std::vector<DBEpsilonEvent> GetEvents()=0;
public:
	int GetDBGroupSize();
	DBGroup GetDBGroup(int nIndex);
};

