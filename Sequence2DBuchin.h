#pragma once
#include "Sequence2D.h"
class Sequence2DBuchin :
	public Sequence2D
{
public:
	Sequence2DBuchin();
	virtual ~Sequence2DBuchin();
public:
	// detect the trend
	virtual void TrendDetect();
	virtual void TrendDetectDB();
	// get the event list
	virtual std::vector<DBEpsilonEvent> GetEvents() {
		return _vecFreeDBEvents;
	};
private:
	// find all the dbEvents
	// new version, check each pair of event
	void findDBEvents();
	// old version, find event for memebers in the same group at a given time, check whether there are near enough in the adjacent time
	void findDBEvents_1();
	// check if nIndex1 and nIndex2 is free at time dbTime
	// 1:free event, index1 is less; -1: free event, index2 is less; 0:not free event
	int isDBEventFree(double dbTime, int nIndex1, int nIndex2);
	// detect groups from Events
	void detectGroupsFromDBEvents();
private:
	std::vector<DBEpsilonEvent> _vecFreeDBEvents;
};

