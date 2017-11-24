#pragma once
#include "Sequence2D.h"
class Sequence2DVan :
	public Sequence2D
{
public:
	Sequence2DVan();
	virtual ~Sequence2DVan();
protected:
	std::vector<DBEpsilonEvent> _vecStartDBEvents;
	std::vector<DBEpsilonEvent> _vecEndDBEvents;
protected:
	// detect trend based on the found events
	void trendDetectBasedOnDBEvents();
	// detect trend based on the found events
	void trendDetectBasedOnDBEventsImproved();
	// find all the events
	void findDBEvents();
public:
	// detect the trend
	virtual void TrendDetectDB();
	// get the event list
	virtual std::vector<DBEpsilonEvent> GetEvents() {
		return _vecStartDBEvents;
	};
};

