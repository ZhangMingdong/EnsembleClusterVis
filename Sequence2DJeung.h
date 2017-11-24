#pragma once
#include "Sequence2D.h"
#include <set>
#include <vector>

class Sequence2DJeung :
	public Sequence2D
{
public:
	Sequence2DJeung();
	virtual ~Sequence2DJeung();
public:
	// detect the trend
	virtual void TrendDetectDB();
	// detect the trend
	virtual void TrendDetect();
	// get the event list
	virtual std::vector<DBEpsilonEvent> GetEvents();
protected:
	// group at each time step
	void groupTimeSteps();
	// detect trends based on the groups
	void trendDetectBasedOnGroups();
protected:
	// cluster at each timestep
	std::vector<std::vector<std::set<int>>> _vecClusters;
};

