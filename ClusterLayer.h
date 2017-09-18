#pragma once
#include "MeteLayer.h"

/*
	layer used to draw the cluster chart
	mingdong
	2017/09/15
*/
class ClusterLayer :
	public MeteLayer
{
public:
	ClusterLayer();
	virtual ~ClusterLayer(); 
public:

	virtual void SetUncertaintyAreas(int nUCAreas) { _nUncertaintyRegions = nUCAreas; };

	virtual void SetFocusedCluster(int nFocusedCluster) { _nFocusedCluster = nFocusedCluster; };
protected:
	virtual void draw(DisplayStates states);
	virtual void init();

private:
	// draw the chart of variance
	void drawVarChart();

	// draw the pca points
	void drawPCAPoints();

	// draw the cluster bars
	void drawClusterBars();

	// draw the relationships of uncertainty areas
	void drawUCAreaRelation();
private:
	// number of uncertainty areas to show
	int _nUncertaintyRegions = 4;
	int _nFocusedCluster = 0;
};

