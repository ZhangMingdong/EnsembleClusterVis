#pragma once
#include "MeteModel.h"
/*
	implementation of contour band depth, used for contour boxplots visualization
	Mingdong
	2017/06/26
*/
class ContourBandDepthModel :
	public MeteModel
{
public:
	ContourBandDepthModel();
	virtual ~ContourBandDepthModel();
protected:
	virtual void initializeModel();				// specialized model initialization
private:
	std::vector<int> _vecContourBandDepth;		// contour band depth of each ensemble member
private:
	void calculateCBD(double dbIso);			// calculate contour band depth
	bool testDepth(double dbIso
		, const double* fieldI
		, const double* fieldJ
		, const double* fieldK);				// test whether i in j and k

};

