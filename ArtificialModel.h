#pragma once
#include "MeteModel.h"
class ArtificialModel :
	public MeteModel
{
public:
	ArtificialModel();
	~ArtificialModel();
protected:
	// specialized model initialization
	virtual void initializeModel();
private:
	void regenerateData();		// regenerate data from contours
};

