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
	QList<QList<ContourLine>> regenerateData();		// regenerate data from contours
	void generateDataFromContour();					// regenerate data from given contours
	/*
		regenerate data from given field
		pattern1:
			5 high value and 5 low value plus 1 outlier control
	*/
	void generateDataFromField();				// regenerate data from given contours
	/*
		regenerate data from given field
		pattern2:
			one core
	*/
	void generateDataFromField_2();
	/*
		modified from generateDataFromField, just two high core and two low core
	*/
	void generateDataFromField_3();

	/*
		from 3
		reverse the modes
	*/
	void generateDataFromField_4();
};

