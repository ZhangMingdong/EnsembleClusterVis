#pragma once
#include "MeteModel.h"

/*
	model of the ensemble data
	Mingdong
	2017/07/17
*/
class EnsembleModel :
	public MeteModel
{
public:
	EnsembleModel();
	virtual ~EnsembleModel();
protected:
	// projected 2D points of the data
	std::vector<DPoint3> _vecPoints;
	// matrix of the data. d1:day,d2: step
	std::vector<std::vector<DataField* >> _matrixData;
	int _nDay = 31;
	int _nStep = 61;
//	int _nWidth = 31;
//	int _nHeight = 31;
	int _nEnsembleLen = 50;
protected:

	// read ensemble data from text file
	virtual void readDataFromText();
public:
	virtual GLubyte* generateTextureRange(int nIndex);
};

