#pragma once
#include "MeteModel.h"
/*
	model of the reanalysis data
	Mingdong
	2017/05/05
*/
class ReanalysisModel :
	public MeteModel
{
public:
	ReanalysisModel();
	virtual ~ReanalysisModel();
public:
	// initialize the model
	virtual void InitModel(int nEnsembleLen, int nWidth, int nHeight, int nFocusX, int nFocusY, int nFocusW, int nFocusH
		, QString strFile, bool bBinary = false
		, int nWest = -179, int nEast = 180, int nSouth = -90, int nNorth = 90
		, int nFocusWest = -179, int nFocusEast = 180, int nFocusSouth = -90, int nFocusNorth = 90, bool bFiltered = false);
protected:

	// read ensemble data from text file
	virtual void readDataFromText();
	// calculate relationship between adjacent gridpoint
	virtual void calculateRelation();
protected:
	// relationship between adjacent grid points
	DataField* _pRelationH = NULL;
	DataField* _pRelationV = NULL;
};

