#pragma once
#include <QList>
#include <QPointF>
#include "BasicStruct.h"


// struct to represent a closed contour
struct ClosedContour
{
	// point list
	QList<QPointF> _listPt;
	double _fX0;
	double _fX1;
	double _fY0;
	double _fY1;
	void GenerateBounding();									// generate the bounding box
	bool Contain(const ClosedContour& contour);					// check if contain contour
	bool rayCasting_on(const ClosedContour& contour);			// ray casting algorithm based on on the gridline
	bool rayCasting_intersection(const ClosedContour& contour);	// ray casting algorithm based on intersection
};

// struct to represent the uncertainty area
class UnCertaintyArea
{
public:
	UnCertaintyArea();
	~UnCertaintyArea();

	// closed contour of this area
	ClosedContour _contour;
	// list of embedded area
	QList<UnCertaintyArea*> _listEmbeddedArea;
	/*
		state of this area: 
			-1:lower;
			0:uncertain;
			1:higher;
			-2:with max border
			2:with min border
	*/
	int _nState;
	// the list of the children area
	QList<UnCertaintyArea*> _listChildren;
	// the parent area
	UnCertaintyArea* _pParent;
public:
	// check if this area is embedded in this area. Insert it and return 1. Otherwise return 0
	int InsertContour(UnCertaintyArea* area);
	// return the size of all the contours
	int ContourSize();
	// release resource
	void Release();
};

