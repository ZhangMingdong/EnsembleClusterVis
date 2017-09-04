/*
	used to generate contour strip
	Mingdong
	2017/03/28
*/
#pragma once

#include <QQueue>
#include "BasicStruct.h"
#include "UnCertaintyArea.h"
class ContourStripGenerator
{
public:
	ContourStripGenerator();
	~ContourStripGenerator();
public:
	// generate the uncertainty areas
	void Generate(QList<UnCertaintyArea*>& areas
		,const QList<ContourLine>& contourMin
		, const QList<ContourLine>& contourMax
		, int nWidth
		, int nHeight
		, int nFocusX
		, int nFocusY
		, int nFocusW
		, int nFocusH
		);
private:
	// generate uncertainty areas from the unclosed contours
	void generateFromUnClosedContours(QList<UnCertaintyArea*>& areas);
	// insert the closed contours as embedded areas
	void insertClosedContours(QList<UnCertaintyArea*>& areas);
	// generate embedded areas from the built structure
	void generateEmbeddedArea(QList<UnCertaintyArea*>& areas);
	// generate the list of endpoints
	void generateEndPointsList();
	// handle an area in queue, if queue is empty, return false;
	bool handleAnAreaInQueue(QQueue<UnCertaintyArea*>& queueAreas, QList<UnCertaintyArea*>& areas);
	// generate the contour tree
	ContourTreeNode generateContourTree();

	/*
		generate contour areas 
		params:
			root: the current root of the tree
			areas: store the generated area list
			nParentState: the state of the parent:0:lower;1:uncertainty;2:higher;-1:not clear
			ptParent: one of the point of its parent, used for the leaf to decide which side is its area
			bRoot: if this is the root of the tree
			return the state of the generated area
	*/
	int generateArea(const ContourTreeNode& root, QList<UnCertaintyArea*>& areas, UnCertaintyArea* pParent=NULL);

	// reset the state of the segmentations
	void resetState(QList<UnCertaintyArea*>& areas);

	// build the contour tree
	ContourTreeNode buildTree(QList<ContourTreeNode>& nodes);

	/*
		generate a closed contour from a point list, find the corresponding border
		params:
			contour: store the result contour
			pts: the given point list
			ptRef: the reference point, which is outside this area and on the border.
	*/
	void generateLeafArea(ClosedContour& contour, const QList<QPointF>& pts, QPointF ptRef);

	/*
		generate the uncertainty area which is not leaf
	*/
	void generateNoLeafArea(ClosedContour& contour, const ContourTreeNode& node);


	// check which border of the point on: 0,1,2,3
	int CheckBorderPoint(const QPointF& pt);

	// calculate the relative distance of the point from the corner point 0
	double PointDistance(const QPointF& pt);

	// return the point in front using our topology order,0 or 1
	int PointCompare(const QPointF& pt0, const QPointF& pt1);
private:
	QList<ContourLine> _listClosed;				// list of closed contours
	QList<ContourLine> _listUnClosed;			// list of unclosed contours
	QList<int> _listEndPointsContourIndex;		// index of contours
	QList<int> _listEndPointsIndex;				// index of endpoints
	QList<QPointF> _listEndPoints;				// list of endpoints

	int _nWidth;
	int _nHeight;
	int _nFocusX;
	int _nFocusY;
	int _nFocusW;
	int _nFocusH;


	// list of corner points
	QPointF arrCorners[4];

	// represent an end point
	struct EndPoint{
		int _nContourIndex;		// index of the contour it belong to
		int _nIndex;			// index of the EndPoint: 0:first, 1:last
		double _fDis;			// distance from the base corner
		bool operator<(const EndPoint& pt) const{
			return _fDis < pt._fDis;
		}
	};
};