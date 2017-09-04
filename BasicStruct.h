#pragma once
#include <QList>
#include <QPointF>
#include <QGLWidget>
/*
	Basic struct used in this project
*/



// a contour line
struct ContourLine
{
	int _nState;			// -1:min;1:max, used in area calculating
	int _nClosed;			// 0:unclosed;1:closed
	QList<QPointF> _listPt;	// point list
};
// the line segment used in marching squares
struct LineSeg
{
	// i11,j11,i12,j12,i21,j21,i22,j22
	int _arrIndices[8];
	LineSeg(int i11, int j11, int i12, int j12, int i21, int j21, int i22, int j22);
};

// compare the point in the borders
int PointCompare(const QPointF& pt0, const QPointF& pt1);

// a node in the contour tree
struct ContourTreeNode 
{
	// list of its children
	QList<ContourTreeNode> _children;
	// index in the contour list
	int _nContourIndex;
	// 0: not reverse;1:reverse
	int _nReverse;
	// 0:unmatched;1:matched
	int _nMatched;
};
