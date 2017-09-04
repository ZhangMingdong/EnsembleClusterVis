#include "ContourStripGenerator.h"
#include <QDebug>
#include "def.h"



// check which border of the point on: 0,1,2,3
int ContourStripGenerator::CheckBorderPoint(const QPointF& pt){
	double fX = pt.x();
	double fY = pt.y();
	if (abs(fX - _nFocusX) < g_fDifference) return 0;
	else if (abs(fY - _nFocusY) < g_fDifference) return 1;
	else if (abs(fX - _nFocusX - _nFocusW + 1) < g_fDifference) return 2;
	else if (abs(fY - _nFocusY - _nFocusH + 1) < g_fDifference) return 3;
	else{		
		qDebug() << "error";
		return -1;
	}
}

// calculate the relative distance of the point from the corner point 0
double ContourStripGenerator::PointDistance(const QPointF& pt){
	double fDis = 0;
	double fX = pt.x() - _nFocusX;
	double fY = pt.y() - _nFocusY;
	double fBiasX = fX / (_nFocusW - 1);
	double fBiasY = fY / (_nFocusH - 1);
	if (fBiasX < g_fDifference) return 1 - fBiasY;
	else if (fBiasY < g_fDifference) return 1 + fBiasX;
	else if (abs(fBiasX - 1) < g_fDifference) return 2 + fBiasY;
	else if (abs(fBiasY - 1) < g_fDifference) fDis = 4 - fBiasX;
	else{
		qDebug() << "error";
		return -1;
	}
}

// return the point in front using our topology order,0 or 1
int ContourStripGenerator::PointCompare(const QPointF& pt0, const QPointF& pt1){
	// first check which border they are
	int nState0 = CheckBorderPoint(pt0);
	int nState1 = CheckBorderPoint(pt1);
	double fX0 = pt0.x();
	double fY0 = pt0.y();
	double fX1 = pt1.x();
	double fY1 = pt1.y();

	if (nState0 < nState1) return 0;
	else if (nState0>nState1) return 1;
	else{
		switch (nState0)
		{
		case 0:
			if (fY0 > fY1) return 0;
			else return 1;
		case 1:
			if (fX0 < fX1) return 0;
			else return 1;
		case 2:
			if (fY0 < fY1) return 0;
			else return 1;
		case 3:
			if (fX0 > fX1) return 0;
			else return 1;
		default:
			return -1;
		}
	}
}

ContourStripGenerator::ContourStripGenerator()
{
}

ContourStripGenerator::~ContourStripGenerator()
{
}

void ContourStripGenerator::Generate(QList<UnCertaintyArea*>& areas
	, const QList<ContourLine>& contourMin
	, const QList<ContourLine>& contourMax
	, int nWidth
	, int nHeight
	, int nFocusX
	, int nFocusY
	, int nFocusW
	, int nFocusH){
	// 0.set value
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nFocusX = nFocusX;
	_nFocusY = nFocusY;
	_nFocusW = nFocusW;
	_nFocusH = nFocusH;

	arrCorners[0] = QPointF(_nFocusX, _nFocusY + _nFocusH - 1);
	arrCorners[1] = QPointF(_nFocusX, _nFocusY);
	arrCorners[2] = QPointF(_nFocusX + nFocusW - 1, _nFocusY);
	arrCorners[3] = QPointF(_nFocusX + nFocusW - 1, _nFocusY + _nFocusH - 1);

	// 1.put the contours into two lists
	for each (ContourLine contour in contourMin)
	{
		contour._nState = -1;
		if (contour._nClosed) _listClosed.push_back(contour);
		else _listUnClosed.push_back(contour);
	}
	for each (ContourLine contour in contourMax)
	{
		contour._nState = 1;
		if (contour._nClosed) _listClosed.push_back(contour);
		else _listUnClosed.push_back(contour);
	}
	// 2.handle the unclosed contours
	generateFromUnClosedContours(areas);
	// 3.handle the closed contours
	insertClosedContours(areas);
	// 4.handle the embedded area
	generateEmbeddedArea(areas);

	resetState(areas);
}

void ContourStripGenerator::generateFromUnClosedContours(QList<UnCertaintyArea*>& areas){
	// 1.generate the sorted endpoint lists
	generateEndPointsList();
	// 2.generate the contour tree
	ContourTreeNode tree = generateContourTree();
	// 3.generate the area
	int nState=generateArea(tree, areas);
	if (nState > -2 && nState < 2) resetState(areas);
}

void ContourStripGenerator::resetState(QList<UnCertaintyArea*>& areas){
	bool bFinished = false;
	while (!bFinished)
	{
		bFinished = true;
		for (int i = 0, len = areas.length(); i < len; i++)
		{
			if (areas[i]->_nState == -2 || areas[i]->_nState == 2)
			{
				// judged by parent
				if (areas[i]->_pParent != NULL && areas[i]->_pParent->_nState>-2 && areas[i]->_pParent->_nState < 2)
				{
					if (areas[i]->_pParent->_nState == 0)
					{
						areas[i]->_nState /= 2;
					}
					else{
						areas[i]->_nState = 0;
					}
					bFinished = false;
					continue;
				}


				// judged by children
				for each (UnCertaintyArea* pChild in areas[i]->_listChildren)
				{
					if (pChild->_nState > -2 && pChild->_nState < 2)
					{
						if (pChild->_nState == 0)
						{
							areas[i]->_nState /= 2;
						}
						else{
							areas[i]->_nState = 0;
						}
						bFinished = false;
						continue;
					}
				}


				// judged by embedded area
				for each (UnCertaintyArea* pChild in areas[i]->_listEmbeddedArea)
				{
					if (pChild->_nState > -2 && pChild->_nState < 2)
					{
						if (pChild->_nState == 0)
						{
							areas[i]->_nState /= 2;
						}
						else{
							areas[i]->_nState = 0;
						}
						bFinished = false;
						continue;
					}
				}
			}

		}
	}
}

void ContourStripGenerator::insertClosedContours(QList<UnCertaintyArea*>& areas){
	// handle each closed contour
	for each (ContourLine contour in _listClosed)
	{
		// generate a new area
		UnCertaintyArea* pAreaNew=new UnCertaintyArea;
		for each (QPointF pt in contour._listPt)
		{
			pAreaNew->_contour._listPt.append(pt);
		}
		pAreaNew->_contour._listPt.removeLast();	// the last one is duplicate
		pAreaNew->_contour.GenerateBounding();
		pAreaNew->_nState = contour._nState*2;		// set the state of the new area to the state of the contour, this will be used when this area generate to a standalone area
		// search every area (which can cover the whole space), put the closed contour into the one it belong to
		for (int i = 0, len = areas.length(); i < len; i++)
		{
			if (areas[i]->InsertContour(pAreaNew) == 1) break;
		}
	}
}

void ContourStripGenerator::generateEmbeddedArea(QList<UnCertaintyArea*>& areas){
	QQueue<UnCertaintyArea*> queueAreas;
	// put all the areas in a queue
	for each (UnCertaintyArea* area in areas)
	{
		queueAreas.enqueue(area);
	}
	// handle the areas in the queue, until it return false, which means it is empty
	while (handleAnAreaInQueue(queueAreas, areas));
}

bool ContourStripGenerator::handleAnAreaInQueue(QQueue<UnCertaintyArea*>& queueAreas, QList<UnCertaintyArea*>& areas){
	// finished if the queue is empty
	if (queueAreas.empty()) return false;
	// dequeue an area and handle its embedded areas 
	UnCertaintyArea* area = queueAreas.dequeue();
	bool bNewUncertain = false;
	for each (UnCertaintyArea* embedded in area->_listEmbeddedArea)
	{
		// if the outer area is uncertainty, this embedded area can be decide by the state of its contour, which has already done
		if (area->_nState == 1 || area->_nState==-1) embedded->_nState = 0;
		else if (area->_nState == 0) embedded->_nState/=2;
		else if (area->_nState != embedded->_nState) bNewUncertain = true;
		areas.append(embedded);
		queueAreas.enqueue(embedded);
	}
	if (bNewUncertain){
		area->_nState = 0;
		for each (UnCertaintyArea* embedded in area->_listEmbeddedArea)
		{
			embedded->_nState /= 2;
		}
	}

	return true;
}

void ContourStripGenerator::generateEndPointsList(){
	for (int i = 0, len = _listUnClosed.size(); i < len; i++)
	{
		// get the contour
		ContourLine contour = _listUnClosed[i];
		// get the endpoints of the contour
		QPointF endpoints[2] = {
			contour._listPt[0]
			, contour._listPt[contour._listPt.size() - 1]
		};
		for (int j = 0; j < 2; j++)
		{
			QPointF pt = endpoints[j];
			int insertPos = _listEndPoints.size();
			while (insertPos>0 && PointCompare(pt, _listEndPoints[insertPos - 1]) == 0) insertPos--;
			_listEndPoints.insert(insertPos, pt);
			_listEndPointsContourIndex.insert(insertPos, i);
			_listEndPointsIndex.insert(insertPos, j);
		}
	}
}

ContourTreeNode ContourStripGenerator::generateContourTree(){
	QList<ContourTreeNode> listTreeNode;
	for (int i = 0, len = _listEndPointsIndex.size(); i < len; i++)
	{
		ContourTreeNode treeNode;
		treeNode._nContourIndex = _listEndPointsContourIndex[i];
		treeNode._nReverse = _listEndPointsIndex[i];
		treeNode._nMatched = 0;
		listTreeNode.push_back(treeNode);
	}
	ContourTreeNode tree = buildTree(listTreeNode);
	return tree;
}

int ContourStripGenerator::generateArea(const ContourTreeNode& root, QList<UnCertaintyArea*>& areas, UnCertaintyArea* pParent){
	if (pParent==NULL){	// for root, another area should be generated
		if (root._children.empty()){	// root as leaf, means there's only two areas

		}
		else{
			// 1.generate the root area
			UnCertaintyArea* pArea=new UnCertaintyArea();
			pArea->_pParent = NULL;
			ContourLine contour = _listUnClosed[root._nContourIndex];
			generateLeafArea(pArea->_contour, contour._listPt, _listUnClosed[root._children[0]._nContourIndex]._listPt[0]);
			pArea->_nState = contour._nState * 2;
			// 2.generate the area of the root
			int nChildState = generateArea(root, areas, pArea);
			// 3.set the state of the root area based on the state of the area of the root
			if (nChildState == 0) pArea->_nState = contour._nState;
			else if (nChildState == 2 || nChildState == -2) pArea->_nState = contour._nState * 2;
			else pArea->_nState = 0;
			areas.push_back(pArea);
			return pArea->_nState;
		}
	}
	else if (root._children.empty())
	{
		// generate leaf area
		UnCertaintyArea* pArea = new UnCertaintyArea();
		pArea->_pParent = pParent;
		pParent->_listChildren.append(pArea);
		ContourLine contour = _listUnClosed[root._nContourIndex];
		generateLeafArea(pArea->_contour, contour._listPt,pParent->_contour._listPt[0]);
		if (pParent->_nState == 0) pArea->_nState = contour._nState;
		else if (pParent->_nState == 2 || pParent->_nState == -2) pArea->_nState = contour._nState * 2;
		else pArea->_nState = 0;
		areas.push_back(pArea);
		return pArea->_nState;
	}
	else{
		// generate an uncertainty area between more than one contours
		// 1.generate the area
		UnCertaintyArea* pArea = new UnCertaintyArea();
		pArea->_pParent = pParent;
		pParent->_listChildren.append(pArea);
		generateNoLeafArea(pArea->_contour, root);
		// 2.calculate state
		ContourLine contour0 = _listUnClosed[root._nContourIndex];
		bool bUncertain = false;
		for each (ContourTreeNode child in root._children)
		{
			if (_listUnClosed[child._nContourIndex]._nState != contour0._nState)
			{
				bUncertain = true;
			}
		}
		if (bUncertain) pArea->_nState = 0;								// uncertain if two contours have different state
		else if (pParent->_nState == 0) pArea->_nState = contour0._nState;	// as the state of the contour is the parent is uncertain
		else if (pParent->_nState == 2 || pParent->_nState == -2) pArea->_nState = contour0._nState * 2;
		else pArea->_nState = 0;											// can not decide yet
		// 3.generate the children area, and check whether there's child with uncertainty state
		bool bUncertainChild = false;
		for each (ContourTreeNode node in  root._children)
		{
			int nChildState = generateArea(node, areas, pArea);
			if (nChildState == 0) bUncertainChild = true;
		}
		if (pArea->_nState == -2||pArea->_nState==2)
		{
			if (bUncertainChild)
			{
				// modify state according to children
				pArea->_nState = contour0._nState;
			}
		}

		areas.append(pArea);
		return pArea->_nState;
	}
	return -1;
}

ContourTreeNode ContourStripGenerator::buildTree(QList<ContourTreeNode>& nodes){
	// 1.match nodes
	for (int ii = nodes.size() - 1; ii > 0; ii--)
	{
		if (nodes[ii]._nContourIndex == nodes[ii - 1]._nContourIndex)
		{
			// match the two nodes
			nodes[ii - 1]._nMatched = 1;
			nodes.removeAt(ii);
		}
	}
	// 2.generate children
	int nBegin = -1;
	int nState = 0;		// 1 means looking for the end point
	int i = 0;

	while (i<nodes.length())
	{
		// 2.1.check the current state
		if (nState==0)
		{
			// no begin node
			if (nodes[i]._nMatched == 0)
			{
				// find the left node
				nBegin = i;
				nState = 1;
			}
		}
		else{
			// already find begin node
			if (nodes[i]._nMatched == 0)
			{
				// find an unmatched node
				if (nodes[i]._nContourIndex==nodes[nBegin]._nContourIndex)
				{
					// the nodes matched
					// we have nodes nBegin...i
					// set the first to matched, and put the nodes between as its children, and remove the last one
					nodes[nBegin]._nMatched = 1;
					for (int j = nBegin + 1; j < i;j++)
					{
						nodes[nBegin]._children.append(nodes[j]);
					}
					for (int j = i; j > nBegin; j--)
					{
						nodes.removeAt(j);
					}
					i = 0;
					nBegin = -1;
					nState = 0;
					continue;
				}
				else{
					// can not match, restart
					nBegin = i;
				}
			}
		}
		i++;
	}
	// set the other nodes as the children of the first leaf nodes
// 	if (nodes.length()==1)
// 	{
// 		return nodes[0];
// 	}
// 	ContourTreeNode root;
// 	for (int i = 0; i < nodes.length(); i++)
// 	{
// 		if (nodes[i]._children.empty())
// 		{
// 			root = nodes[i];
// 			nodes.removeAt(i);
// 			break;
// 		}
// 	}
// 	for (int i = 0; i < nodes.length();i++)
// 	{
// 		root._children.append(nodes[i]);
// 	}
// 	return root;


	while (nodes.length()>1)
	{
		// put the children of the first node into a new list
		QList<ContourTreeNode> listNew;
		for (int i = nodes[0]._children.length() - 1; i >= 0; i--)
		{
			listNew.append(nodes[0]._children[i]);
			nodes[0]._children.removeAt(i);
		}
		// put all the nodes except the first one into the children list of the first node
		for (int i = nodes.length() - 1; i > 0; i--)
		{
			nodes[0]._children.append(nodes[i]);
			nodes.removeAt(i);
		}
		// append the first node to the new list
		listNew.append(nodes[0]);
		// use the new list to replace the original list
		nodes = listNew;
	}
	return nodes[0];
}

void ContourStripGenerator::generateNoLeafArea(ClosedContour& contour, const ContourTreeNode& node){
	// 1.generate the list of all the end points
	QList<EndPoint> listEnd;
	ContourLine contour0 = _listUnClosed[node._nContourIndex];
	EndPoint pt0;
	pt0._nIndex = 0;
	pt0._nContourIndex = node._nContourIndex;
	pt0._fDis = PointDistance(contour0._listPt[0]);
	listEnd.push_back(pt0);
	EndPoint pt1;
	pt1._nIndex = 1;
	pt1._nContourIndex = node._nContourIndex;
	pt1._fDis = PointDistance(contour0._listPt[contour0._listPt.length() - 1]);
	listEnd.push_back(pt1);
	for (int i = 0, len = node._children.length(); i < len; i++)
	{
		ContourLine contour1 = _listUnClosed[node._children[i]._nContourIndex];
		EndPoint pt0;
		pt0._nIndex = 0;
		pt0._nContourIndex = node._children[i]._nContourIndex;
		pt0._fDis = PointDistance(contour1._listPt[0]);
		listEnd.push_back(pt0);
		EndPoint pt1;
		pt1._nIndex = 1;
		pt1._nContourIndex = node._children[i]._nContourIndex;
		pt1._fDis = PointDistance(contour1._listPt[contour1._listPt.length() - 1]);
		listEnd.push_back(pt1);
	}
	// 2.Sort the end point according to their distance from the base point
	qSort(listEnd);
	int i = 0;
	int len = listEnd.length();
	// if the last and the first is a pair, put the last in the front.
	EndPoint ptLast = listEnd[len - 1];
	if (ptLast._nContourIndex == listEnd[0]._nContourIndex)
	{
		listEnd.removeLast();
		listEnd.insert(listEnd.begin(), ptLast);
	}
	// 3.connect the contour lines
	for (int i = 0; i < len; i += 2)
	{
		Q_ASSERT(listEnd[i]._nContourIndex == listEnd[i + 1]._nContourIndex);
		if (i>0)
		{
			for (int j = int(listEnd[i - 1]._fDis)+1; j <= int(listEnd[i]._fDis); j++)
			{
				contour._listPt.append(arrCorners[j]);
			}
		}
		// find an contour line, insert its points
		if (listEnd[i]._nIndex == 0)
		{
			for each (QPointF pt in _listUnClosed[listEnd[i]._nContourIndex]._listPt)
			{
				contour._listPt.append(pt);
			}
		}
		else{
			for (int j = _listUnClosed[listEnd[i]._nContourIndex]._listPt.length() - 1; j >= 0; j--)
			{
				contour._listPt.append(_listUnClosed[listEnd[i]._nContourIndex]._listPt[j]);
			}
		}
	}
	// 4.insert he last corner point
	int nLastPoint = listEnd[len - 1]._fDis;
	int nFirstPoint = listEnd[0]._fDis;
	if (nLastPoint != nFirstPoint)
	{
		int j = (nLastPoint + 1) % 4;
		while (true)
		{
			contour._listPt.append(arrCorners[j]);
			if (j == nFirstPoint) break;
			j = (j + 1) % 4;
		}
	}

	contour.GenerateBounding();
}

void ContourStripGenerator::generateLeafArea(ClosedContour& contour, const QList<QPointF>& pts, QPointF ptRef){
	QPointF pt0 = pts[0];
	QPointF pt1 = pts[pts.length() - 1];

	int nReverse = PointCompare(pt0, pt1);
	if (nReverse==1){
		QPointF ptTemp = pt0;
		pt0 = pt1;
		pt1 = ptTemp;
	}
	int n0 = PointCompare(ptRef, pt0);
	int n1 = PointCompare(ptRef, pt1);
	int nBorder0 = CheckBorderPoint(pt0);
	int nBorder1 = CheckBorderPoint(pt1);
	// add the point list into the contour
	if (nReverse==0)
	{
		for each (QPointF pt in pts)
		{
			contour._listPt.push_back(pt);
		}
	}
	else{
		for (int i = pts.length() - 1; i >= 0;i--)
		{
			contour._listPt.push_back(pts[i]);
		}
	}
	// add the border points
	if (n1==n0)	{
		for (int i = nBorder1; i > nBorder0; i--) contour._listPt.push_back(arrCorners[i]);
	}
	else{
		if (nBorder0!=nBorder1)
		{
			for (int i = nBorder1 + 1; i <= nBorder0 + 4; i++) contour._listPt.push_back(arrCorners[i % 4]);
		}
	}
	contour.GenerateBounding();
}

void ClosedContour::GenerateBounding(){
	_fX0 = _fX1 = _listPt[0].x();
	_fY0 = _fY1 = _listPt[0].y();
	for each (QPointF pt in _listPt)
	{
		if (pt.x() < _fX0) _fX0 = pt.x();
		if (pt.x() > _fX1) _fX1 = pt.x();
		if (pt.y() < _fY0) _fY0 = pt.y();
		if (pt.y() > _fY1) _fY1 = pt.y();
	}
}