#include "ClusterLayer.h"

#include "MeteModel.h"
#include "LayerLayout.h"
#include "VarAnalysis.h"

ClusterLayer::ClusterLayer()
{
}

ClusterLayer::~ClusterLayer()
{
}

void ClusterLayer::draw(DisplayStates states) {
	drawVarChart();

	// draw pca point
	drawPCAPoints();

	// draw cluster bars
	drawClusterBars();

	drawUCAreaRelation();
}

void ClusterLayer::init() {

}

void ClusterLayer::drawVarChart() {
	
	double dbVarThreshold = _pModel->GetVarThreshold();			// threshold value
	std::vector<double> vecVar = _pModel->GetVariance();		// list of the variances
	int nLen = vecVar.size();									// length of the Variance vectors
	double dbVarMin = vecVar[0];								// min of variance
	double dbVarMax = vecVar[nLen - 1];							// max of variance
	int nVarVLayers = 10;// dbVarMax + 1;											// Vertical layers of the chart
	double dbVarVLayer = _pLayout->_dbVarChartHeight / nVarVLayers;			// distance between two layers



	// 1.draw chart framework
	glLineWidth(2.0f);
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	// bottom axis
	glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartTop);
	glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartTop);
	// top axis
	glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartBottom);
	glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartBottom);
	// left axis
	glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartBottom);
	glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartTop);
	// right axis
	glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartBottom);
	glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartTop);
	glEnd();

	// inner grid lines
	glLineWidth(1.0f);
	glColor3f(0, 0, .5);
	// horizontal 
	for (size_t i = 1; i < nVarVLayers; i++)
	{
		glBegin(GL_LINES);
		glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartBottom + i * dbVarVLayer);
		glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartBottom + i * dbVarVLayer);
		glEnd();
		char buf[10];
		sprintf_s(buf, "%d", i);
		_pCB->DrawText(buf, _pLayout->_dbVarChartLeft + .01, _pLayout->_dbVarChartBottom + i * dbVarVLayer + .01);
	}

	// vertical
	for (size_t i = 10; i < _pLayout->_nVarHLayers; i+=10)
	{
		glBegin(GL_LINES);
		glVertex2d(_pLayout->_dbVarChartLeft + i*_pLayout->_dbVarHLayer, _pLayout->_dbVarChartTop);
		glVertex2d(_pLayout->_dbVarChartLeft + i*_pLayout->_dbVarHLayer, _pLayout->_dbVarChartBottom);
		glEnd();
		char buf[10];
		sprintf_s(buf, "%d", i);
		_pCB->DrawText(buf, _pLayout->_dbVarChartLeft + i*_pLayout->_dbVarHLayer, _pLayout->_dbVarChartBottom - dbVarVLayer/2.0);
	}


	// 2.draw the chart line
	glLineWidth(3.0f);
	glColor3f(0.0, .5, .5);
	double dbStep = (_pLayout->_dbRight - _pLayout->_dbVarChartLeft) / (nLen - 1);
	double dbThresholdX = -1001;
	int nStep = nLen / 1000;
	int i = 0;
	glBegin(GL_LINE_STRIP);
	while(true)
	{
		double dbYBias = vecVar[i] * dbVarVLayer;
		if (dbYBias>_pLayout->_dbVarChartHeight)
		{
			dbYBias = _pLayout->_dbVarChartHeight;
		}
		glVertex2d(_pLayout->_dbVarChartLeft + dbStep*i, _pLayout->_dbVarChartBottom + dbYBias);
		if (dbThresholdX<-999 && vecVar[i]>dbVarThreshold) dbThresholdX = _pLayout->_dbVarChartLeft + dbStep*i;

		if (i==nLen-1)
		{
			break;
		}
		i += nStep;
		if (i > nLen - 1) i = nLen - 1;
	}
	glEnd();

	// 3.draw the threshold line
	if (dbThresholdX > _pLayout->_dbVarChartLeft) {
		glColor3f(1.0, 0, 0);
		glBegin(GL_LINES);
		glVertex2d(dbThresholdX, _pLayout->_dbVarChartBottom);
		glVertex2d(dbThresholdX, _pLayout->_dbVarChartTop);
		glEnd();

		char buf[10];
		sprintf_s(buf, "%.2f", dbVarThreshold);
		_pCB->DrawText(buf, dbThresholdX, _pLayout->_dbVarChartTop + .01);
	}
}

void ClusterLayer::drawPCAPoints() {

	glPointSize(3.0f);


	int nUncertaintyAreas = std::min(_pModel->GetUncertaintyAreas(), _nUncertaintyRegions);

	std::vector<DPoint3>* vecPoints = _pModel->GetPCAPoints();
	double dbMaxBias = 0;
	// calculate scale
	for (size_t clusterIndex = 0; clusterIndex < nUncertaintyAreas; clusterIndex++)
	{
		for (size_t i = 0, length = vecPoints[clusterIndex].size(); i < length; i++)
		{
			double dbX = abs(vecPoints[clusterIndex][i].x);
			double dbY = abs(vecPoints[clusterIndex][i].y);
			if (dbX > dbMaxBias) dbMaxBias = dbX;
			if (dbY > dbMaxBias) dbMaxBias = dbY;
		}
	}

	double dbScale = _pLayout->_dbPCAChartRadius*.98 / dbMaxBias;			// scale for the points

	for (size_t clusterIndex = 0; clusterIndex < nUncertaintyAreas; clusterIndex++)
	{
		double dbX = _pLayout->_dbPCAChartLeft + (clusterIndex * 2 + 1)*_pLayout->_dbPCAChartRadius;
		double dbY = _pLayout->_dbPCAChartBottom + _pLayout->_dbPCAChartRadius;

		// draw border

		// set color according to the color of the spatial cluster
		SetGroupColor(clusterIndex);

		glBegin(GL_LINE_LOOP);
		glVertex3f(dbX - (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), dbY - (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), 0);
		glVertex3f(dbX + (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), dbY - (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), 0);
		glVertex3f(dbX + (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), dbY + (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), 0);
		glVertex3f(dbX - (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), dbY + (_pLayout->_dbPCAChartRadius - _pLayout->_dbSpace), 0);
		glEnd();

		// draw the points
		glBegin(GL_POINTS);
		for (size_t i = 0, length = vecPoints[clusterIndex].size(); i < length; i++)
		{
			int nZ = vecPoints[clusterIndex][i].z;
			SetGroupColor(nZ);
			glVertex3f(dbX + vecPoints[clusterIndex][i].x*dbScale, dbY + vecPoints[clusterIndex][i].y*dbScale, 0);
		}
		glEnd();
	}
}

void ClusterLayer::drawClusterBars() {
	const ClusterResult* pCR = _pModel->GetClusterResults();
	if (pCR->_nK == 0) return;

	int nUncertaintyAreas = std::min(_pModel->GetUncertaintyAreas(), _nUncertaintyRegions);

	double arrBaseY[50];			// the y position of the element in the first cluster

	double dbBarLen = _pLayout->_dbPCAChartRadius / 2.0;		// length of the bar
	for (size_t nIndex = 0; nIndex < nUncertaintyAreas; nIndex++)
	{
		double dbY = _pLayout->_dbBarChartBottom + _pLayout->_dbSpace;
		double dbX0 = _pLayout->_dbBarChartLeft + nIndex * 2 * _pLayout->_dbPCAChartRadius - dbBarLen;
		double dbX1 = _pLayout->_dbBarChartLeft + nIndex * 2 * _pLayout->_dbPCAChartRadius + dbBarLen;
		double dbX2 = _pLayout->_dbBarChartLeft + nIndex * 2 * _pLayout->_dbPCAChartRadius + 3 * dbBarLen;
		// draw the link
		if (nIndex>0)
		{
			for (size_t i = 0; i < pCR[nIndex]._nK; i++)
			{
				for (size_t j = 0; j < pCR[nIndex]._vecItems[i].size(); j++) {
					SetGroupColor(i);

					glBegin(GL_LINES);
					glVertex3f(dbX0, arrBaseY[pCR[nIndex]._vecItems[i][j]], 0);
					glVertex3f(dbX1, dbY, 0);
					glEnd();
					dbY += _pLayout->_dbSpace;
				}
				dbY += _pLayout->_dbSpaceII;
			}
		}
		dbY = _pLayout->_dbBarChartBottom + _pLayout->_dbSpace;
		// draw this region
		for (size_t i = 0; i < pCR[nIndex]._nK; i++)
		{
			SetGroupColor(i);

			for (size_t j = 0; j < pCR[nIndex]._vecItems[i].size(); j++) {
				glBegin(GL_LINES);
				glVertex3f(dbX1, dbY, 0);
				glVertex3f(dbX2, dbY, 0);
				glEnd();
				arrBaseY[pCR[nIndex]._vecItems[i][j]] = dbY;
				dbY += _pLayout->_dbSpace;
			}
			dbY += _pLayout->_dbSpaceII;
		}
	}
}

void ClusterLayer::drawUCAreaRelation() {
	// draw border
	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);
	glVertex3f(_pLayout->_dbUCARChartLeft, _pLayout->_dbUCARChartTop, 0);
	glVertex3f(_pLayout->_dbUCARChartRight, _pLayout->_dbUCARChartTop, 0);
	glVertex3f(_pLayout->_dbUCARChartRight, _pLayout->_dbUCARChartBottom, 0);
	glVertex3f(_pLayout->_dbUCARChartLeft, _pLayout->_dbUCARChartBottom, 0);
	glEnd();


	int nUncertaintyAreas = std::min(_pModel->GetUncertaintyAreas(), _nUncertaintyRegions);
	if (nUncertaintyAreas < 2) return;

	// calculate node position
	DPoint3 arrNodePosition[g_nUncertaintyAreaMax];
	double dbNodeRadius = .1;
	double dbChartRadius = _pLayout->dbChartRadius - 2 * dbNodeRadius;
	for (size_t i = 0; i < nUncertaintyAreas; i++)
	{
		double dbAngle = 2.0*M_PI*i / nUncertaintyAreas;
		double dbX = _pLayout->_dbUCARChartCenterX + dbChartRadius*cos(dbAngle);
		double dbY = _pLayout->_dbUCARChartCenterY + dbChartRadius*sin(dbAngle);

		arrNodePosition[i].x = dbX;
		arrNodePosition[i].y = dbY;
	}

	// draw lines
	SetGroupColor(_nFocusedCluster);
	for (size_t i = 0; i < nUncertaintyAreas; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			glLineWidth(_pModel->GetRegionSimilarity(i, j, _nFocusedCluster));
			glBegin(GL_LINES);
			glVertex3f(arrNodePosition[i].x, arrNodePosition[i].y, 0);
			glVertex3f(arrNodePosition[j].x, arrNodePosition[j].y, 0);
			glEnd();
		}
	}

	// draw nodes
	for (size_t i = 0; i < nUncertaintyAreas; i++)
	{
		SetGroupColor(i);
		drawCircle(arrNodePosition[i].x, arrNodePosition[i].y, dbNodeRadius, true);
	}

}