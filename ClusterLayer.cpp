#include "ClusterLayer.h"

#include"MeteModel.h"
#include"LayerLayout.h"
#include "VarAnalysis.h"

ClusterLayer::ClusterLayer()
{
}


ClusterLayer::~ClusterLayer()
{
}

void ClusterLayer::draw(DisplayStates states) {
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	drawVarChart();

	// draw pca point
	drawPCAPoints();

	// draw cluster bars
	drawClusterBars();

	drawUCAreaRelation();

	glPopAttrib();
}

void ClusterLayer::init() {

}

void ClusterLayer::drawVarChart() {
	// draw chart
	double dbVarThreshold = _pModel->GetVarThreshold();

	// draw chart framework
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

	glLineWidth(1.0f);
	glColor3f(0, 0, .5);
	// scale
	for (size_t i = 1; i < _pLayout->_nVarLayers; i++)
	{
		glBegin(GL_LINES);
		glVertex2d(_pLayout->_dbVarChartLeft, _pLayout->_dbVarChartBottom + i * _pLayout->_dbVarLayer);
		glVertex2d(_pLayout->_dbVarChartRight, _pLayout->_dbVarChartBottom + i * _pLayout->_dbVarLayer);
		glEnd();
		char buf[10];
		sprintf_s(buf, "%d", i);
		_pCB->DrawText(buf, _pLayout->_dbVarChartLeft + .01, _pLayout->_dbVarChartBottom + i * _pLayout->_dbVarLayer + .01);
	}

	// draw the chart line
	glLineWidth(3.0f);
	glColor3f(0.0, .5, .5);
	std::vector<double> vecVar = _pModel->GetVariance();
	int nLen = vecVar.size();
	double dbStep = (_pLayout->_dbRight - _pLayout->_dbVarChartLeft) / (nLen - 1);
	double dbThresholdX = -1;
	glBegin(GL_LINE_STRIP);
	for (size_t i = 0; i < nLen; i += 100)
	{
		double dbYBias = vecVar[i] * _pLayout->_dbVarLayer;
		if (dbYBias>_pLayout->_dbVarChartHeight)
		{
			dbYBias = _pLayout->_dbVarChartHeight;
		}
		glVertex2d(_pLayout->_dbVarChartLeft + dbStep*i, _pLayout->_dbVarChartBottom + dbYBias);
		if (dbThresholdX<0 && vecVar[i]>dbVarThreshold) dbThresholdX = _pLayout->_dbVarChartLeft + dbStep*i;
	}
	glEnd();

	// draw the threshold line
	if (dbThresholdX > 0) {
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



	// draw nodes
	DPoint3 arrNodePosition[g_nUncertaintyAreaMax];

	int nUncertaintyAreas = std::min(_pModel->GetUncertaintyAreas(), _nUncertaintyRegions);
	double dbNodeRadius = .1;
	double dbChartRadius = _pLayout->dbChartRadius - 2 * dbNodeRadius;
	if (nUncertaintyAreas < 2) return;
	for (size_t i = 0; i < nUncertaintyAreas; i++)
	{
		double dbAngle = 2.0*M_PI*i / nUncertaintyAreas;
		double dbX = _pLayout->_dbUCARChartCenterX + dbChartRadius*cos(dbAngle);
		double dbY = _pLayout->_dbUCARChartCenterY + dbChartRadius*sin(dbAngle);
		SetGroupColor(i);
		drawCircle(dbX, dbY, dbNodeRadius, true);
		arrNodePosition[i].x = dbX;
		arrNodePosition[i].y = dbY;
	}
	// draw lines

	SetGroupColor(_nFocusedCluster);
	for (size_t i = 0; i < nUncertaintyAreas; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			glLineWidth(_pModel->GetRegionSimilarity(i,j,_nFocusedCluster));
			glBegin(GL_LINES);
			glVertex3f(arrNodePosition[i].x, arrNodePosition[i].y, 0);
			glVertex3f(arrNodePosition[j].x, arrNodePosition[j].y, 0);
			glEnd();
		}
	}

}