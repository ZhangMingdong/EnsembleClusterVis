#include "myglwidget.h"
#include <gl/GLU.h>

#include <QMouseEvent>
#include <QImage>

#include <iostream>
#include <QDebug>
#include "ContourStripGenerator.h"

#include "MeteLayer.h"
#include "CostLineLayer.h"
#include "EnsembleLayer.h"
#include "MeteModel.h"



using namespace std;
double g_arrColors[20][3] = {
// 	{ 0, 0, 0 },			// R
// 	{ 0, 1, 0 },			// G
// 	{ 1, 1, 1 },			// B

	{ 1, 0, 0 },			// R
	{ 0, 1, 0 },			// G
	{ 0, 0, 1 },			// B

	{ 1, 1, 0 },			// yellow
	{ 0, 1, 1 },			// purple
	{ 1, 0, 1 },			// 

	{ .6, .3, 1 },			// R
	{ 1, .6, .3 },			// R
	{ .3, 1, .6 },			// R

	{ .3, .6, 1 },			// R
	{ 1, .3, .6 },			// R
	{ .6, 1, .3 },			// R

	{ .5, 1, 1 },			// R
	{ 1, .5, 1 },			// R
	{ 1, 1, .5 },			// R

	{ .5, 1, 0 },			// R
	{ 0, .5, 1 },			// R
	{ 1, 0, .5 },			// R

	{ .2, .7, .2 },			// R
	{ .2, .2, .7 },			// R

};

class GridStatistic
{
public:
	double _fMin;
	double _fMean;
	double _fMax;
	bool operator<(const GridStatistic& t) const{
		return _fMean < t._fMean;
	}
	bool operator>(const GridStatistic& t) const{
		return _fMean > t._fMean;
	}
	GridStatistic(double fMin, double fMax) :_fMin(fMin), _fMax(fMax),_fMean((fMax+fMin)/2){};
};




void vertorMult(const double* m1, const double* m2, double* result)
{
	result[0] = m1[1] * m2[2] - m2[1] * m1[2];//y1*z2 - y2*z1;
	result[1] = m1[2] * m2[0] - m2[2] * m1[0];//z1*x2 - z2*x1;
	result[2] = m1[0] * m2[1] - m2[0] * m1[1];//x1*y2 - x2*y1;
}
// result = pt0pt1*pt0pt2
void vectorMult(const double* pt0, const double* pt1, const double* pt2, double*result)
{
	double m1[3], m2[3];

	m1[0] = pt1[0] - pt0[0];
	m1[1] = pt1[1] - pt0[1];
	m1[2] = pt1[2] - pt0[2];
	m2[0] = pt2[0] - pt0[0];
	m2[1] = pt2[1] - pt0[1];
	m2[2] = pt2[2] - pt0[2];
	vertorMult(m1, m2, result);
}

MyGLWidget::MyGLWidget(QWidget *parent)
: QGLWidget(parent)
, m_dbRadius(1000)
, m_nMouseState(0)
, m_nRank(3)
, c_dbPI(3.14159265)
, c_dbSideLen(1.0)
, _pDataT(NULL)
, _pDataB(NULL)
, _nSelectedX(10)
, _nSelectedY(10)
, _bShowGridLines(false)
, _bShowIntersection(false)
, _bShowUnionB(false)
, _bShowClusterBS(false)
, _bShowClusterBV(false)
, _bShowLineChart(false)
, _bShowContourLineTruth(false)
, _nSelectedLeft(-1)
, _nSelectedRight(-1)
, _nSelectedTop(-1)
, _nSelectedBottom(-1)
{
	_fLeft = -1.80;
	_fRight = 1.80;
	_fBottom = -0.90;
	_fTop = 0.90;

	_fScaleW = 0.01;
	_fScaleH = 0.01;

	_fChartLeft = 0.54;
	_fChartW = .5;
	_fChartRight = _fChartLeft + _fChartW;
	
	_fChartLRight = -0.54;
	_fChartLW = .5;
	_fChartLLeft = _fChartLRight-_fChartLW;

	startTimer(100);
}

MyGLWidget::~MyGLWidget()
{
// 	if (_dataTexture)
// 	{
// 		delete[]_dataTexture;
// 	}
	// release the uncertainty area
// 	for (int i = 0; i < g_temperatureLen; i++)
// 		for each (UnCertaintyArea* pArea in _listUnionAreaE[i])
// 			delete pArea;
	for (int i = 0; i < g_temperatureLen; i++)
		for each (UnCertaintyArea* pArea in _listIntersectionAreaB[i]) 
			delete pArea;
	for (int i = 0; i < g_gradient_l; i++)
		for each (UnCertaintyArea* pArea in _listUnionAreaB[i])
			delete pArea;

	// release the layers
	for each (MeteLayer* layer in _vecLayers)
	{
		delete layer;
	}
}

void MyGLWidget::initializeGL()
{
// 	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearColor(1.0, 1.0, 1.0, 1.0);
//	glEnable(GL_DEPTH_TEST);
// 	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// used for border
	glLineWidth(2.0);

	// enable blending
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glLineWidth(1.5);	



//	MeteLayer* pLayer = new EnsembleLayer();
//	pLayer->SetModel(_pModelT);
//	pLayer->InitLayer(_fLeft, _fRight, _fTop, _fBottom, _fScaleW, _fScaleH);
//	_vecLayers.push_back(pLayer);	

	MeteLayer* pLayer = new EnsembleLayer();
	pLayer->SetModel(_pModelE);
	pLayer->InitLayer(_fLeft, _fRight, _fTop, _fBottom, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);

	// create cost line layer
	pLayer = new CostLineLayer();
	pLayer->InitLayer(_fLeft, _fRight, _fTop, _fBottom, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);


	font.InitFont(wglGetCurrentDC(), L"Arial", 15);
	MeteLayer::_pCB = this;
}

void MyGLWidget::paintGL(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw selected point
	glLineWidth(2.0f);
	glColor4f(1.0, .3, 0.0, 1.0);
	double fX0 = _fLeft + _nSelectedLeft*_fScaleW;
	double fY0 = _fBottom + _nSelectedBottom*_fScaleH;
	double fX1 = _fLeft + _nSelectedRight*_fScaleW;
	double fY1 = _fBottom + _nSelectedTop*_fScaleH;
	glBegin(GL_LINE_LOOP);
	glVertex2f(fX0,fY0);
	glVertex2f(fX0,fY1);
	glVertex2f(fX1,fY1);
	glVertex2f(fX1,fY0);
	glEnd();


	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->DrawLayer(_displayStates);
	}

// 	font.PrintText("123", 0, 0);
	// Render quad


	if (_bShowUnionB)
	{
// 		glColor4f(1.0, 1.0, 0.0, 0.05);
		for (int i = 0; i < g_gradient_l;i++)
		{
			glColor4f(g_arrColors[0][0], g_arrColors[0][1], g_arrColors[0][2], 0.05);
			glCallList(_gllist + i * 3 + 0);
			glColor4f(g_arrColors[1][0], g_arrColors[1][1], g_arrColors[1][2], 0.05);
			glCallList(_gllist + i * 3 + 1);
			glColor4f(g_arrColors[2][0], g_arrColors[2][1], g_arrColors[2][2], 0.05);
			glCallList(_gllist + i * 3 + 2);
		}
	}

	if (_bShowIntersection)
	{
		double fTransparency = .5;
		for (int i = 0; i < g_temperatureLen; i++)
		{
			glColor4f(g_arrColors[0][0], g_arrColors[0][1], g_arrColors[0][2], fTransparency);
			glCallList(_gllist + g_gradient_l * 3 + g_temperatureLen * 3 + i * 3 + 0);
			glColor4f(g_arrColors[1][0], g_arrColors[1][1], g_arrColors[1][2], fTransparency);
			glCallList(_gllist + g_gradient_l * 3 + g_temperatureLen * 3 + i * 3 + 1);
			glColor4f(g_arrColors[2][0], g_arrColors[2][1], g_arrColors[2][2], fTransparency);
			glCallList(_gllist + g_gradient_l * 3 + g_temperatureLen * 3 + i * 3 + 2);
		}

	}

	// contour line of truth data
	if (_bShowContourLineTruth)
	{
		glLineWidth(2);
		glColor4f(0.0, 0.0, 0.0, 1.0); 
		for (int i = 0; i < g_temperatureLen; i++)
		{
			drawContourLine(_listContourTruth[i]);
		}
		glLineWidth(1);
	}


// 	double fStepW = 1.0 / (g_focus_w - 1);
// 	double fStepH = (_fTop-_fBottom) / (g_focus_h - 1);

	// grid lines
	if (_bShowGridLines){
		
		glColor4f(.5, .5, .5, .5);
		glBegin(GL_LINES);
		int nStep = 10;
		for (int i = 0; i < g_globalW; i += nStep)
		{
			glVertex2f(_fLeft + i*_fScaleW, _fBottom);
			glVertex2f(_fLeft + i*_fScaleW, _fTop);
		}
		for (int j = 0; j < g_globalH; j += nStep)
		{
			glVertex2f(_fLeft, _fBottom + j*_fScaleH);
			glVertex2f(_fRight, _fBottom + j*_fScaleH);
		}
		glEnd();
		// show detail
		if (false) {	
			glColor4f(.5, .5, .5, .3);
			glBegin(GL_LINES);
			int nStep = 1;
			for (int i = 0; i < g_globalW; i += nStep)
			{
				glVertex2f(_fLeft + i*_fScaleW, _fBottom);
				glVertex2f(_fLeft + i*_fScaleW, _fTop);
			}
			for (int j = 0; j < g_globalH; j += nStep)
			{
				glVertex2f(_fLeft, _fBottom + j*_fScaleH);
				glVertex2f(_fRight, _fBottom + j*_fScaleH);
			}
			glEnd();
		}
	}



	double fMargin = .02;

	if (_bShowLineChart)
	{

		// draw selected point
		glLineWidth(2.0f);
		glColor4f(1.0, 1.0, 1.0, 1.0);
		double fX = _fLeft + _nSelectedX*_fScaleW;
		double fY = _fBottom + _nSelectedY*_fScaleH;
		glBegin(GL_LINE_LOOP);
		glVertex2f(fX - 0.5*_fScaleW, fY - 0.5*_fScaleH);
		glVertex2f(fX + 0.5*_fScaleW, fY - 0.5*_fScaleH);
		glVertex2f(fX + 0.5*_fScaleW, fY + 0.5*_fScaleH);
		glVertex2f(fX - 0.5*_fScaleW, fY + 0.5*_fScaleH);
		glEnd();

		// draw the statistic stack graph
		QList<GridStatistic> list;
		for (int ll = 0; ll < g_lenEnsembles; ll++)
		{
			const double *data = _pDataB + ll * 4 * g_focus_l + (_nSelectedY*g_focus_w + _nSelectedX) * 4;
			double mb = data[0];
			double mv = data[1];
			double m = data[2];
			double v = data[3];
			// 			qDebug() << mb << "\t" << mb << "\t" << m << "\t" << v;
			double fMin = m - v;
			double fMax = m + v;
			double fMean = m;
			fMin = fMin*mv + mb;
			fMax = fMax*mv + mb;
			fMean = fMean*mv + mb;
			list.append(GridStatistic(fMin, fMax));
			// 			qDebug() << fMin << "\t" << fMax << "\t" << fMean;

		}
		qSort(list);
		int nMinT = -70;
		int nMaxT = 30;
		double fMinT = nMinT;
		double fMaxT = nMaxT;
		for (int ll = 0; ll < g_lenEnsembles; ll++)
		{
			double fMin = list[ll]._fMin;
			double fMax = list[ll]._fMax;
			double fX1 = _fChartLeft + (fMin - fMinT) / (fMaxT - fMinT)*_fChartW;
			double fX2 = _fChartLeft + (fMax - fMinT) / (fMaxT - fMinT)*_fChartW;

			glColor4f(0, 0.8, 0, 0.5);
			glBegin(GL_QUADS);
			glVertex2f(fX1, _fBottom + (ll + 1)*(_fTop - _fBottom) / 50.0);
			glVertex2f(fX2, _fBottom + (ll + 1)*(_fTop - _fBottom) / 50.0);
			glVertex2f(fX2, _fBottom + ll*(_fTop - _fBottom) / 50.0);
			glVertex2f(fX1, _fBottom + ll*(_fTop - _fBottom) / 50.0);

			// 			glVertex2f(fMin*0.01,0.51+(ll-1)*0.01);
			// 			glVertex2f(fMax*0.01,0.51+ (ll - 1)*0.01);
			// 			glVertex2f(fMax*0.01,0.51+ ll*0.01);
			// 			glVertex2f(fMin*0.01,0.51+ ll*0.01);
			glEnd();
			glColor4f(0, 0.8, 0, 1.0);
			glBegin(GL_LINES);
			glVertex2f((fX1 + fX2) / 2, _fBottom + (ll + 1)*(_fTop - _fBottom) / 50.0);
			glVertex2f((fX1 + fX2) / 2, _fBottom + ll*(_fTop - _fBottom) / 50.0);
			glEnd();
		}
		// draw the frame
		glBegin(GL_LINES);
		glLineWidth(2.0);
		glColor3f(.5, .5, .5);
		glVertex2f(_fChartLeft, _fBottom);
		glVertex2f(_fChartRight, _fBottom);
		glVertex2f(_fChartLeft, _fTop);
		glVertex2f(_fChartLeft, _fBottom);
		glEnd();
		// x coordinate
		for (int i = 0, j = nMinT; j <= nMaxT; i += 10, j += 10)
		{
			char buf[5];
			sprintf_s(buf, "%d", j);
			font.PrintText(buf, _fChartLeft + _fChartW*i / (nMaxT - nMinT), _fBottom - fMargin);
		}
		// y coordinate
		for (int i = 10; i <= g_lenEnsembles; i += 10)
		{
			char buf[5];
			sprintf_s(buf, "%d", i);
			font.PrintText(buf, _fChartLeft - fMargin, _fBottom + (_fTop - _fBottom)*i / g_lenEnsembles);

		}

		// draw intersection line chart
		EnsembleIntersections intersection = _arrIntersections[_nSelectedY*g_focus_w + _nSelectedX];
		int nLen = intersection._listIntersection.length();
		// 		double fMin = intersection._listIntersection[0]._fMin;
		// 		double fMax = intersection._listIntersection[nLen - 1]._fMax;
		// 		double fRange = fMax - fMin;
		// 
		// 		glColor3f(1.0, 0, 0);
		// 		glBegin(GL_LINE_STRIP);
		// 		// 		qDebug() << "==================";
		// 		double fLastOverlap = 0;
		// 		glVertex2f(_fChartLeft, _fBottom + (_fTop - _fBottom)*intersection._listIntersection[0]._nOverlap / g_lenEnsembles);
		// 		for (int i = 0; i < nLen; i++)
		// 		{
		// 			double fValue = intersection._listIntersection[i]._fMax;
		// 			double fOverlap = intersection._listIntersection[i]._nOverlap;
		// 			double fX = _fChartLeft + (fValue - fMin) / fRange;
		// 			if (i>0) glVertex2f(fX, _fBottom + (_fTop - _fBottom)*fLastOverlap / g_lenEnsembles);
		// 
		// 			glVertex2f(fX, _fBottom + (_fTop - _fBottom)*fOverlap / g_lenEnsembles);
		// 			fLastOverlap = fOverlap;
		// 			qDebug() << fOverlap;
		// 		}
		// 		glEnd();

		double fMin = nMinT;
		double fMax = nMaxT;
		double fRange = fMax - fMin;

		glColor3f(1.0, 0, 0);
		glBegin(GL_LINE_STRIP);
		// 		qDebug() << "==================";
		// left end
		glVertex2f(_fChartLeft, _fBottom);
		glVertex2f(_fChartLeft + (intersection._listIntersection[0]._fMin - fMin) / fRange*_fChartW, _fBottom);
		// 		glVertex2f(_fChartLeft + (intersection._listIntersection[0]._fMax - fMin) / fRange, _fBottom);
		// 
		// 		qDebug() << intersection._listIntersection[0]._fMin;
		// 		qDebug() << intersection._listIntersection[0]._fMax;

		for (int i = 0; i < nLen; i++)
		{
			double fValue1 = intersection._listIntersection[i]._fMin;
			double fValue2 = intersection._listIntersection[i]._fMax;
			double fOverlap = intersection._listIntersection[i]._nOverlap;
			double fX1 = _fChartLeft + (fValue1 - fMin) / fRange*_fChartW;
			double fX2 = _fChartLeft + (fValue2 - fMin) / fRange*_fChartW;
			glVertex2f(fX1, _fBottom + (_fTop - _fBottom)*fOverlap / g_lenEnsembles);
			glVertex2f(fX2, _fBottom + (_fTop - _fBottom)*fOverlap / g_lenEnsembles);
			// 			qDebug() << fOverlap;
		}

		glVertex2f(_fChartLeft + (intersection._listIntersection[intersection._listIntersection.length() - 1]._fMax - fMin) / fRange*_fChartW, _fBottom);
		// right end
		glVertex2f(_fChartRight, _fBottom);
		glEnd();


		// draw the threshold line
		glColor3f(.5, 0, 0);
		glBegin(GL_LINE_STRIP);
		glVertex2f(_fChartLeft, _fBottom + (_fTop - _fBottom)*g_nThreshold / g_lenEnsembles);
		glVertex2f(_fChartRight, _fBottom + (_fTop - _fBottom)*g_nThreshold / g_lenEnsembles);
		glEnd();
	}
	// draw clustering
	if (_bShowClusterBS)
	{
		double colors[20][3] = {
			{ 1, 0, 0 },			// R
			{ 0, 1, 0 },			// G
			{ 0, 0, 1 },			// B
	
			{ 1, 1, 0 },			// yellow
			{ 0, 1, 1 },			// purple
			{ 1, 0, 1 },			// 

			{ .6, .3, 1 },			// R
			{ 1, .6, .3 },			// R
			{ .3, 1, .6 },			// R

			{ .3, .6, 1 },			// R
			{ 1, .3, .6 },			// R
			{ .6, 1, .3 },			// R
		
			{ .5, 1, 1 },			// R
			{ 1, .5, 1 },			// R
			{ 1, 1, .5 },			// R
		
			{ .5, 1, 0 },			// R
			{ 0, .5, 1 },			// R
			{ 1, 0, .5 },			// R
			
			{ .2, .7, .2 },			// R
			{ .2, .2, .7 },			// R

		};

		glBegin(GL_QUADS);
		for (int i = 0; i < g_focus_h; i++)
		{
			for (int j = 0; j < g_focus_w; j++)
			{

				double fX = _fLeft + j*_fScaleW;
				double fY = _fBottom + i*_fScaleH;
				double* color = colors[_arrLabels[i*g_focus_w + j]];
				glColor4f(color[0] * .5, color[1] * .5, color[2] * .5, 1);
				glVertex2f(fX - 0.5*_fScaleW, fY - 0.5*_fScaleH);
				glVertex2f(fX + 0.5*_fScaleW, fY - 0.5*_fScaleH);
				glVertex2f(fX + 0.5*_fScaleW, fY + 0.5*_fScaleH);
				glVertex2f(fX - 0.5*_fScaleW, fY + 0.5*_fScaleH);
// 				qDebug() << _arrLabels[i*g_focus_h + j];
			}
		}
		glEnd();


		// draw the cluster chart



		// draw the frame
		glBegin(GL_LINES);
		glLineWidth(2.0);
		glColor3f(.5, .5, .5);
		glVertex2f(_fChartLLeft, _fBottom);
		glVertex2f(_fChartLRight, _fBottom);
		glVertex2f(_fChartLLeft, _fTop);
		glVertex2f(_fChartLLeft, _fBottom);
		glLineWidth(.5);
// 		glColor3f(.2, .2, .2);
		glVertex2f(_fChartLLeft, _fBottom + (_fTop - _fBottom)/2);
		glVertex2f(_fChartLRight, _fBottom + (_fTop - _fBottom) / 2);
		glVertex2f(_fChartLLeft + _fChartLW/2, _fTop);
		glVertex2f(_fChartLLeft + _fChartLW/2, _fBottom);
		glEnd();



// 		double fMargin = .02;
// 		// x coordinate
// 		for (int i = 0, j = nMinT; j <= nMaxT; i += 10, j += 10)
// 		{
// 			char buf[5];
// 			sprintf_s(buf, "%d", j);
// 			font.PrintText(buf, _fChartLeft + _fChartW*i / (nMaxT - nMinT), _fBottom - fMargin);
// 		}
// 		// y coordinate
// 		for (int i = 10; i <= g_lenEnsembles; i += 10)
// 		{
// 			char buf[5];
// 			sprintf_s(buf, "%d", i);
// 			font.PrintText(buf, _fChartLeft - fMargin, _fBottom + (_fTop - _fBottom)*i / g_lenEnsembles);
// 
// 		}

		font.PrintText("0", _fChartLLeft - fMargin, _fBottom - fMargin);
		font.PrintText("1", _fChartLLeft + _fChartLW, _fBottom - fMargin);
		font.PrintText("1", _fChartLLeft - fMargin, _fBottom + (_fTop - _fBottom));


		font.PrintText(".5", _fChartLLeft + _fChartLW/2, _fBottom - fMargin);
		font.PrintText(".5", _fChartLLeft - fMargin, _fBottom + (_fTop - _fBottom)/2);



// 		double fVMMin = 1000;
// 		double fMVMin = 1000;
// 		double fVMMax = -1000;
// 		double fMVMax = -1000;
// 		for (int i = 0; i < g_focus_l; i++){
// 			if (_gridDataMeanVarB[i] > fMVMax) fMVMax = _gridDataMeanVarB[i];
// 			if (_gridDataMeanVarB[i] < fMVMin) fMVMin = _gridDataMeanVarB[i];
// 			if (_gridDataVarMeanB[i] > fVMMax) fVMMax = _gridDataVarMeanB[i];
// 			if (_gridDataVarMeanB[i] < fVMMin) fVMMin = _gridDataVarMeanB[i];
// 		}
// 		double fVMRange = fVMMax - fVMMin;
		// 		double fMVRange = fMVMax - fMVMin;
		glPointSize(2.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < g_focus_l;i++)
		{
			double* color = colors[_arrLabels[i]];
// 			glColor3f(color[0], color[1], color[2]);
// 			glVertex2f(_fChartLLeft + (_gridDataMeanVarB[i] - fMVMin)*_fChartLW / fMVRange,
// 				_fBottom + (_gridDataVarMeanB[i] - fVMMin)*(_fTop-_fBottom) / fVMRange);

			glColor3f(color[0] * .5, color[1] * .5, color[2] * .5);
			glVertex2f(_fChartLLeft + _gridDataMeanVarB[i]*_fChartLW,
				_fBottom + _gridDataVarMeanB[i]*(_fTop - _fBottom));

		}
		// draw selected point
// 		glPointSize(4.0);
// 		int index = _nSelectedY*g_focus_w + _nSelectedX;
// 		double* color = colors[_arrLabels[index]];
// 
// 		glColor3f(color[0], color[1], color[2]);
// 		glVertex2f(_fChartLLeft + _gridDataMeanVarB[index] * _fChartLW,
// 			_fBottom + _gridDataVarMeanB[index] * (_fTop - _fBottom));
		glEnd();

		// draw the selected point
		glColor3f(1, 1, 1);
		glLineWidth(.5);
		int index = _nSelectedY*g_focus_w + _nSelectedX;
		double fX = _fChartLLeft + _gridDataMeanVarB[index] * _fChartLW;
		double fY = _fBottom + _gridDataVarMeanB[index] * (_fTop - _fBottom);
		glBegin(GL_LINES);
		glVertex2f(_fChartLLeft, fY);
		glVertex2f(fX, fY);
		glVertex2f(fX, fY);
		glVertex2f(fX, _fBottom);
		glEnd();

	}



	// draw the surface
	if (false)
	{
		double fMin = 100000;
		double fMax = -100000;

		for (int i = 0; i < g_focus_l; i++)
		{
			if (_pDataT[i]>fMax) fMax = _pDataT[i];
			if (_pDataT[i] < fMin) fMin = _pDataT[i];
		}
		double fRange = fMax - fMin;
		qDebug() << fMin;
		qDebug() << fMax;
		qDebug() << fRange;

		glColor4f(0, 0, .5, .5);
		glBegin(GL_QUADS);
		for (int i = 0; i < g_focus_h - 1; i++)
		{
			for (int j = 0; j < g_focus_w - 1; j++)
			{
				double fX1 = _fLeft + j*_fScaleW;
				double fY1 = _fBottom + i*_fScaleH;
				double fX2 = _fLeft + (j + 1)*_fScaleW;
				double fY2 = _fBottom + (i + 1)*_fScaleH;
				int nZ11 = i*g_focus_w + j;
				int nZ12 = (i + 1)*g_focus_w + j;
				int nZ21 = i*g_focus_w + j + 1;
				int nZ22 = (i + 1)*g_focus_w + j + 1;
				double z11 = (_pDataT[nZ11] - fMin)*.1 / fRange;
				double z12 = (_pDataT[nZ12] - fMin)*.1 / fRange;
				double z21 = (_pDataT[nZ21] - fMin)*.1 / fRange;
				double z22 = (_pDataT[nZ22] - fMin)*.1 / fRange;
				glVertex3f(fX1, fY1, z11);
				glVertex3f(fX1, fY2, z12);
				glVertex3f(fX2, fY2, z22);
				glVertex3f(fX2, fY1, z21);
				// 			qDebug() << z11 << z12 << z21 << z22;
			}
		}
		glEnd();
	}



}

void MyGLWidget::drawContourLine(const QList<ContourLine>& contours){

	for each (ContourLine contour in contours)
	{
		glBegin(GL_LINE_STRIP);
		for each (QPointF pt in contour._listPt)
		{
			double x = _fLeft + pt.x() * _fScaleW;
			double y = _fTop - pt.y() * _fScaleH;
			glVertex2f(x, y);
		}
		glEnd();
	}
}

void MyGLWidget::resizeGL(int width, int height){
	_nWidth = width;
	_nHeight = height;
	// 1.viewport
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);

	// 2.projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(30.0, (GLfloat)width / (GLfloat)height, .1, 100.0);

	reset();

}

void MyGLWidget::reset()
{
	// 3.view
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, -0.8*m_nRank);
	/*
	glRotated(45, 1, 0, 0);
	glRotated(45, 0, 1, 0);
	*/
	// 4.record the model view matrix
	glGetDoublev(GL_MODELVIEW_MATRIX, m_modelViewMatrixGlobal);
	// 5.init mapper and axis


}

void MyGLWidget::timerEvent(QTimerEvent* event)
{

}

void MyGLWidget::mousePressEvent(QMouseEvent * event)
{
	if (event->button() == Qt::MidButton)
	{
	}
	else if (event->button()==Qt::LeftButton)
	{
		updateTrackBallPos(event->pos(), m_dbTrackBallPos);
		m_nMouseState = 1;
		_ptLast = event->pos();
	}
	else if (event->button() == Qt::RightButton)
	{
		m_nMouseState = 2;
		select(_nSelectedLeft, _nSelectedBottom,event->pos());
		_nSelectedRight = _nSelectedLeft;
		_nSelectedTop = _nSelectedBottom;
	}
	updateGL();
}

void MyGLWidget::mouseReleaseEvent(QMouseEvent * event)
{
	if (m_nMouseState == 2)
	{
		select(_nSelectedRight, _nSelectedTop, event->pos());
		if (_nSelectedRight < _nSelectedLeft)
		{
			int nTemp = _nSelectedRight;
			_nSelectedRight = _nSelectedLeft;
			_nSelectedLeft = nTemp;
		}
		if (_nSelectedTop < _nSelectedBottom)
		{
			int nTemp = _nSelectedTop;
			_nSelectedTop = _nSelectedBottom;
			_nSelectedBottom = nTemp;
		}

		if (_nSelectedRight > _nSelectedLeft&&_nSelectedTop>_nSelectedBottom)
		{
			for each (MeteLayer* pLayer in _vecLayers)
			{
				pLayer->Brush(_nSelectedLeft, _nSelectedRight, _nSelectedTop, _nSelectedBottom);
			}
		}
	}
	m_nMouseState = 0;
	updateGL();

}

void MyGLWidget::mouseMoveEvent(QMouseEvent *event)
{
	if (m_nMouseState==2)
	{
		select(_nSelectedRight, _nSelectedTop, event->pos());
	}
	else if (m_nMouseState==1){	// left button down
		if (false){

			// calculate current pos
			double currentPos[3];
			updateTrackBallPos(event->pos(), currentPos);
			// calculate axis
			double axis[3];
			double center[3] = { 0, 0, 0 };
			vectorMult(center, currentPos, m_dbTrackBallPos, axis);

			// calculate angle
			double d[3];
			for (int i = 0; i < 3; i++)
			{
				d[i] = m_dbTrackBallPos[i] - currentPos[i];
			}
			double b = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			double angle = 2 * asin(b / 2 / m_dbRadius) * 180 / c_dbPI;

			// Rotate
			rotate(axis, angle);
			for (int i = 0; i < 3; i++)
			{
				m_dbTrackBallPos[i] = currentPos[i];
			}
		}
		else{
			move(_ptLast, event->pos());
			_ptLast = event->pos();
		}
		updateGL();
	}
}


void MyGLWidget::select(int& nX, int&nY, const QPoint& pt) {


	// last pickIndex
	PickIndex pick;
	// 1.用glSelectBuffer()函数指定用于返回点击记录的数组
	GLuint selectBuf[BUFSIZE];
	glSelectBuffer(BUFSIZE, selectBuf);

	// 2.用glRenderMode()指定GL_SELECT，进入选择模式
	(void)glRenderMode(GL_SELECT);

	// 3.使用glInitNames()和glPushName()对名字堆栈进行初始化
	glInitNames();
	glPushName(0);

	// 4.定义用于选择的视景体
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	// create 5*5 pixel picking region near cursor location
	gluPickMatrix((GLdouble)pt.x(), (GLdouble)(viewport[3] - pt.y()), 1.0, 1.0, viewport);
	gluPerspective(30.0, (GLfloat)width() / (GLfloat)height(), .1, 100.0);
	glMatrixMode(GL_MODELVIEW);
	// 5.交替调用绘制图元的函数和操纵名字栈的函数，为每个相关的图元分配一个适当的名称
	drawPlaceHolder();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glFlush();

	// 6. 退出选择模式，并处理返回的选择数据
	GLint hits = glRenderMode(GL_RENDER);
	if (processHits(hits, selectBuf, pick))
	{
		if (pick._nX < g_focus_w) {
			nX = pick._nX;
			nY = pick._nY;
			qDebug() << pick._nX << "\t" << pick._nY << "\t" << _pDataT[nY*g_focus_w + nX];
		}
		else {
			int nIndex = pick._nX - g_focus_w;
			for each (MeteLayer* pLayer in _vecLayers)
			{
				pLayer->OnSelectVar(nIndex);
			}
		}
		
		updateGL();
	}
}
void MyGLWidget::mouseDoubleClickEvent(QMouseEvent *event){
	select(_nSelectedX, _nSelectedRight, event->pos());
}

void MyGLWidget::wheelEvent(QWheelEvent * event)
{
	double fScale;
	if (event->delta()>0)
	{
		fScale = 1.1;
	}
	else
	{
		fScale =0.9;
	}
			
	double t[3] = { m_modelViewMatrixGlobal[12], m_modelViewMatrixGlobal[13], m_modelViewMatrixGlobal[14] };//得到x、y、z轴的平移量

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(t[0], t[1], t[2]);
	glScalef(fScale, fScale, fScale);
	glTranslatef(-t[0], -t[1], -t[2]);
	glMultMatrixd(m_modelViewMatrixGlobal);
	glGetDoublev(GL_MODELVIEW_MATRIX, m_modelViewMatrixGlobal);

	updateGL();
}

void MyGLWidget::rotate(const double* axis, double angle)
{

	double t[3] = { m_modelViewMatrixGlobal[12], m_modelViewMatrixGlobal[13], m_modelViewMatrixGlobal[14] };//得到x、y、z轴的平移量

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(t[0], t[1], t[2]);
	glRotated(angle, axis[0], axis[1], axis[2]);
	glTranslatef(-t[0], -t[1], -t[2]);
	glMultMatrixd(m_modelViewMatrixGlobal);
	glGetDoublev(GL_MODELVIEW_MATRIX, m_modelViewMatrixGlobal);
}

// move screen from pt1 to pt2
void MyGLWidget::move(QPoint pt1, QPoint pt2){
	int x = pt2.x() - pt1.x();
	int y = pt2.y() - pt1.y();
	double fX = x*1.0 / _nWidth;
	double fY = -y*1.0 / _nHeight;

	double t[3] = { m_modelViewMatrixGlobal[12], m_modelViewMatrixGlobal[13], m_modelViewMatrixGlobal[14] };//得到x、y、z轴的平移量

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(t[0], t[1], t[2]);
	glTranslatef(fX, fY,0.0f);
	glTranslatef(-t[0], -t[1], -t[2]);
	glMultMatrixd(m_modelViewMatrixGlobal);
	glGetDoublev(GL_MODELVIEW_MATRIX, m_modelViewMatrixGlobal);
}

void MyGLWidget::updateTrackBallPos(QPoint pt, double* result)
{

	double x = this->width() / 2 - pt.x();
	double y = pt.y() - this->height() / 2;
	if (x > m_dbRadius) x = m_dbRadius;
	if (y > m_dbRadius) y = m_dbRadius;
	double x2y2 = x*x + y*y;

	result[0] = x;
	result[1] = y;
	result[2] = sqrt(m_dbRadius*m_dbRadius - x2y2);

}

void MyGLWidget::SetDataT(const double* data){
	_pDataT = data;

}

void MyGLWidget::SetLabels(const int* labels){
	_arrLabels = labels;
}

void MyGLWidget::SetDataB(const double* data){
	_pDataB = data;
}

void MyGLWidget::SetModelE(MeteModel* pModelE){
	_pModelE = pModelE;
}

void MyGLWidget::SetModelT(MeteModel* pModelT) {
	_pModelT = pModelT;
}

void MyGLWidget::SetVar(const double* pVM, const double *pMV){
	_gridDataMeanVarB = pMV;
	_gridDataVarMeanB = pVM;
}

void MyGLWidget::SetContourTruth(QList<ContourLine>* listContourTruth){
	for (int i = 0; i < g_temperatureLen; i++)
	{
		_listContourTruth[i] = listContourTruth[i];
	}
}

void MyGLWidget::SetContourMean(QList<ContourLine>* listContour){
	for (int i = 0; i < g_temperatureLen; i++)
	{
		_listContourMean[i] = listContour[i];
	}

}

void MyGLWidget::SetContourIntervalB(QList<ContourLine>* listContourMin, QList<ContourLine>* listContourMax){
	for (int i = 0; i < g_temperatureLen; i++)
	{
		_listContourMaxB[i] = listContourMax[i];
		_listContourMinB[i] = listContourMin[i];
	}
}

void MyGLWidget::SetMultiStatistic(QList<ContourLine>* listContourMin, QList<ContourLine>* listContourMax){
	for (int i = 0; i < g_gradient_l;i++)
	{
		_listContourMin_3[i] = listContourMin[i];
		_listContourMax_3[i] = listContourMax[i];
	}
}

void MyGLWidget::generateContour(){
	// data preparation
	for (int i = 0; i < g_gradient_l;i++)
	{
// 		generateContourImp_New(_listContourMin_3[i], _listContourMax_3[i], _listUnionAreaB[i]);
	}
	for (int i = 0; i < g_temperatureLen;i++)
	{
// 		generateContourImp_New(_listContourMinE[i], _listContourMaxE[i], _listUnionAreaE[i]);
	}

	for (int i = 0; i < g_temperatureLen; i++)
	{
// 		generateContourImp_New(_listContourMinB[i], _listContourMaxB[i], _listIntersectionAreaB[i]);
	}
}

bool MyGLWidget::processHits(GLint hits, GLuint buffer[], PickIndex& pick)
{
	if (hits > 0)
	{
		GLuint* ptr = (GLuint*)buffer;
		double fLastDepth = 0;
		if (hits>0)
		{
			// name
			GLuint names = *ptr;
			ptr++;
			PickIndex pickI;
			// depth1
			double depth1 = (double)*ptr;
			ptr++;
			// depth2
			double depth2 = (double)*ptr;
			ptr++;
			// x
			pickI._nX = *ptr;
			ptr++;
			// y
			pickI._nY = *ptr;
			ptr++;
			pick = pickI;
		}
		return true;
	}
	else
	{
		return false;
	}
	return false;
}

void MyGLWidget::drawPlaceHolder()
{
	for (int i = 0; i < g_focus_w;i++)
	{
		glLoadName(i);
		for (int j = 0; j < g_focus_h; j++){
			glPushName(j);
			double fX = _fLeft + i*_fScaleW;
			double fY = _fBottom + j*_fScaleH;
			glBegin(GL_QUADS);
			glVertex2f(fX - 0.5*_fScaleW, fY - 0.5*_fScaleH);
			glVertex2f(fX + 0.5*_fScaleW, fY - 0.5*_fScaleH);
			glVertex2f(fX + 0.5*_fScaleW, fY + 0.5*_fScaleH);
			glVertex2f(fX - 0.5*_fScaleW, fY + 0.5*_fScaleH);
			glEnd();
			glPopName();
		}
	}

	// for variance chart
	int nWidth = 100;
	double dbStep = _fRight / nWidth;

	for (size_t i = 0; i < nWidth; i++)
	{
		glLoadName(i + g_focus_w);
		double fX0 = 0 + i*dbStep;
		double fX1 = 0 + (i+1)*dbStep;
		double fY = _fTop;
		glBegin(GL_QUADS);
		glVertex2f(fX0, fY );
		glVertex2f(fX1, fY );
		glVertex2f(fX1, fY + 10);
		glVertex2f(fX0, fY + 10);
		glEnd();
	}
}

void MyGLWidget::viewShowGrid(bool on){
	_bShowGridLines = on;
	updateGL();
}

void MyGLWidget::viewShowBackground(bool on){
	_displayStates._bShowBackground = on; updateGL();
};
void MyGLWidget::onCheckShowBeliefEllipse(bool bChecked) {
	_displayStates._bShowBeliefEllipse = bChecked; updateGL();
}
void MyGLWidget::viewShowIntersection(bool on){
	_bShowIntersection = on; updateGL(); };
void MyGLWidget::viewShowUnionB(bool on){ 
	_bShowUnionB = on; updateGL();
};
void MyGLWidget::viewShowUnionE(bool on){
	_displayStates._bShowUnionE = on; 
	updateGL();
};
void MyGLWidget::viewShowGradientE(bool on) {
	_displayStates._bShowGradientE = on; 
	updateGL();
};
void MyGLWidget::viewShowLineChart(bool on){
//	_bShowLineChart = on;
	updateGL(); 
};
void MyGLWidget::viewShowContourLineTruth(bool on){
	_bShowContourLineTruth = on; updateGL(); };
void MyGLWidget::viewShowContourLine(bool on){
	_displayStates._bShowContourLine = on; updateGL(); };
void MyGLWidget::viewShowContourLineMin(bool on){
	_displayStates._bShowContourLineMin = on; updateGL();
};
void MyGLWidget::viewShowContourLineMax(bool on){
	_displayStates._bShowContourLineMax = on; updateGL();
};
void MyGLWidget::viewShowContourLineMean(bool on){ 
	_displayStates._bShowContourLineMean = on; updateGL();
};

void MyGLWidget::viewShowClusterBS(bool on){
	_bShowClusterBS = on; updateGL();
}
void MyGLWidget::viewShowClusterBV(bool on){
	_bShowClusterBV = on; updateGL();
}

void MyGLWidget::DrawText(char* pText, double fX, double fY){
	font.PrintText(pText, fX, fY);

}

// reload texture
void MyGLWidget::ReloadTexture() {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->ReloadTexture();
	}
	updateGL();
}


void MyGLWidget::updateMinPts(int minPts) {
	qDebug() << "set minpts";
	if(_pModelE)
		_pModelE->SetMinPts(minPts);

	ReloadTexture();
}

void MyGLWidget::updateEps(double eps) {
	qDebug() << "set eps";
	if (_pModelE)
		_pModelE->SetEps(eps);

	ReloadTexture();
}