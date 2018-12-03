#include "MyMapWidget.h"
#include <gl/GLU.h>

#include <QMouseEvent>
#include <QDebug>


#include <iostream>


#include "MeteLayer.h"
#include "CostLineLayer.h"
#include "EnsembleLayer.h"
#include "ClusterLayer.h"
#include "MeteModel.h"
#include "VarAnalysis.h"
#include "LayerLayout.h"

#define BUFSIZE 512

using namespace std;



MyMapWidget::MyMapWidget(QWidget *parent)
: MyGLWidget(parent)
, _bShowGridLines(false)
, _bShowIntersection(false)
, _bShowUnionB(false)
, _bShowClusterBS(false)
, _bShowClusterBV(false)
, _bShowLineChart(false)
, _bShowContourLineTruth(false)
{
	_pLayout = new LayerLayout();
	
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

MyMapWidget::~MyMapWidget()
{
	// release the layers
	for each (MeteLayer* layer in _vecLayers)
	{
		delete layer;
	}

	delete _pLayout;
}

void MyMapWidget::init()
{
	MeteLayer* pLayer = new EnsembleLayer();
	pLayer->SetModel(_pModelE);
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);

	/*
	pLayer = new ClusterLayer();
	pLayer->SetModel(_pModelE);
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);
	*/
	// create cost line layer
	pLayer = new CostLineLayer();
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);


	_font.InitFont(wglGetCurrentDC(), L"Arial", 22);
	MeteLayer::_pCB = this;
}


void MyMapWidget::timerEvent(QTimerEvent* event)
{

}

void MyMapWidget::mousePressEvent(QMouseEvent * event)
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

void MyMapWidget::mouseReleaseEvent(QMouseEvent * event)
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

void MyMapWidget::mouseMoveEvent(QMouseEvent *event)
{
	if (m_nMouseState==2)
	{
		select(_nSelectedRight, _nSelectedTop, event->pos());
	}
	else if (m_nMouseState==1){	// left button down

		double scale = GetSWScale();
		m_pt3Eye.x -= (event->pos().x() - _ptLast.x()); //cx*scale;//(x*cos(m_rotation.z)-y*sin(m_rotation.z));
		m_pt3Eye.y += (event->pos().y() - _ptLast.y()); //cy*scale;//(y*cos(m_rotation.z)+x*sin(m_rotation.z));

		_ptLast = event->pos();
	}
	updateGL();
}

void MyMapWidget::mouseDoubleClickEvent(QMouseEvent *event){
	select(_nSelectedX, _nSelectedRight, event->pos());
}

void MyMapWidget::drawContourLine(const QList<ContourLine>& contours) {

	for each (ContourLine contour in contours)
	{
		glBegin(GL_LINE_STRIP);
		for each (QPointF pt in contour._listPt)
		{
			double x = _pLayout->_dbLeft + pt.x() * _fScaleW;
			double y = _pLayout->_dbTop - pt.y() * _fScaleH;
			glVertex2f(x, y);
		}
		glEnd();
	}
}

void MyMapWidget::paint() {

	// draw selected point
	glLineWidth(2.0f);
	glColor4f(1.0, .3, 0.0, 1.0);
	double fX0 = _pLayout->_dbLeft + _nSelectedLeft*_fScaleW;
	double fY0 = _pLayout->_dbBottom + _nSelectedBottom*_fScaleH;
	double fX1 = _pLayout->_dbLeft + _nSelectedRight*_fScaleW;
	double fY1 = _pLayout->_dbBottom + _nSelectedTop*_fScaleH;
	glBegin(GL_LINE_LOOP);
	glVertex2f(fX0, fY0);
	glVertex2f(fX0, fY1);
	glVertex2f(fX1, fY1);
	glVertex2f(fX1, fY0);
	glEnd();

	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->DrawLayer(_displayStates);
	}

	// grid lines
	if (_bShowGridLines) {

		glColor4f(.5, .5, .5, .5);
		glBegin(GL_LINES);
		int nStep = 10;
		for (int i = 0; i < g_globalW; i += nStep)
		{
			glVertex2f(_pLayout->_dbLeft + i*_fScaleW, _pLayout->_dbBottom);
			glVertex2f(_pLayout->_dbLeft + i*_fScaleW, _pLayout->_dbTop);
		}
		for (int j = 0; j < g_globalH; j += nStep)
		{
			glVertex2f(_pLayout->_dbLeft, _pLayout->_dbBottom + j*_fScaleH);
			glVertex2f(_pLayout->_dbRight, _pLayout->_dbBottom + j*_fScaleH);
		}
		glEnd();
		// show detail
		if (false) {
			glColor4f(.5, .5, .5, .3);
			glBegin(GL_LINES);
			int nStep = 1;
			for (int i = 0; i < g_globalW; i += nStep)
			{
				glVertex2f(_pLayout->_dbLeft + i*_fScaleW, _pLayout->_dbBottom);
				glVertex2f(_pLayout->_dbLeft + i*_fScaleW, _pLayout->_dbTop);
			}
			for (int j = 0; j < g_globalH; j += nStep)
			{
				glVertex2f(_pLayout->_dbLeft, _pLayout->_dbBottom + j*_fScaleH);
				glVertex2f(_pLayout->_dbRight, _pLayout->_dbBottom + j*_fScaleH);
			}
			glEnd();
		}
	}
}

void MyMapWidget::select(int& nX, int&nY, const QPoint& pt) {


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
			qDebug() << pick._nX << "\t" << pick._nY;
		}
		else {
			// bias of the index of the variance is beyond g_focus_w
			int nIndex = pick._nX - g_focus_w;
			for each (MeteLayer* pLayer in _vecLayers)
			{
				pLayer->OnSelectVar(nIndex);
			}
		}

		updateGL();
	}
}

void MyMapWidget::SetModelE(MeteModel* pModelE){
	_pModelE = pModelE;
}

bool MyMapWidget::processHits(GLint hits, GLuint buffer[], PickIndex& pick)
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

void MyMapWidget::drawPlaceHolder()
{
	for (int i = 0; i < g_focus_w;i++)
	{
		glLoadName(i);
		for (int j = 0; j < g_focus_h; j++){
			glPushName(j);
			double fX = _pLayout->_dbLeft + i*_fScaleW;
			double fY = _pLayout->_dbBottom + j*_fScaleH;
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
	int nWidth = VarAnalysis::_nVarBars;
	double dbStep = (_pLayout->_dbRight-_pLayout->_dbLeft) / nWidth;

	for (size_t i = 0; i < nWidth; i++)
	{
		glLoadName(i + g_focus_w);
		double fX0 = _pLayout->_dbLeft + i*dbStep;
		double fX1 = _pLayout->_dbLeft + (i+1)*dbStep;
		double fY = _pLayout->_dbTop;
		glBegin(GL_QUADS);
		glVertex2f(fX0, fY );
		glVertex2f(fX1, fY );
		glVertex2f(fX1, fY + 10);
		glVertex2f(fX0, fY + 10);
		glEnd();
	}
}

void MyMapWidget::viewShowGrid(bool on){
	_bShowGridLines = on;
	updateGL();
}

void MyMapWidget::viewShowBackground(bool on){
	_displayStates._bShowBackground = on; updateGL();
};

void MyMapWidget::onCheckShowBeliefEllipse(bool bChecked) {
	_displayStates._bShowBeliefEllipse = bChecked; updateGL();
}

void MyMapWidget::viewShowIntersection(bool on){
	_bShowIntersection = on; updateGL(); };

void MyMapWidget::viewShowUnionB(bool on){ 
	_bShowUnionB = on; updateGL();
};

void MyMapWidget::viewShowUnionE(bool on){
	_displayStates._bShowUnionE = on; 
	updateGL();
};

void MyMapWidget::viewShowGradientE(bool on) {
	_displayStates._bShowGradientE = on; 
	updateGL();
};

void MyMapWidget::viewShowLineChart(bool on){
//	_bShowLineChart = on;
	updateGL(); 
};

void MyMapWidget::viewShowContourLineTruth(bool on){
	_bShowContourLineTruth = on; updateGL(); };

void MyMapWidget::viewShowContourLine(bool on){
	_displayStates._bShowContourLine = on; updateGL(); };

void MyMapWidget::viewShowContourLineSorted(bool on) {
	_displayStates._bShowContourLineSorted = on; updateGL();
};

void MyMapWidget::viewShowContourLineSortedSDF(bool on) {
	_displayStates._bShowContourLineSortedSDF = on; updateGL();
};
void MyMapWidget::viewShowContourLineSDF(bool on) {
	_displayStates._bShowContourLineSDF = on; updateGL();
};
void MyMapWidget::viewShowContourLineMin(bool on){
	_displayStates._bShowContourLineMin = on; updateGL();
};

void MyMapWidget::viewShowContourLineMax(bool on){
	_displayStates._bShowContourLineMax = on; updateGL();
};

void MyMapWidget::viewShowContourLineMean(bool on){ 
	_displayStates._bShowContourLineMean = on; updateGL();
};

void MyMapWidget::viewShowClusterBS(bool on){
	_bShowClusterBS = on; updateGL();
}

void MyMapWidget::viewShowClusterBV(bool on){
	_bShowClusterBV = on; updateGL();
}

void MyMapWidget::DrawText(char* pText, double fX, double fY){
	_font.PrintText(pText, fX, fY);

}

// reload texture
void MyMapWidget::onTextureReloaded() {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->ReloadTexture();
	}

	updateGL();
}

void MyMapWidget::setUncertaintyAreas(int nUCAreas) {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->SetUncertaintyAreas(nUCAreas);
	}

	updateGL();
}

void MyMapWidget::updateMinPts(int minPts) {
	qDebug() << "set minpts";
	/*
	if(_pModelE)
		_pModelE->SetMinPts(minPts);

	onTextureReloaded();
	*/
}

void MyMapWidget::updateEps(double eps) {
	qDebug() << "set eps";
	/*
	if (_pModelE)
		_pModelE->SetEps(eps);

	onTextureReloaded(); 
	*/
}

void MyMapWidget::updateVarSmooth(int nSmooth) {
	if (_pModelE)
		_pModelE->SetSmooth(nSmooth);
	onTextureReloaded();
}

void MyMapWidget::updateBgFunction(int nBgFunction) {
	if (_pModelE)
		_pModelE->SetBgFunction((MeteModel::enumBackgroundFunction)nBgFunction);
	onTextureReloaded();
}

void MyMapWidget::updateUncertaintyAreas(int nAreas) {
	if (_pModelE)
		_pModelE->SetUncertaintyAreas(nAreas);
	setUncertaintyAreas(nAreas);
}

void MyMapWidget::updateFocusedCluster(int nFocusedCluster) {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->SetFocusedCluster(nFocusedCluster);
	}

	updateGL();
}

void MyMapWidget::updateFocusedRegion(int nFocusedRegion) {
//	qDebug() << "updateFocusedRegion";
	_pModelE->SetFocusedRegion(nFocusedRegion);

	updateGL();
}

void MyMapWidget::updateEOF(int nEOF) {
	if (_pModelE)
		_pModelE->SetEOF(nEOF);
	onTextureReloaded();
}

void MyMapWidget::updateMember(int nMember) {
	if (_pModelE)
		_pModelE->SetMember(nMember);
	onTextureReloaded();
}

void MyMapWidget::updateEnsCluster(int nEnsCluster) {
	if (_pModelE)
		_pModelE->SetEnsCluster(nEnsCluster);
	onTextureReloaded();
}


void MyMapWidget::updateContourLevel(int nLevel) {
	if (_pModelE)
		_pModelE->SetContourLevel(nLevel);
	onTextureReloaded();
}