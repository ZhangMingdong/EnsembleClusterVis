#include "myglwidget.h"
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

MyGLWidget::~MyGLWidget()
{
	// release the layers
	for each (MeteLayer* layer in _vecLayers)
	{
		delete layer;
	}

	delete _pLayout;
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
//	pLayer->InitLayer(_pLayout->_dbLeft, _pLayout->_dbRight, _pLayout->_dbTop, _pLayout->_dbBottom, _fScaleW, _fScaleH);
//	_vecLayers.push_back(pLayer);	

	MeteLayer* pLayer = new EnsembleLayer();
	pLayer->SetModel(_pModelE);
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);

	pLayer = new ClusterLayer();
	pLayer->SetModel(_pModelE);
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);

	// create cost line layer
	pLayer = new CostLineLayer();
	pLayer->InitLayer(_pLayout, _fScaleW, _fScaleH);
	_vecLayers.push_back(pLayer);


	_font.InitFont(wglGetCurrentDC(), L"Arial", 22);
	MeteLayer::_pCB = this;
}

void MyGLWidget::paintGL(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw selected point
	glLineWidth(2.0f);
	glColor4f(1.0, .3, 0.0, 1.0);
	double fX0 = _pLayout->_dbLeft + _nSelectedLeft*_fScaleW;
	double fY0 = _pLayout->_dbBottom + _nSelectedBottom*_fScaleH;
	double fX1 = _pLayout->_dbLeft + _nSelectedRight*_fScaleW;
	double fY1 = _pLayout->_dbBottom + _nSelectedTop*_fScaleH;
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

	// grid lines
	if (_bShowGridLines){
		
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

void MyGLWidget::drawContourLine(const QList<ContourLine>& contours){

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

void MyGLWidget::SetModelE(MeteModel* pModelE){
	_pModelE = pModelE;
}

void MyGLWidget::SetModelT(MeteModel* pModelT) {
	_pModelT = pModelT;
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
	_font.PrintText(pText, fX, fY);

}

// reload texture
void MyGLWidget::onTextureReloaded() {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->ReloadTexture();
	}

	updateGL();
}
void MyGLWidget::setUncertaintyAreas(int nUCAreas) {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->SetUncertaintyAreas(nUCAreas);
	}

	updateGL();
}

void MyGLWidget::updateMinPts(int minPts) {
	qDebug() << "set minpts";
	if(_pModelE)
		_pModelE->SetMinPts(minPts);

	onTextureReloaded();
}

void MyGLWidget::updateEps(double eps) {
	qDebug() << "set eps";
	if (_pModelE)
		_pModelE->SetEps(eps);

	onTextureReloaded();
}
void MyGLWidget::updateVarSmooth(int nSmooth) {
	if (_pModelE)
		_pModelE->SetSmooth(nSmooth);
	onTextureReloaded();
}
void MyGLWidget::updateBgFunction(int nBgFunction) {
	if (_pModelE)
		_pModelE->SetBgFunction((MeteModel::enumBackgroundFunction)nBgFunction);
	onTextureReloaded();
}

void MyGLWidget::updateUncertaintyAreas(int nAreas) {
	if (_pModelE)
		_pModelE->SetUncertaintyAreas(nAreas);
	setUncertaintyAreas(nAreas);
}

void MyGLWidget::updateFocusedCluster(int nFocusedCluster) {
	for each (MeteLayer* pLayer in _vecLayers)
	{
		pLayer->SetFocusedCluster(nFocusedCluster);
	}

	updateGL();
}

void MyGLWidget::updateFocusedRegion(int nFocusedRegion) {
//	qDebug() << "updateFocusedRegion";
	_pModelE->SetFocusedRegion(nFocusedRegion);

	updateGL();
}