#include "MyGLWidget.h"

#include <gl/GLU.h>

#include <QMouseEvent>

MyGLWidget::MyGLWidget(QWidget *parent):QGLWidget(parent)
, m_dbRadius(1000)
, m_nMouseState(0)
, c_dbPI(3.14159265)
, _nSelectedX(10)
, _nSelectedY(10)
, _nSelectedLeft(-1)
, _nSelectedRight(-1)
, _nSelectedTop(-1)
, _nSelectedBottom(-1)
, m_pt3Eye(0, 0, 300)
, m_dbFovy(45)
, m_dbZNear(0.0001)
, m_dbZFar(1500.0)
{
}


MyGLWidget::~MyGLWidget()
{
}


void MyGLWidget::paintGL(){
	//1.viewport setting
	glViewport(0, 0, _nWidth, _nHeight);
	//2.projection mode setting
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-_nWidth / 2.0, _nWidth / 2.0, -_nHeight / 2.0, _nHeight / 2.0, m_dbZNear, m_dbZFar);
	// 	gluPerspective(m_dbFovy, (m_nHeight == 0) ? 1.0 : (double)m_nWidth / (double)m_nHeight, m_dbZNear, m_dbZFar);

	//3.clear background
	glDrawBuffer(GL_BACK);
	// 	glClearColor(m_clearColor.r, m_clearColor.g, m_clearColor.b, m_clearColor.a);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//4.Model View Setting
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//4.1 View
	gluLookAt(m_pt3Eye.x, m_pt3Eye.y, m_pt3Eye.z, m_pt3Eye.x, m_pt3Eye.y, 0, 0, 1, 0);

	glScaled(m_dbScale, m_dbScale, 1);

	paint();
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

	init();




}

void MyGLWidget::resizeGL(int width, int height) {
	_nWidth = width;
	_nHeight = height;
	// 1.viewport
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);

	// 2.projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(30.0, (GLfloat)width / (GLfloat)height, .1, 100.0);
	
}

void MyGLWidget::timerEvent(QTimerEvent* event)
{

}

void MyGLWidget::mousePressEvent(QMouseEvent * event)
{
	if (event->button() == Qt::MidButton)
	{
	}
	else if (event->button() == Qt::LeftButton)
	{
		updateTrackBallPos(event->pos(), m_dbTrackBallPos);
		m_nMouseState = 1;
		_ptLast = event->pos();
	}
	updateGL();
}

void MyGLWidget::mouseReleaseEvent(QMouseEvent * event)
{

	m_nMouseState = 0;
	updateGL();

}

void MyGLWidget::mouseMoveEvent(QMouseEvent *event)
{
	if (m_nMouseState == 1) {	// left button down

		double scale = GetSWScale();
		m_pt3Eye.x -= (event->pos().x() - _ptLast.x()); //cx*scale;//(x*cos(m_rotation.z)-y*sin(m_rotation.z));
		m_pt3Eye.y += (event->pos().y() - _ptLast.y()); //cy*scale;//(y*cos(m_rotation.z)+x*sin(m_rotation.z));

		_ptLast = event->pos();
	}
	updateGL();
}

void MyGLWidget::mouseDoubleClickEvent(QMouseEvent *event) {
}

void MyGLWidget::wheelEvent(QWheelEvent * event)
{
	if (event->delta()>0)
	{
		m_dbScale *= 1.1;
	}
	else
	{
		m_dbScale *= .9;
	}
	updateGL();
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