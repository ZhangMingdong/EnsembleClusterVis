#include "mychartwidget.h"
#include <gl/GLU.h>



#include <iostream>


#include "MeteModel.h"
#include "LayerLayout.h"
#include "DataField.h"
#include "ColorMap.h"




using namespace std;


MyChartWidget::MyChartWidget(QWidget *parent)
	: MyGLWidget(parent)
{
	m_dbScale = 20;
}

MyChartWidget::~MyChartWidget()
{
}
void SetClusterColor(int nIndex, double dbOpacity = .8) {
	glColor4f(ColorMap::GetCategory10D(nIndex, 0)
		, ColorMap::GetCategory10D(nIndex, 1)
		, ColorMap::GetCategory10D(nIndex, 2)
		, dbOpacity);
}
void MyChartWidget::paint() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



	glLineWidth(1.0f);
	glColor3f(1.0, 0, 0);

	glPushMatrix();

	//	
	/*
	glScaled(.1, .1, 1);
	glTranslatef(-5, -5, 0);
	glBegin(GL_POLYGON);
	glVertex3f(0, 0,0);
	glVertex3f(0, 10,0);
	glVertex3f(10, 10,0);
	glVertex3f(10, 0,0);
	glEnd();
	*/
	/*
	// rectangle
	glColor3f(0, 0, 0);
	glBegin(GL_LINE_LOOP);
	glVertex2d(-10, -10);
	glVertex2d(-10, 10);
	glVertex2d(10, 10);
	glVertex2d(10, -10);	
	glEnd();
	*/

	// axis
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(0, -10000);
	glVertex2d(0, 10000);
	glVertex2d(-10000, 0);
	glVertex2d(10000, 0);
	glEnd();
	

	glPointSize(5);
	glBegin(GL_POINTS);

	int nEnsemble = _pModelE->GetEnsembleLen();
	for (int i = 0; i < nEnsemble; i++)
	{
		SetClusterColor(_pModelE->GetLabel(i));
		glVertex2d(_pModelE->GetPC(i, 0), _pModelE->GetPC(i, 1));

	}
	glEnd();



	glPopMatrix();
}

void MyChartWidget::init() {

}

void MyChartWidget::SetModelE(MeteModel* pModelE) {	
	_pModelE = pModelE;	
}

void MyChartWidget::mouseDoubleClickEvent(QMouseEvent *event) {

}
