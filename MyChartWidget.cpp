#include "mychartwidget.h"
#include <gl/GLU.h>



#include <iostream>


#include "MeteLayer.h"
#include "CostLineLayer.h"
#include "EnsembleLayer.h"
#include "ClusterLayer.h"
#include "MeteModel.h"
#include "VarAnalysis.h"
#include "LayerLayout.h"
#include "DataField.h"
#include "ColorMap.h"

#include "Sequence2D.h"

#define BUFSIZE 512


using namespace std;


MyChartWidget::MyChartWidget(QWidget *parent)
	: MyGLWidget(parent)
{
	_pSequence = new Sequence2D();
}

MyChartWidget::~MyChartWidget()
{
	delete _pSequence;
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

	double dbMin = _pSequence->GetMin();
	double dbMax = _pSequence->GetMax();
	int nLen = _pSequence->GetLength();
	glScaled(.05, .03, 1);
	glTranslatef(-nLen / 2.0, -(dbMin+dbMax) / 2.0, 0);

	// draw spaghetti
	drawSpaghetti();

	// draw grid lines
	drawGridLines();

	// draw selected group
//	drawSelectedGroup();
	drawSelectedDBGroup();

	// draww abstracted groups
//	drawGroups();
	drawDBGroups();

	glPopMatrix();
}

void MyChartWidget::init() {

}

void MyChartWidget::SetModelE(MeteModel* pModelE) {	
	// set model
	_pModelE = pModelE;

	// generate the sequences
	generateSequences();

	// detect the trends
//	_pSequence->TrendDetectRootOnly();
//	_pSequence->TrendDetect();
//	_pSequence->TrendDetectDB();
	_pSequence->TrendDetectDBImproved();
}

void MyChartWidget::mouseDoubleClickEvent(QMouseEvent *event) {
	_nCurrentGroup++;
	if (_pSequence->GetGroupSize())
	{
		_nCurrentGroup %= _pSequence->GetGroupSize();
	}
	else {
		_nCurrentGroup %= _pSequence->GetDBGroupSize();
	}
	updateGL();
}

void MyChartWidget::generateSequences() {
	int nWidth = _pModelE->GetW();
	int nHeight = _pModelE->GetH();
	int nEns = _pModelE->GetEnsembleLen();
	nEns = 30;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < nWidth; j++)
		{
#ifdef REVERSE
			double dbValue = _pModelE->GetData()->GetData(i, nHeight - 1, nWidth - 1 - j);
#else
			double dbValue = _pModelE->GetData()->GetData(i, nHeight - 1, j);
#endif
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;

			seq.push_back(dbValue);
		}
		//		if(i==0|| i == 1 || i == 2 || i == 3 || i == 4 || i==5||i==8)
		_pSequence->AddSequence(seq);
	}

//	double dbEpsilon = (dbMax - dbMin) *.03;

	double dbEpsilon = 1.3;

	_pSequence->Init(nWidth, dbMin, dbMax, dbEpsilon);
}

void MyChartWidget::drawGroups() {
	for (size_t i = 0; i < _pSequence->GetGroupSize(); i++)
	{
		Group g = _pSequence->GetGroup(i);
		glColor4f(0, 1, 1, .5);
		double dbSize = g._member.size();
		glLineWidth(dbSize);

		glBegin(GL_LINE_STRIP);
		for (size_t k = g._nS; k <= g._nE; k++)
		{
			double dbY = 0;
			for (size_t j = 0; j < g._member.size(); j++)
			{
				dbY += _pSequence->GetValue(g._member[j], k);
			}

			glVertex3f(k, dbY / dbSize, 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawSelectedGroup() {

	glColor3f(0, 1, 0);
	glLineWidth(2.0f);
	glColor3f(0, 1, 0);
	Group g = _pSequence->GetGroup(_nCurrentGroup);
//	cout << g._nS << "\t" << g._nE << endl;
	for (size_t j = 0; j < g._member.size(); j++)
	{
//		cout << g._member[j] << endl;
		glBegin(GL_LINE_STRIP);
		for (size_t k = g._nS; k <= g._nE; k++)
		{
			glVertex3f(k, _pSequence->GetValue(g._member[j],k), 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawGridLines() {
	int nLen = _pSequence->GetLength();
	double dbMin = _pSequence->GetMin();
	double dbMax = _pSequence->GetMax();
	double dbEpsilon = _pSequence->GetEpsilon();

	glColor3f(0, 0, .5);
	glBegin(GL_LINES);

	// horizontal lines
	double dbY = dbMin;
	while (dbY<dbMax+.01)
	{
		glVertex2d(0, dbY);
		glVertex2d(nLen - 1, dbY);

		dbY += dbEpsilon;
	}
	
	// vertical lines
	for (size_t i = 0; i < nLen; i += 10)
	{
		glVertex2d(i, _pSequence->GetMin());
		glVertex2d(i, _pSequence->GetMax());
	}
	glColor4f(0, 0, 0,.2);
	// vertical lines
	for (size_t i = 0; i < nLen; i ++)
	{
		glVertex2d(i, _pSequence->GetMin());
		glVertex2d(i, _pSequence->GetMax());
	}
	glEnd();
}

void MyChartWidget::drawSpaghetti() {
	int nEns = _pSequence->GetEns();
	int nLen = _pSequence->GetLength();
	for (size_t i = 0; i < nEns; i++)
	{
		glBegin(GL_LINE_STRIP);
		for (size_t j = 0; j < nLen; j++)
		{
			glVertex3f(j, _pSequence->GetValue(i, j), 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawDBGroups() {
	for (size_t i = 0; i < _pSequence->GetDBGroupSize(); i++)
	{
		DBGroup g = _pSequence->GetDBGroup(i);
		int nS = g._dbS;
		int nE = g._dbE;
		glColor4f(0, 1, 1, .5);
		double dbSize = g._member.size();
		glLineWidth(dbSize);


		bool bFirstSeg = false;
		if (g._dbS - nS >.001) {
			bFirstSeg = true;
			nS++;
		}

		glBegin(GL_LINE_STRIP);
		// first segment		
		if (bFirstSeg)
		{
			double dbY = 0;
			for (size_t j = 0; j < g._member.size(); j++)
			{
				dbY += _pSequence->GetDBValue(g._dbS, g._member[j]);
			}

			glVertex3f(g._dbS, dbY / dbSize, 0);
		}
		// the segments of full timesteps
		for (size_t k = nS; k <= nE; k++)
		{
			double dbY = 0;
			for (size_t j = 0; j < g._member.size(); j++)
			{
				dbY += _pSequence->GetValue(g._member[j], k);
			}

			glVertex3f(k, dbY / dbSize, 0);
		}

		// last segment
		if (g._dbE - nE >.001) {
			double dbY = 0;
			for (size_t j = 0; j < g._member.size(); j++)
			{
				dbY += _pSequence->GetDBValue(g._dbE, g._member[j]);
			}

			glVertex3f(g._dbE, dbY / dbSize, 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawSelectedDBGroup() {
	glColor3f(0, 1, 0);
	glLineWidth(2.0f);
	glColor3f(0, 1, 0);
	DBGroup g = _pSequence->GetDBGroup(_nCurrentGroup);
	int nS = g._dbS;
	int nE = g._dbE;
	bool bFirstSeg = false;

	if (g._dbS - nS >.001) {
		bFirstSeg = true;
		nS++;
	}

	for (size_t j = 0; j < g._member.size(); j++)
	{
		//		cout << g._member[j] << endl;

		glBegin(GL_LINE_STRIP);

		// first segment
		if (bFirstSeg)
		{
			glVertex3f(g._dbS, _pSequence->GetDBValue(g._dbS, g._member[j]), 0);
		}
		// the segments of full timesteps
		for (size_t k = nS; k <= nE; k++)
		{
			glVertex3f(k, _pSequence->GetValue(g._member[j], k), 0);
		}
		// last segment
		if (g._dbE - nE >.001) {
			glVertex3f(g._dbE, _pSequence->GetDBValue(g._dbE, g._member[j]), 0);
		}

		glEnd();
	}
}