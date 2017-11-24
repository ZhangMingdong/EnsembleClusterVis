#include "mychartwidget.h"
#include <gl/GLU.h>



#include <iostream>


#include "MeteModel.h"
#include "LayerLayout.h"
#include "DataField.h"
#include "ColorMap.h"

#include "Sequence2D.h"

#include <qfile.h>
#include <qtextstream.h>
#include <random>

#define BUFSIZE 512

//#define USE_ARTIFICIAL
//#define USE_GDP


#ifdef USE_ARTIFICIAL
const double g_dbScaleW = 2;
const double g_dbScaleH = 1;
const double g_dbEpsilon = 5;
#else
const double g_dbScaleW = .05;
const double g_dbScaleH = .03;
const double g_dbEpsilon = 1;
#endif



using namespace std;


MyChartWidget::MyChartWidget(QWidget *parent)
	: MyGLWidget(parent)
{
//	_pSequence = Sequence2D::GenerateInstance(Sequence2D::ST_Buchin);
//	_pSequence = Sequence2D::GenerateInstance(Sequence2D::ST_Van);
	_pSequence = Sequence2D::GenerateInstance(Sequence2D::ST_Jeung);
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
	glScaled(g_dbScaleW, g_dbScaleH, 1);
	glTranslatef(-nLen / 2.0, -(dbMin+dbMax) / 2.0, 0);

	// draw spaghetti
	drawSpaghetti();

	// draw grid lines
	drawGridLines();

	// draw selected group
	drawSelectedGroup();
	drawSelectedDBGroup();

	// draww abstracted groups
	drawGroups();
	drawDBGroups();

	// draw the event line
	drawEvents();

	glPopMatrix();
}

void MyChartWidget::init() {

}

void MyChartWidget::SetModelE(MeteModel* pModelE) {	
#ifdef USE_ARTIFICIAL
//	generateSequenceArtificial1();
//	generateSequenceArtificial2();
	generateSequenceArtificial3();
#else	
#ifdef USE_GDP	
	// use gdp data
	QFile file("data\\Industry (p of GDP).csv"); // this is a name of a file text1.txt sent from main method
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		return;
	}
	QTextStream in(&file);
	QString line = in.readLine();
	vector<vector<double>> vecData;
	while (!line.isNull())
	{
		vector<double> seq;
		QStringList list = line.split(",");
		for (size_t i = 0, length = list.size() - 1; i < length; i++)
		{
			if (list[i].isEmpty()) seq.push_back(-1);
			else seq.push_back(list[i].toDouble());
		}
		vecData.push_back(seq);
		line = in.readLine();
	}

	generateSequences(vecData);
#else 

	// use ensembles
	// set model
	_pModelE = pModelE;

	// generate the sequences
	generateSequences();
#endif

#endif

	// detect the trends
//	_pSequence->TrendDetectRootOnly();
//	_pSequence->TrendDetect();
//	_pSequence->TrendDetectDB();
//	_pSequence->TrendDetectDBImproved();

	_pSequence->TrendDetectDB();

//	_pSequence->TrendDetect();
}

void MyChartWidget::mouseDoubleClickEvent(QMouseEvent *event) {
	_nCurrentGroup++;
	if (_pSequence->GetGroupSize())
	{
		_nCurrentGroup %= _pSequence->GetGroupSize();
	}
	else if (_pSequence->GetDBGroupSize()) {
		_nCurrentGroup %= _pSequence->GetDBGroupSize();
	}
	updateGL();
}

void MyChartWidget::generateSequences() {
	int nWidth = _pModelE->GetW();
	int nHeight = _pModelE->GetH();
	int nEns = _pModelE->GetEnsembleLen();
//	nEns = 5;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < nWidth; j++)
		{
			double dbValue = _pModelE->GetData()->GetData(i, nHeight - 1, j);
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;

			seq.push_back(dbValue);
		}
		//		if(i==0|| i == 1 || i == 2 || i == 3 || i == 4 || i==5||i==8)
		_pSequence->AddSequence(seq);
	}


	_pSequence->Init(nWidth, dbMin, dbMax, g_dbEpsilon);
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
			for(int nIndex: g._member)
			{
				dbY += _pSequence->GetValue(nIndex, k);
			}

			glVertex3f(k, dbY / dbSize, 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawSelectedGroup() {
	if (_pSequence->GetGroupSize() == 0) return;
	glColor3f(0, 1, 0);
	glLineWidth(2.0f);
	glColor3f(0, 1, 0);
	Group g = _pSequence->GetGroup(_nCurrentGroup);
//	cout << g._nS << "\t" << g._nE << endl;
	for(int nIndex : g._member){
		glBegin(GL_LINE_STRIP);
		for (size_t k = g._nS; k <= g._nE; k++)
		{
			glVertex3f(k, _pSequence->GetValue(nIndex,k), 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawGridLines() {
	int nLen = _pSequence->GetLength();
	double dbMin = _pSequence->GetMin();
	double dbMax = _pSequence->GetMax();
	double dbEpsilon = _pSequence->GetEpsilon();

	glColor4f(0, 0, 1,.3);
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
//		double dbSize = g._dbR;
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
			for each (int nMemberIndex in  g._member)
			{
				dbY += _pSequence->GetDBValue(g._dbS, nMemberIndex);
			}
//			for (size_t j = 0; j < g._member.size(); j++)
//			{
//				dbY += _pSequence->GetDBValue(g._dbS, g._member[j]);
//			}

			glVertex3f(g._dbS, dbY / dbSize, 0);
		}
		// the segments of full timesteps
		for (size_t k = nS; k <= nE; k++)
		{
			double dbY = 0;
			for each (int nMemberIndex in  g._member)
			{
				dbY += _pSequence->GetValue(nMemberIndex, k);
			}
//			for (size_t j = 0; j < g._member.size(); j++)
//			{
//				dbY += _pSequence->GetValue(g._member[j], k);
//			}
//
			glVertex3f(k, dbY / dbSize, 0);
		}

		// last segment
		if (g._dbE - nE >.001) {
			double dbY = 0;
			for each (int nMemberIndex in  g._member)
			{
				dbY += _pSequence->GetDBValue(g._dbE, nMemberIndex);
			}
	//		for (size_t j = 0; j < g._member.size(); j++)
	//		{
	//			dbY += _pSequence->GetDBValue(g._dbE, g._member[j]);
	//		}

			glVertex3f(g._dbE, dbY / dbSize, 0);
		}
		glEnd();
	}
}

void MyChartWidget::drawSelectedDBGroup() {
	if (_pSequence->GetDBGroupSize() == 0) return;
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
	for each (int nMemberIndex in  g._member)
	{

		glBegin(GL_LINE_STRIP);

		// first segment
		if (bFirstSeg)
		{
			glVertex3f(g._dbS, _pSequence->GetDBValue(g._dbS, nMemberIndex), 0);
		}
		// the segments of full timesteps
		for (size_t k = nS; k <= nE; k++)
		{
			glVertex3f(k, _pSequence->GetValue(nMemberIndex, k), 0);
		}
		// last segment
		if (g._dbE - nE >.001) {
			glVertex3f(g._dbE, _pSequence->GetDBValue(g._dbE, nMemberIndex), 0);
		}

		glEnd();
	}
}

void MyChartWidget::generateSequences(std::vector<std::vector<double>> vecData) {
	int nLen = vecData[0].size();
	int nEns = vecData.size();

	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (int i = nEns-1; i >=0 ; i--)
	{
		int nFirstValidIndex = -1;
		int nValidCount = 0;
		for (size_t j = 0; j < nLen; j++)
		{
			double dbValue = vecData[i][j];
			if (dbValue>0)
			{
				nValidCount++;
				if (nFirstValidIndex < 0) nFirstValidIndex = j;
				if (dbValue > dbMax) dbMax = dbValue;
				if (dbValue < dbMin) dbMin = dbValue;
			}
			else if(nFirstValidIndex>-1)
				vecData[i][j] = vecData[i][j - 1];
		}
		if (nValidCount < 45)
			vecData.erase(vecData.begin() + i); // remove the sequence with no valide data
		else {
			for (size_t j = 0; j < nFirstValidIndex; j++)
			{
				vecData[i][j] = vecData[i][nFirstValidIndex];
			}
		}
	}
	nEns = vecData.size();
	cout << "number of sequences: " << nEns << endl;

	for each (vector<double> seq in vecData)
	{
		_pSequence->AddSequence(seq);
	}



	_pSequence->Init(nLen, dbMin, dbMax, g_dbEpsilon);
}

void MyChartWidget::generateSequenceArtificial1() {
	double arrData[4][8] = {
		{0.1,0.1,0.1,0.1,1.1,1.1,1.1,1.1 },
		{0.2,0.2,2.2,2.2,2.2,2.2,2.2,2.2 },
		{1.3,1.3,2.3,2.3,1.3,1.3,1.3,1.3 },
		{3.4,3.4,2.4,2.4,1.4,1.4,2.4,2.4 }
	};
	int nEns = 4;
	int nLen = 8;
	double dbEpsilon = .4;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < nLen; j++)
		{
			double dbValue = arrData[i][j];
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;

			seq.push_back(dbValue);
		}
		//		if(i==0|| i == 1 || i == 2 || i == 3 || i == 4 || i==5||i==8)
		_pSequence->AddSequence(seq);
	}


	_pSequence->Init(nLen, dbMin, dbMax, dbEpsilon);
}

void MyChartWidget::generateSequenceArtificial2() {
	double arrData[5][8] = {
		{ 0.01,0.01,0.01,0.01,1.01,1.01,1.01,1.01 },
		{ 0.02,0.02,2.02,2.02,2.02,2.02,2.02,2.02 },
		{ 1.03,1.03,2.03,2.03,1.03,1.03,1.03,1.03 },
		{ 3.04,3.04,2.04,2.04,1.04,1.04,2.04,2.04 },
		{0,0,1,1,2,2,3,3}

	};
	int nEns = 5;
	int nLen = 8;
	double dbEpsilon = .1;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < nLen; j++)
		{
			double dbValue = arrData[i][j];
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;

			seq.push_back(dbValue);
		}
		//		if(i==0|| i == 1 || i == 2 || i == 3 || i == 4 || i==5||i==8)
		_pSequence->AddSequence(seq);
	}


	_pSequence->Init(nLen, dbMin, dbMax, dbEpsilon);
}

void MyChartWidget::generateSequenceArtificial3() {
	double arrData[10][20];

	int nEns = 10;
	int nLen = 20;
	int nGradeMin = 0;
	int nGradeMax = 10;

//	nEns = 5;
//	nLen = 5;
//	nGradeMax = 2;

//	nEns = 10;
//	nLen = 10;
//	nGradeMax = 10;






	int nLastGrade = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		for (size_t j = 0; j < nLen; j++)
		{
			double dbDist = rand()*.2 / RAND_MAX;
			int nGradeRandom = rand();
			if (nLastGrade == nGradeMin) {
				if (nGradeRandom % 2) nLastGrade++;
			}
			else if (nLastGrade == nGradeMax) {
				if (nGradeRandom % 2) nLastGrade--;
			}
			else {
				switch (nGradeRandom%3)
				{
				case 0:
					nLastGrade++;
					break;
				case 1:
					nLastGrade--;
					break;
				default:
					break;
				}
			}
			arrData[i][j] = nLastGrade + dbDist;
		}
	}
	double dbEpsilon = .1;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < nLen; j++)
		{
			double dbValue = arrData[i][j];
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;

			seq.push_back(dbValue);
		}
		//		if(i==0|| i == 1 || i == 2 || i == 3 || i == 4 || i==5||i==8)
		_pSequence->AddSequence(seq);
	}


	_pSequence->Init(nLen, dbMin, dbMax, dbEpsilon);
}

void MyChartWidget::drawEvents() {
	glBegin(GL_LINES);
	vector<DBEpsilonEvent> vecEvents = _pSequence->GetEvents();
	for (DBEpsilonEvent evt : vecEvents) {
		if (evt._nType)
		{
			glColor4f(1, 0, 0, .5);
		}
		else {
			glColor4f(0, 1, 0, .5);
		}
		double dbX = evt._dbTime;
		double dbY1 = _pSequence->GetDBValue(evt._dbTime, evt._nIndex1);
		double dbY2 = _pSequence->GetDBValue(evt._dbTime, evt._nIndex2);
		glVertex3f(dbX,dbY1,0);
		glVertex3f(dbX,dbY2,0);
	}
	glEnd();
}