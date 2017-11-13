#include "mychartwidget.h"
#include <gl/GLU.h>

#include <QDebug>


#include <iostream>
#include <algorithm>


#include "MeteLayer.h"
#include "CostLineLayer.h"
#include "EnsembleLayer.h"
#include "ClusterLayer.h"
#include "MeteModel.h"
#include "VarAnalysis.h"
#include "LayerLayout.h"
#include "DataField.h"
#include "ColorMap.h"

#define BUFSIZE 512

#define REVERSE

using namespace std;

MyChartWidget::MyChartWidget(QWidget *parent)
	: MyGLWidget(parent)
{
	startTimer(100);
}

MyChartWidget::~MyChartWidget()
{

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


	glScaled(.05, .03, 1);
	glTranslatef(-_nWidth / 2.0, -(_dbMin + _dbMax) / 2.0, 0);

	// draw spaghetti
	for (size_t i = 0; i < _nEns; i++)
	{
		glBegin(GL_LINE_STRIP);
		for (size_t j = 0; j < _nWidth; j++)
		{
			glVertex3f(j, _vecSequences[i][j], 0);
		}
		glEnd();
	}
	// draw grid lines
	glColor3f(0, 0, .5);	
	glBegin(GL_LINES);
	for (size_t i = 0; i <= _nEns; i++)
	{
		glVertex2d(0, _dbMin+i*_dbEpsilon);
		glVertex2d(_nWidth-1, _dbMin + i*_dbEpsilon);
	}
	for (size_t i = 0; i < _nWidth; i+=10)
	{
		glVertex2d(i, _dbMin);
		glVertex2d(i, _dbMax);
	}
	glEnd();

	/*
	glColor3f(0, 0, .5);
	glBegin(GL_POLYGON);
	glVertex2d(0, _dbMin);
	glVertex2d(_nWidth, _dbMin);
	glVertex2d(_nWidth, _dbMax);
	glVertex2d(0, _dbMax);
	glEnd();
	*/

	// draw selected group
	/*
	glColor3f(0, 1, 0);
	glLineWidth(2.0f);
	glColor3f(0, 1, 0);
	for (size_t j = 0; j < _vecGroups[_nCurrentGroup]._member.size(); j++)
	{
		glBegin(GL_LINE_STRIP);
		for (size_t k = _vecGroups[_nCurrentGroup]._nS; k <= _vecGroups[_nCurrentGroup]._nE; k++)
		{
			glVertex3f(k, _vecSequences[_vecGroups[_nCurrentGroup]._member[j]][k], 0);
		}
		glEnd();
	}
	*/
	// draww abstracted groups
	drawGroups();

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
	trendDetect();
//	trendDetectRootOnly();
}

void MyChartWidget::generateSequences() {
	_nWidth = _pModelE->GetW();
	_nHeight = _pModelE->GetH();
	_nEns = _pModelE->GetEnsembleLen();
	_nEns = 12;
	// calculate the max and min height
	double dbMin = 1000;
	double dbMax = 0;
	for (size_t i = 0; i < _nEns; i++)
	{
		vector<double> seq;
		for (size_t j = 0; j < _nWidth; j++)
		{
#ifdef REVERSE
			double dbValue = _pModelE->GetData()->GetData(i, _nHeight - 1, _nWidth - 1 - j);
#else
			double dbValue = _pModelE->GetData()->GetData(i, _nHeight - 1, j);
#endif
			if (dbValue > dbMax) dbMax = dbValue;
			if (dbValue < dbMin) dbMin = dbValue;
			seq.push_back(dbValue);
		}
		_vecSequences.push_back(seq);
	}
	_dbMin = dbMin;
	_dbMax = dbMax;
	_dbEpsilon = (dbMax - dbMin) *1.2 / _nEns;
}
void MyChartWidget::trendDetectRootOnly() {
	// candidates of groups
	Group group0;
	group0._nS = 0;
	group0._nE = 0;
	for (size_t i = 0; i < _nEns; i++)
	{
		group0._member.push_back(i);
	}

	// start from the first step
	trendDetectStep(group0, 0);
}

void MyChartWidget::trendDetect() {
	for (size_t i = 0; i < _nWidth-_nDelta; i++)
	{
		Group group0;
		group0._nS = i;
		for (size_t i = 0; i < _nEns; i++)
		{
			group0._member.push_back(i);
		}

		// start from the first step
		trendDetectStep(group0, i);
	}
	// remove duplicated
	for (int i = _vecGroups.size() - 1; i >= 0; i--)
	{
		for (int j = 0, length = _vecGroups.size(); j < length; j++) {
			if (i!=j&& _vecGroups[j].Contain(_vecGroups[i]))
			{
				_vecGroups.erase(_vecGroups.begin() + i);
				break;
			}
		}
	}

	for each (Group g in _vecGroups)
	{
#ifdef REVERSE
		cout << _nWidth - 1 - g._nE << "\t" << _nWidth - 1 - g._nS << "\t";
#else
		cout << g._nS << "\t" << g._nE << "\t";
#endif		
		for each (int nMember in g._member)
		{
			cout << nMember << "\t";
		}
		cout << endl;
	}
}

void MyChartWidget::trendDetect_1() {
	for (size_t i = 0; i < _nWidth - _nDelta; i++)
	{
		Group group0;
		group0._nS = i;
		for (size_t i = 0; i < _nEns; i++)
		{
			group0._member.push_back(i);
		}

		// start from the first step
		trendDetectStep(group0, i);
	}
	// remove duplicated
	for (int i = _vecGroups.size() - 1; i >= 0; i--)
	{
		for (int j = 0, length = _vecGroups.size(); j < length; j++) {
			if (i != j&& _vecGroups[j].Contain(_vecGroups[i]))
			{
				_vecGroups.erase(_vecGroups.begin() + i);
				break;
			}
		}
	}
	for each (Group g in _vecGroups)
	{
#ifdef REVERSE
		cout << _nWidth - 1 - g._nE << "\t" << _nWidth - 1 - g._nS << "\t";
#else
		cout << g._nS << "\t" << g._nE << "\t";
#endif		
		for each (int nMember in g._member)
		{
			cout << nMember << "\t";
		}
		cout << endl;
	}
}
// cluster comparison function, used for sort
bool IndexAndValueCompare(IndexAndValue iv1, IndexAndValue iv2) {
	return iv1._dbValue < iv2._dbValue;
}

void MyChartWidget::trendDetectStep(Group candidate, int nStep) {
	// 0.finished if nStep is the end
	if (nStep==_nWidth)
	{
		candidate._nE = nStep - 1;
		addGroup(candidate);
		return;
	}

	// 1.get the point in the new timestep
	vector<IndexAndValue> vecNewPoints;
	for (size_t i = 0,length=candidate._member.size(); i < length; i++)
	{
		vecNewPoints.push_back(IndexAndValue(candidate._member[i], _vecSequences[candidate._member[i]][nStep]));
	}

	// 2.sort them according to the value
	sort(vecNewPoints.begin(), vecNewPoints.end(), IndexAndValueCompare);

	// 3.split the points in this timestep
	vector<vector<IndexAndValue>> vecNewPointArray;
	int nCurrentStartPos = 0;	// start position of current group;
	for (size_t i = 0; i < vecNewPoints.size() - 1; i++)
	{
		if (vecNewPoints[i + 1]._dbValue - vecNewPoints[i]._dbValue>_dbEpsilon) {
			// space larger than epsilon, split
			vector<IndexAndValue> newArray;
			for (size_t j = nCurrentStartPos; j <= i; j++)
			{
				newArray.push_back(vecNewPoints[j]);
			}
			vecNewPointArray.push_back(newArray);
			nCurrentStartPos = i + 1;
		}
	}
	// add the last segment
	vector<IndexAndValue> newArray;
	for (size_t j = nCurrentStartPos; j <vecNewPoints.size(); j++)
	{
		newArray.push_back(vecNewPoints[j]);
	}
	vecNewPointArray.push_back(newArray);

	// 4.handle the splited points
	if (vecNewPointArray.size()==1)
	{
		// no splitting
		trendDetectStep(candidate, nStep + 1);
	}
	else {
		// generate a new group of the candidate
		candidate._nE = nStep - 1;
		addGroup(candidate);
		// generate new candidates
		for (size_t i = 0,length=vecNewPointArray.size(); i < length; i++)
		{
			int length2 = vecNewPointArray[i].size();
			if (length2 > _nM) {
				Group newCandidate;
				newCandidate._nS = candidate._nS;
				for (size_t j = 0; j < length2; j++)
				{
					newCandidate._member.push_back(vecNewPointArray[i][j]._nIndex);
				}
				// reset the start time of the new candidate
				newCandidate._nS = calculateStart(newCandidate, nStep);
				trendDetectStep(newCandidate, nStep + 1);
			}
		}
	}
}

void MyChartWidget::mouseDoubleClickEvent(QMouseEvent *event) {
	_nCurrentGroup++;
	_nCurrentGroup %= _vecGroups.size();
	updateGL();
}

void MyChartWidget::drawGroups() {
	for (size_t i = 0; i < _vecGroups.size(); i++)
	{
		glColor4f(0, 1, 1,.5);
		double dbSize = _vecGroups[i]._member.size();
		glLineWidth(dbSize);

		glBegin(GL_LINE_STRIP);
		for (size_t k = _vecGroups[i]._nS; k <= _vecGroups[i]._nE; k++)
		{
			double dbY = 0;
			for (size_t j = 0; j < _vecGroups[i]._member.size(); j++)
			{
				dbY += _vecSequences[_vecGroups[i]._member[j]][k];
			}

			glVertex3f(k, dbY/dbSize, 0);
		}
		glEnd();


	}
}

void MyChartWidget::addGroup(Group candidate) {
	if (candidate._nE - candidate._nS>_nDelta && candidate._member.size()>_nM)
	{
		sort(candidate._member.begin(), candidate._member.end());
		_vecGroups.push_back(candidate);
	}
}

int MyChartWidget::calculateStart(Group g, int nCurrent) {
	int nS = nCurrent - 1;
	while (nS>=g._nS)
	{
		// 1.get the point in the new timestep
		vector<IndexAndValue> vecNewPoints;
		for (size_t i = 0, length = g._member.size(); i < length; i++)
		{
			vecNewPoints.push_back(IndexAndValue(g._member[i], _vecSequences[g._member[i]][nS]));
		}

		// 2.sort them according to the value
		sort(vecNewPoints.begin(), vecNewPoints.end(), IndexAndValueCompare);


		int nCurrentStartPos = 0;	// start position of current group;
		for (size_t i = 0; i < vecNewPoints.size() - 1; i++)
		{
			if (vecNewPoints[i + 1]._dbValue - vecNewPoints[i]._dbValue>_dbEpsilon)
				return nS+1;
		}

		nS--;
	}

	return g._nS;
}