#pragma once


#include "def.h"
#include "ContourGenerator.h"
#include "EnsembleIntersections.h"
#include "GLFont.h"
#include "UnCertaintyArea.h"
#include "MeteLayer.h"

#include "MyGLWidget.h"

#include<vector>
#include<unordered_map>

class MeteModel;

class Sequence2D;

class MyChartWidget : public MyGLWidget
{
	Q_OBJECT
public:
	MyChartWidget(QWidget *parent = 0);
	~MyChartWidget();
protected:
	MeteModel* _pModelE;
	Sequence2D* _pSequence;
	int _nCurrentGroup = 0;
public:
	void SetModelE(MeteModel* pModelE);
protected:
	virtual void mouseDoubleClickEvent(QMouseEvent *event);
protected:
	// paint the content
	virtual void paint();
	virtual void init();
protected:
	// generate sequences for trend detection
	void generateSequences();

	// draw all the groups
	void drawGroups();

	void drawSelectedGroup();
	void drawGridLines();
	void drawSpaghetti();
};

