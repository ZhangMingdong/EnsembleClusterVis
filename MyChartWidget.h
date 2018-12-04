#pragma once


#include "def.h"
#include "EnsembleIntersections.h"
#include "GLFont.h"

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
	MeteModel* _pModelE = NULL;
	Sequence2D* _pSequence = NULL;
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

	// draw all the groups
	void drawGroups();
	// draw the selected groups
	void drawSelectedGroup();
	void drawGridLines();
	void drawSpaghetti();
	void drawEvents();

	// draw all the groups
	void drawDBGroups();
	// draw the selected groups
	void drawSelectedDBGroup();

	// generate sequences for trend detection -- ensembles
	void generateSequences();

	// genereate sequences from the data -- GDPs
	void generateSequences(std::vector<std::vector<double>> vecData);

	// generate artificial data
	void generateSequenceArtificial1();
	void generateSequenceArtificial2();
	void generateSequenceArtificial3();
};

