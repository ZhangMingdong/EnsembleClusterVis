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
	double _dbStepZ = 0;			// step length in Z, changed by .1
	double _dbBaseZ = 0;			// base height in Z, changed by 10
public:
	void SetModelE(MeteModel* pModelE);
protected:
	virtual void mouseDoubleClickEvent(QMouseEvent *event);
protected:
	// paint the content
	virtual void paint();
	virtual void init();
	virtual void wheelEvent(QWheelEvent * event);
private:
	void drawPoints();				// draw the mds of data points
	void drawLines();
	void drawHierarchy();
	void drawClusters();			// draw hulls of each cluster

};

