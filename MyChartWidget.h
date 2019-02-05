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
public:
	void SetModelE(MeteModel* pModelE);
protected:
	virtual void mouseDoubleClickEvent(QMouseEvent *event);
protected:
	// paint the content
	virtual void paint();
	virtual void init();

};

