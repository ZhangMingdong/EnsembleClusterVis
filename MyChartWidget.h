#pragma once


#include "def.h"
#include "ContourGenerator.h"
#include "EnsembleIntersections.h"
#include "GLFont.h"
#include "UnCertaintyArea.h"
#include "MeteLayer.h"

#include "MyGLWidget.h"

#include<vector>

class MeteModel;

struct Group {
	std::vector<int> _member;	// members
	int _nS;					// start time
	int _nE;					// end time
	// check whether containing g
	bool Contain(const Group& g) {
		int nLen = _member.size();
		if (g._member.size()==nLen)
		{
			for (size_t i = 0; i < nLen; i++)
				if (g._member[i] != _member[i]) return false;
			return _nS <= g._nS&&_nE >= g._nE;
		}
		return false;
	}
};
// pair of index and value
struct IndexAndValue {
	int _nIndex = 0;
	double _dbValue = 0;
	IndexAndValue(int nIndex, double dbValue) :_nIndex(nIndex), _dbValue(dbValue) {};
};


class MyChartWidget : public MyGLWidget
{
	Q_OBJECT
public:
	MyChartWidget(QWidget *parent = 0);
	~MyChartWidget();
protected:
	MeteModel* _pModelE;
	// the sequences, same length
	std::vector<std::vector<double>> _vecSequences;
	int _nWidth = 0;
	int _nHeight = 0;
	int _nEns = 0;
	double _dbMin = 0;
	double _dbMax = 0;
	// groups
	std::vector<Group> _vecGroups;
	// trend detecting parameters
	int _nM = 2;
	double _dbEpsilon = 0;
	int _nDelta = 10;
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
	// detect the trend from the root
	void trendDetectRootOnly();
	// detect the trend
	void trendDetect();
	// detect the trend
	// old version, find all and then remove duplicated
	void trendDetect_1();
	// detect the trend at a single step
	// recursive version, only find the trend start from beginning
	void trendDetectStep(Group candidate, int nStep);
private:
	// draw all the groups
	void drawGroups();

	// add a new group
	void addGroup(Group g);

	// calculate the start time of g
	// which whould between g.start and nCurrent
	int calculateStart(Group g, int nCurrent);
};

