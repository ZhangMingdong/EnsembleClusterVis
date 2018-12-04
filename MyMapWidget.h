#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H


#include "def.h"
#include "ContourGenerator.h"
#include "EnsembleIntersections.h"
#include "GLFont.h"
#include "UnCertaintyArea.h"
#include "MeteLayer.h"

#include "MyGLWidget.h"

class MeteLayer;
class MeteModel;
class LayerLayout;


class MyMapWidget : public MyGLWidget,ILayerCB
{
	Q_OBJECT

public:
	// draw text
	virtual void DrawText(char* pText, double fX, double fY);
public:
	MyMapWidget(QWidget *parent = 0);
	~MyMapWidget();
protected:
//	virtual void resizeGL(int width, int height);
	virtual void timerEvent(QTimerEvent* event);
	virtual void mousePressEvent(QMouseEvent * event);
	virtual void mouseDoubleClickEvent(QMouseEvent *event);
	virtual void mouseReleaseEvent(QMouseEvent * event);
	virtual void mouseMoveEvent(QMouseEvent *event);
//	virtual void wheelEvent(QWheelEvent * event);
private:

	GLFont _font;                /**< 定义一个字体实例 */

	// layout of the layers
	LayerLayout* _pLayout = NULL;

	// scale: 0.01 (180/1.8)
	double _fScaleW;
	double _fScaleH;

	double _fChartW;
	double _fChartLeft;
	double _fChartRight;

	double _fChartLW;
	double _fChartLLeft;
	double _fChartLRight;

	struct PickIndex
	{
		int _nX;
		int _nY;
		PickIndex(){}
		PickIndex(int x, int y) :_nX(x), _nY(y){}
	};

private:
	void drawPlaceHolder();
	// process the hits and trigger rotateLayer
	bool processHits(GLint hits, GLuint buffer[], PickIndex& pick);

	// draw a contour line
	void drawContourLine(const QList<ContourLine>& contours);

	// generate the background texture
	void generateBackground();

	// select the grid point when clicked at pt
	void select(int& nX, int&nY, const QPoint& pt);


	// called when texture reloaded
	void onTextureReloaded();
	// set new uncertainty areas
	void setUncertaintyAreas(int nUCAreas);
public:
	void SetModelE(MeteModel* pModelE);
private:// state
	MeteLayer::DisplayStates _displayStates;
	bool _bShowGridLines;
	bool _bShowIntersection;
	bool _bShowUnionB;
	bool _bShowLineChart;
	bool _bShowContourLineTruth;
	bool _bShowClusterBS;
	bool _bShowClusterBV;
public slots:
	void viewShowGrid(bool on);
	void viewShowBackground(bool on);
	void viewShowIntersection(bool on);
	void viewShowUnionB(bool on);
	void viewShowUnionE(bool on);
	void viewShowGradientE(bool on);
	void viewShowLineChart(bool on);
	void viewShowContourLineTruth(bool on);
	void viewShowContourLine(bool on);
	void viewShowContourLineSorted(bool on);
	void viewShowContourLineSortedSDF(bool on);
	void viewShowContourLineSDF(bool on);
	void viewShowContourLineMin(bool on);
	void viewShowContourLineMax(bool on);
	void viewShowContourLineMean(bool on);
	void viewShowContourLineMedian(bool on);
	void viewShowContourLineOutlier(bool on);
	void viewShowClusterBS(bool on);
	void viewShowClusterBV(bool on);
	void onCheckShowBeliefEllipse(bool bChecked);
	// for clustering
	void updateMinPts(int minPts);
	void updateEps(double eps);
	// for clustered variance
	void updateVarSmooth(int nSmooth);
	void updateBgFunction(int nBgFunction);
	void updateUncertaintyAreas(int nAreas);
	void updateFocusedCluster(int nFocusedCluster);
	void updateFocusedRegion(int nFocusedRegion);
	void updateEOF(int nEOF);
	void updateMember(int nMember);
	void updateEnsCluster(int nEnsCluster);
	void updateContourLevel(int nLevel);
private:
	// vector of layers to render
	std::vector<MeteLayer*> _vecLayers;
	MeteModel* _pModelE;



protected:
	// paint the content
	virtual void paint(); 
	virtual void init();

};

#endif // MYGLWIDGET_H
