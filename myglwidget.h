#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H
//#include <gl/glew.h>
#include <QGLWidget>

#include "def.h"
#include "ContourGenerator.h"
#include "EnsembleIntersections.h"

#include "GLFont.h"
#include "UnCertaintyArea.h"
#include "MeteLayer.h"

class MeteLayer;
class MeteModel;

#define BUFSIZE 512

class GLUtesselator;
class MyGLWidget : public QGLWidget,ILayerCB
{
	Q_OBJECT

public:

	virtual void DrawText(char* pText, double fX, double fY);
public:
	MyGLWidget(QWidget *parent = 0);
	~MyGLWidget();


protected:
	virtual void initializeGL();
	virtual void paintGL();
	virtual void resizeGL(int width, int height);
	virtual void timerEvent(QTimerEvent* event);
	virtual void mousePressEvent(QMouseEvent * event);
	virtual void mouseDoubleClickEvent(QMouseEvent *event);
	virtual void mouseReleaseEvent(QMouseEvent * event);
	virtual void mouseMoveEvent(QMouseEvent *event);
	virtual void wheelEvent(QWheelEvent * event);
private:
	// state of the mouse
	// 0-noButton
	// 1-L Button down
	// 2-R button down and picking
	int m_nMouseState;
	// rank of the magic cube:2-9;
	int m_nRank;
private:		// mouse trackball
	// Trackball radius
	double m_dbRadius;
	// last trackball position
	double m_dbTrackBallPos[3];
private:		// modelview matrix transformation
	// global modelview matrix
	double m_modelViewMatrixGlobal[16];

private:
	const double c_dbPI;
	const double c_dbSideLen;
private:
	// draw one cube
	void updateTrackBallPos(QPoint pt, double* result);
	// rotate an angle
	void rotate(const double* axis, double angle);
	// move screen from pt1 to pt2
	void move(QPoint pt1, QPoint pt2);
	// reset every thing
	void reset();

	// 2017/02/28
private:

	// truth data
	const double* _pDataT;
	const double* _pDataB;
// 	const double* _pDataE;
	// labels of each grid
	const int *_arrLabels;
	// mean of the variance of Bayesian
	const double* _gridDataMeanVarB;
	// variance of the means of Bayesian
	const double* _gridDataVarMeanB;
	// truth contour
	QList<ContourLine> _listContourTruth[g_temperatureLen];
	QList<ContourLine> _listContourMean[g_temperatureLen];
	QList<ContourLine> _listContourMinB[g_temperatureLen];
	QList<ContourLine> _listContourMaxB[g_temperatureLen];
	QPoint _ptLast;
	// viewport size
	int _nWidth;
	int _nHeight;

	// list of the uncertainty area of intersection of B
	QList<UnCertaintyArea*> _listIntersectionAreaB[g_temperatureLen];

	// list of the uncertainty area of union of B
	QList<UnCertaintyArea*> _listUnionAreaB[g_gradient_l];

	QList<ContourLine> _listContourMin_3[g_gradient_l];
	QList<ContourLine> _listContourMax_3[g_gradient_l];


	GLUtesselator *_tobj;                    /**< 网格化对象指针 */
	GLuint _gllist;                            /**< 显示列表索引 */

	// the selected x and y
	int _nSelectedX;
	int _nSelectedY;

	// selected area, used for brushing
	int _nSelectedLeft;
	int _nSelectedRight;
	int _nSelectedTop;
	int _nSelectedBottom;

	

	const EnsembleIntersections* _arrIntersections;


	GLFont font;                /**< 定义一个字体实例 */

	// the map area
	double _fLeft;
	double _fRight;
	double _fBottom;
	double _fTop;
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


public:
	void SetLabels(const int* labels);
	void SetDataT(const double* data);
	void SetDataB(const double* data);
	void SetModelE(MeteModel* pModelE);
	void SetModelT(MeteModel* pModelT);
	void SetVar(const double* pVM, const double *pMV);
	void SetContourTruth(QList<ContourLine>* listContourTruth);
	void SetContourMean(QList<ContourLine>* listContour);
	void SetContourIntervalB(QList<ContourLine>* listContourMin, QList<ContourLine>* listContourMax);

	// generate the contours from the min and max contour lines
	void generateContour();

	void SetMultiStatistic(QList<ContourLine>* listContourMin, QList<ContourLine>* listContourMax);

	void SetIntersections(const EnsembleIntersections* arrIntersections){ _arrIntersections = arrIntersections; }


	// reload texture
	void ReloadTexture();
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
	void viewShowContourLineMin(bool on);
	void viewShowContourLineMax(bool on);
	void viewShowContourLineMean(bool on);
	void viewShowClusterBS(bool on);
	void viewShowClusterBV(bool on);

	void onCheckShowBeliefEllipse(bool bChecked);

private:
	// vector of layers to render
	std::vector<MeteLayer*> _vecLayers;
	MeteModel* _pModelE;
	MeteModel* _pModelT;


	

};

#endif // MYGLWIDGET_H
