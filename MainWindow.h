#ifndef CITYSCAPE_H
#define CITYSCAPE_H

#include <QWidget>
#include <QMainWindow>

#include "EnsembleIntersections.h"
#include "ContourGenerator.h"
#include "def.h"

#include "ControlWidget.h"

class QGraphicsScene;
class MainView;
class TableWidget;
class Widget3D;
class MyGLWidget;

class MeteModel;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
	MainWindow();
	~MainWindow();

private:
	MyGLWidget *_view3D;
	// table
	ControlWidget *_pControlWidget;

	// view
	QAction *viewShowGridLinesAction;			
	QAction *viewShowBackgroundAction;
	QAction *viewShowIntersectionAction;
	QAction *viewShowUnionBAction;
	QAction *viewShowUnionEAction;
	QAction *viewShowGradientEAction;
	QAction *viewShowLineChartAction;
	QAction *viewShowContourLineTruthAction;
	QAction *viewShowContourLineAction;
	QAction *viewShowContourLineMinAction;
	QAction *viewShowContourLineMaxAction;
	QAction *viewShowContourLineMeanAction;
	QAction *viewShowClusterBSAction;		// show cluster of Bayesian variance statistics
	QAction *viewShowClusterBVAction;		// show cluster of Bayesian variance
private:
	void createSceneAndView();
	void createActions();
	void createConnections();
	void createDockWidgets();
	void populateMenuAndToolBar(QMenu *menu, QToolBar *toolBar, QList<QAction*> actions);
	void populateMenusAndToolBars();
	// data reading..............
	// read the data of Bayesian Result
	void readDataB();
	// read the data of truth
	void readDataT();

	// calculate the max, min, mean, and variance
	void doStatisticB();
	void doStatistic3B();


	// directly segment, no clustering
	void doSegment();
	// build the distance matrix for clustering using ensemble variance
	void buildDistanceMatrix_ensemble(double** disMatrix);
	// build the distance matrix for clustering using variance and bias
	void buildDistanceMatrix_vb(double** disMatrix);

	// read the labels of clustering
	void readLabels();

	// ==================Raw Data=================

	// the Bayesian result: [mb,vb,m,v]
	double* _dataB;

	// the data of truth
	double* _dataT;
	// ==================~Raw Data=================

	// ==================Statistic Data=================


	// max and min of the union of Bayesian
	double* _gridDataUnionMaxB;
	double* _gridDataUnionMinB;
	// mean of mean of the Bayesian
	double* _gridDataMeanMeanB;
	// mean of the variance of Bayesian
	double* _gridDataMeanVarB;
	// variance of the means of Bayesian
	double* _gridDataVarMeanB;
	// normalized vm and mv
	double* _normalizedVMB;
	double* _normalizedMVB;

	double* _gridIntersectionMaxB;
	double* _gridIntersectionMinB;


	// the intersections of each grid points
	EnsembleIntersections* _arrIntersections;
	// ==================~Statistic Data=================



	// ==================Grid Data=================
	// min,mean,max of every grid for each result
	double* _gridDataMinB	;
	double* _gridDataMeanB	;
	double* _gridDataMaxB	;

	// for gradient strip
	double* _gridDataMinB_3[g_gradient_l];
	double* _gridDataMaxB_3[g_gradient_l];
	double* _gridDataMaxGradientB[g_gradient_l];
	double* _gridDataMinGradientB[g_gradient_l];

	// label of clustering
	int* _gridLabels;
	double* _gridLabelsF;

	// mean Bayesian vars of according to latitude
	double* _arrLatVarB;



	// ==================~Grid Data=================

	// contours of the truth data
	QList<ContourLine> _listContourT[g_temperatureLen];
	// contours of the ensemble
	QList<ContourLine> _listContourE[g_temperatureLen];
	// contours of the ensemble mean
	QList<ContourLine> _listContourEnsembleMeanB[g_temperatureLen];

	QList<ContourLine> _listContourEnsembleMin_3[g_gradient_l];
	QList<ContourLine> _listContourEnsembleMax_3[g_gradient_l];


	// the most overlapped intersection min and max
	QList<ContourLine> _listContourIntersectionMinB[g_temperatureLen];
	QList<ContourLine> _listContourIntersectionMaxB[g_temperatureLen];

	



	// used model
	MeteModel* _pModel;
	


public slots:
	void onMousePressed(int x, int y);

	void onSelectBackgroundMean(bool bChecked);
	void onSelectBackgroundVari(bool bChecked);
	void onSelectBackgroundCluster(bool bChecked);
	void onSelectBackgroundSDF(bool bChecked);
protected:
	virtual void closeEvent(QCloseEvent *event);
};

#endif
