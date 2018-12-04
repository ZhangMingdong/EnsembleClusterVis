#ifndef CITYSCAPE_H
#define CITYSCAPE_H

#include <QWidget>
#include <QMainWindow>

#include "EnsembleIntersections.h"
#include "ContourGenerator.h"
#include "def.h"

class MyMapWidget;
class MyChartWidget;
class MeteModel;
class ArtificialModel;

class DisplayCtrlWidget;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
	MainWindow();
	~MainWindow();

private:	// action
	QAction *viewShowGridLinesAction;			
	QAction *viewShowBackgroundAction;
	QAction *viewShowIntersectionAction;
	QAction *viewShowUnionBAction;
	QAction *viewShowUnionEAction;
	QAction *viewShowGradientEAction;
	QAction *viewShowLineChartAction;
	QAction *viewShowContourLineTruthAction;
	QAction *viewShowContourLineAction;
	QAction *viewShowContourLineSortedAction;
	QAction *viewShowContourLineSortedSDFAction;
	QAction *viewShowContourLineSDFAction;
	QAction *viewShowContourLineMinAction;
	QAction *viewShowContourLineMaxAction;
	QAction *viewShowContourLineMeanAction;
	QAction *viewShowContourLineMedianAction;
	QAction *viewShowContourLineOutlierAction;
	QAction *viewShowClusterBSAction;		// show cluster of Bayesian variance statistics
	QAction *viewShowClusterBVAction;		// show cluster of Bayesian variance
private:		// widget
	MyMapWidget *_view3D;
	MyChartWidget *_viewChart;
	// control widget
	DisplayCtrlWidget * _pDisplayCtrlWidget;


private:		// model
	// used model
	MeteModel* _pModel;
private:
	void createSceneAndView();
	void createActions();
	void createConnections();
	void createDockWidgets();
	void populateMenuAndToolBar(QMenu *menu, QToolBar *toolBar, QList<QAction*> actions);
	void populateMenusAndToolBars();
public slots:
	void onMousePressed(int x, int y);

protected:
	virtual void closeEvent(QCloseEvent *event);
};

#endif
