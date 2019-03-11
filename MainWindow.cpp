#include <QtGui>
#include <QtWidgets/QGridLayout>
#include <QLabel>
#include <QMenu>
#include <QToolBar>
#include <QMenuBar>
#include <QSet>
#include <QMessageBox>
#include <QDockWidget>
#include <QSplitter>




#include "MainWindow.h"
#include "MyMapWidget.h"
#include "MyChartWidget.h"
#include "MeteModel.h"
#include "ArtificialModel.h"

#include "DisplayCtrlWidget.h"

#include "ClusteringWidget.h"

#include "Switch.h"



const QString ShowGridLines("ShowGridLines");
const QString ShowBackground("ShowBackground");
const QString ShowIntersection("ShowIntersection");
const QString ShowUnionB("ShowUnionB");
const QString ShowUnionE("ShowUnionE");
const QString ShowGradientE("ShowGradientE");
const QString ShowLineChart("ShowLineChart");
const QString ShowContourLineTruth("ShowContourLineTruth");
const QString ShowContourLine("ShowContourLine");
const QString ShowContourLineSorted("ShowContourLineSorted");
const QString ShowContourLineSortedSDF("ShowContourLineSortedSDF");
const QString ShowContourLineResampled("ShowContourLineResampled");
const QString ShowContourLineSDF("ShowContourLineSDF");
const QString ShowContourLineMin("ShowContourLineMin");
const QString ShowContourLineMax("ShowContourLineMax");
const QString ShowContourLineMean("ShowContourLineMean");
const QString ShowContourLineMedian("ShowContourLineMedian");
const QString ShowContourLineOutlier("ShowContourLineOutlier");
const QString ShowClusterBS("ShowClusterBS");
const QString ShowClusterBV("ShowClusterBV");
const QString ShowBeliefEllipse("ShowBeliefEllipse");


MainWindow::~MainWindow(){
	if (_pModel) delete _pModel;
}


MainWindow::MainWindow()
{
	// finished read config file
	setWindowState(Qt::WindowMaximized);

	_pModel = MeteModel::CreateModel(g_bArtificialModel);

	createSceneAndView();
	createActions();

	createDockWidgets();
	createConnections();
	populateMenusAndToolBars();

	QSettings settings;
	viewShowGridLinesAction->setChecked(settings.value(ShowGridLines, true).toBool());
	viewShowBackgroundAction->setChecked(settings.value(ShowBackground, true).toBool());
	viewShowIntersectionAction->setChecked(settings.value(ShowIntersection, true).toBool());
	viewShowUnionBAction->setChecked(settings.value(ShowUnionB, true).toBool());
	viewShowUnionEAction->setChecked(settings.value(ShowUnionE, true).toBool());
	viewShowGradientEAction->setChecked(settings.value(ShowGradientE, true).toBool());
	viewShowLineChartAction->setChecked(settings.value(ShowLineChart, true).toBool());
	viewShowContourLineTruthAction->setChecked(settings.value(ShowContourLineTruth, true).toBool());
	viewShowContourLineAction->setChecked(settings.value(ShowContourLine, true).toBool());
	viewShowContourLineSortedAction->setChecked(settings.value(ShowContourLineSorted, true).toBool());
	viewShowContourLineSortedSDFAction->setChecked(settings.value(ShowContourLineSortedSDF, true).toBool());
	viewShowContourLineResampledAction->setChecked(settings.value(ShowContourLineResampled, true).toBool());
	viewShowContourLineSDFAction->setChecked(settings.value(ShowContourLineSDF, true).toBool());
	viewShowContourLineMinAction->setChecked(settings.value(ShowContourLineMin, true).toBool());
	viewShowContourLineMaxAction->setChecked(settings.value(ShowContourLineMax, true).toBool());
	viewShowContourLineMeanAction->setChecked(settings.value(ShowContourLineMean, true).toBool());
	viewShowContourLineMedianAction->setChecked(settings.value(ShowContourLineMedian, true).toBool());
	viewShowContourLineOutlierAction->setChecked(settings.value(ShowContourLineOutlier, true).toBool());
	viewShowClusterBSAction->setChecked(settings.value(ShowClusterBS, true).toBool());
	viewShowClusterBVAction->setChecked(settings.value(ShowClusterBV, true).toBool());

//	_pControlWidget->ui.radioButtonBackgroundCluster->setChecked(settings.value(ShowClusterBV, true).toBool());


    setWindowTitle(tr("Ensembles"));
}

void MainWindow::createSceneAndView(){
	_view3D = new MyMapWidget;
	_view3D->SetModelE(_pModel);

	if (g_bChart) {
		QSplitter* splitter = new QSplitter();
		splitter->setOrientation(Qt::Vertical);
		splitter->addWidget(_view3D);

		_viewChart = new MyChartWidget;
		_viewChart->SetModelE(_pModel);
		splitter->addWidget(_viewChart);

		setCentralWidget(splitter);

	}
	else {
		setCentralWidget(_view3D);
	}


	//
}

void MainWindow::createDockWidgets() {
	setDockOptions(QMainWindow::AnimatedDocks);
	QDockWidget::DockWidgetFeatures features =
		QDockWidget::DockWidgetMovable |
		QDockWidget::DockWidgetFloatable;

	_pDisplayCtrlWidget = new DisplayCtrlWidget();


	QDockWidget *controlDockWidget = new QDockWidget(
		tr("Control"), this);
	controlDockWidget->setFeatures(features);
	controlDockWidget->setWidget(_pDisplayCtrlWidget);
	addDockWidget(Qt::RightDockWidgetArea, controlDockWidget);


	/*
	ClusteringWidget* _pWidgetClustering = new ClusteringWidget();
	QDockWidget *pDockWidgetClustering = new QDockWidget(tr("Clustering"), this);
	pDockWidgetClustering->setFeatures(features);
	pDockWidgetClustering->setWidget(_pWidgetClustering);
	addDockWidget(Qt::LeftDockWidgetArea, pDockWidgetClustering);
	*/
}

void MainWindow::populateMenusAndToolBars()
{
	QAction *separator = 0;
	// file
// 	QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
// 	QToolBar *fileToolBar = addToolBar(tr("File"));
// 	populateMenuAndToolBar(fileMenu, fileToolBar, QList<QAction*>()
// 		<< fileNewAction
// 		<< fileOpenAction
// 		<< fileSaveAction
// 		// 		<< fileExportAction 
// 		<< exportSVGAction
// 		<< separator
// 		<< setBackgroundAction
// 		<< separator
// 		<< filePrintAction);
// 	fileMenu->insertAction(exportSVGAction, fileSaveAsAction);
// 	fileMenu->addSeparator();
// 	fileMenu->addAction(fileQuitAction);
// 
// 	// edit
// 	editMenu = menuBar()->addMenu(tr("&Edit"));
// 	QToolBar* editToolBar = addToolBar(tr("Edit"));
// 	populateMenuAndToolBar(editMenu, editToolBar, QList<QAction*>()
// 		<< addSphereAction
// 		<< addRectAction
// 		<< addLineAction
// 		<< addTextAction
// 		<< addArcAction
// 		<< addSwitchAction
// 		<< addLinkAction
// 		<< editDeleteAction
// 		<< separator
// 		<< editCutAction
// 		<< editCopyAction
// 		<< editPasteAction
// 		<< separator
// 		<< bringToFrontAction
// 		<< sendToBackAction
// 		<< separator
// 		<< editAlignmentAction
// 		<< editClearTransformsAction
// 		<< separator
// 		<< editCombineAction
// 		<< editDetachAction
// 		<< separator
// 		<< propertiesAction
// 		<< editUndoAction
// 		<< editRedoAction);


	// view
	QMenu* viewMenu = menuBar()->addMenu(tr("&View"));
	QToolBar* viewToolBar = addToolBar(tr("View"));
	populateMenuAndToolBar(viewMenu, viewToolBar, QList<QAction*>()
// 		<< viewZoomInAction
// 		<< viewZoomOutAction
		<< separator
		<< viewShowGridLinesAction
		<< viewShowBackgroundAction
//		<< viewShowIntersectionAction
//		<< viewShowUnionBAction
		<< viewShowUnionEAction
		<< viewShowGradientEAction
//		<< viewShowLineChartAction
		<< viewShowContourLineTruthAction
		<< viewShowContourLineAction
		<< viewShowContourLineSortedAction
		<< viewShowContourLineSDFAction
		<< viewShowContourLineSortedSDFAction
		<< viewShowContourLineResampledAction
		<< viewShowContourLineMinAction
		<< viewShowContourLineMaxAction
		<< viewShowContourLineMeanAction
		<< viewShowContourLineMedianAction
		<< viewShowContourLineOutlierAction
//		<< viewShowClusterBSAction
//		<< viewShowClusterBVAction
);

// 	// state
// 	QMenu* stateMenu = menuBar()->addMenu(tr("State"));
// 	QToolBar* stateToolBar = addToolBar(tr("State"));
// 	populateMenuAndToolBar(stateMenu, stateToolBar, QList<QAction*>()
// 		<< drawSphereAction
// 		<< drawRectAction
// 		<< drawLineAction
// 		<< drawTextAction
// 		<< drawLinkAction
// 		<< stateSelectAction);
// 
// 	// layer
// 	QMenu* layerMenu = menuBar()->addMenu(tr("&Layer"));
// 	layerMenu->addAction(newLayerAction);

}

void MainWindow::populateMenuAndToolBar(QMenu *menu, QToolBar *toolBar, QList<QAction*> actions)
{
	foreach(QAction *action, actions) {
		if (!action) {
			menu->addSeparator();
			toolBar->addSeparator();
		}
		else {
			menu->addAction(action);
			toolBar->addAction(action);
		}
	}
}

void MainWindow::createActions()
{
	// view
	viewShowGridLinesAction = new QAction(tr("Show Gridlines"), this);
	viewShowGridLinesAction->setIcon(QIcon(":/images/showgrid.png"));
	viewShowGridLinesAction->setCheckable(true);

	viewShowBackgroundAction = new QAction(tr("Show Background"), this);
	viewShowBackgroundAction->setIcon(QIcon(":/images/showgrid.png"));
	viewShowBackgroundAction->setCheckable(true);

	viewShowIntersectionAction = new QAction(tr("Show Intersection"), this);
	viewShowIntersectionAction->setIcon(QIcon(":/images/ib.png"));
	viewShowIntersectionAction->setCheckable(true);

	viewShowUnionBAction = new QAction(tr("Show Union B"), this);
	viewShowUnionBAction->setIcon(QIcon(":/images/ub.png"));
	viewShowUnionBAction->setCheckable(true);

	viewShowUnionEAction = new QAction(tr("Show Union E"), this);
	viewShowUnionEAction->setIcon(QIcon(":/images/ue.png"));
	viewShowUnionEAction->setCheckable(true);

	viewShowGradientEAction = new QAction(tr("Show Gradient E"), this);
	viewShowGradientEAction->setIcon(QIcon(":/images/ge.png"));
	viewShowGradientEAction->setCheckable(true);

	viewShowLineChartAction = new QAction(tr("Show LineChart"), this);
	viewShowLineChartAction->setIcon(QIcon(":/images/showgrid.png"));
	viewShowLineChartAction->setCheckable(true);

	viewShowContourLineTruthAction = new QAction(tr("Show Contour of T"), this);
	viewShowContourLineTruthAction->setIcon(QIcon(":/images/t.png"));
	viewShowContourLineTruthAction->setCheckable(true);

	viewShowContourLineAction = new QAction(tr("Show Contour of E"), this);
	viewShowContourLineAction->setIcon(QIcon(":/images/e.png"));
	viewShowContourLineAction->setCheckable(true);

	viewShowContourLineSDFAction = new QAction(tr("Show Contour of SDF"), this);
	viewShowContourLineSDFAction->setIcon(QIcon(":/images/e.png"));
	viewShowContourLineSDFAction->setCheckable(true);


	viewShowContourLineSortedAction = new QAction(tr("Show Sorted Contour of E"), this);
	viewShowContourLineSortedAction->setIcon(QIcon(":/images/e.png"));
	viewShowContourLineSortedAction->setCheckable(true);


	viewShowContourLineSortedSDFAction = new QAction(tr("Show Sorted Contour of SDF"), this);
	viewShowContourLineSortedSDFAction->setIcon(QIcon(":/images/e.png"));
	viewShowContourLineSortedSDFAction->setCheckable(true);

	viewShowContourLineResampledAction = new QAction(tr("Show Resampled Contours"), this);
	viewShowContourLineResampledAction->setIcon(QIcon(":/images/e.png"));
	viewShowContourLineResampledAction->setCheckable(true);

	viewShowContourLineMinAction = new QAction(tr("Show Contour of Min"), this);
	viewShowContourLineMinAction->setIcon(QIcon(":/images/min.png"));
	viewShowContourLineMinAction->setCheckable(true);

	viewShowContourLineMaxAction = new QAction(tr("Show Contour of Max"), this);
	viewShowContourLineMaxAction->setIcon(QIcon(":/images/max.png"));
	viewShowContourLineMaxAction->setCheckable(true);

	viewShowContourLineMeanAction = new QAction(tr("Show Contour of Mean"), this);
	viewShowContourLineMeanAction->setIcon(QIcon(":/images/b.png"));
	viewShowContourLineMeanAction->setCheckable(true);

	viewShowContourLineMedianAction = new QAction(tr("Show Contour of Median"), this);
	viewShowContourLineMedianAction->setIcon(QIcon(":/images/b.png"));
	viewShowContourLineMedianAction->setCheckable(true);

	viewShowContourLineOutlierAction = new QAction(tr("Show Contour of Outlier"), this);
	viewShowContourLineOutlierAction->setIcon(QIcon(":/images/b.png"));
	viewShowContourLineOutlierAction->setCheckable(true);


	viewShowClusterBSAction = new QAction(tr("Show Cluster BS"), this);
	viewShowClusterBSAction->setIcon(QIcon(":/images/showgrid.png"));
	viewShowClusterBSAction->setCheckable(true);

	viewShowClusterBVAction = new QAction(tr("Show Cluster BV"), this);
	viewShowClusterBVAction->setIcon(QIcon(":/images/showgrid.png"));
	viewShowClusterBVAction->setCheckable(true);

}

void MainWindow::createConnections(){
	// view
	connect(viewShowGridLinesAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowGrid(bool)));
	connect(viewShowBackgroundAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowBackground(bool)));
	connect(viewShowIntersectionAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowIntersection(bool)));
	connect(viewShowUnionBAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowUnionB(bool)));
	connect(viewShowUnionEAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowUnionE(bool)));
	connect(viewShowGradientEAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowGradientE(bool)));
	connect(viewShowLineChartAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowLineChart(bool)));
	connect(viewShowContourLineTruthAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineTruth(bool)));
	connect(viewShowContourLineAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLine(bool)));
	connect(viewShowContourLineSortedAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineSorted(bool)));
	connect(viewShowContourLineSortedSDFAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineSortedSDF(bool)));
	connect(viewShowContourLineResampledAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineResampled(bool)));
	connect(viewShowContourLineSDFAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineSDF(bool)));
	connect(viewShowContourLineMinAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMin(bool)));
	connect(viewShowContourLineMaxAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMax(bool)));
	connect(viewShowContourLineMeanAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMean(bool)));
	connect(viewShowContourLineMedianAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMedian(bool)));
	connect(viewShowContourLineOutlierAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineOutlier(bool)));
	connect(viewShowClusterBSAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowClusterBS(bool)));
	connect(viewShowClusterBVAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowClusterBV(bool)));

	connect(_pDisplayCtrlWidget, SIGNAL(bgFunctionChanged(int)), _view3D, SLOT(updateBgFunction(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(smoothChanged(int)), _view3D, SLOT(updateVarSmooth(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(areasChanged(int)), _view3D, SLOT(updateUncertaintyAreas(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(focusedClusterChanged(int)), _view3D, SLOT(updateFocusedCluster(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(EOFChanged(int)), _view3D, SLOT(updateEOF(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(MemberChanged(int)), _view3D, SLOT(updateMember(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(EnsClusterChanged(int)), _view3D, SLOT(updateEnsCluster(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(EnsClusterLenChanged(int)), _view3D, SLOT(updateEnsClusterLen(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(ContourLevelChanged(int)), _view3D, SLOT(updateContourLevel(int)));
	connect(_pDisplayCtrlWidget, SIGNAL(TimeStepChanged(int)), _pModel, SLOT(updateTimeStep(int)));
	connect(_pModel, SIGNAL(UpdateView()), _view3D, SLOT(onUpdateView()));


}

void MainWindow::onMousePressed(int x, int y){
}
   
void MainWindow::closeEvent(QCloseEvent *event)
{
// 	if (okToClearData()) 
	{
		QSettings settings;
		settings.setValue(ShowGridLines, viewShowGridLinesAction->isChecked());
		settings.setValue(ShowBackground, viewShowBackgroundAction->isChecked());
		settings.setValue(ShowIntersection, viewShowIntersectionAction->isChecked());
		settings.setValue(ShowUnionB, viewShowUnionBAction->isChecked());
		settings.setValue(ShowUnionE, viewShowUnionEAction->isChecked());
		settings.setValue(ShowGradientE, viewShowGradientEAction->isChecked());
		settings.setValue(ShowLineChart, viewShowLineChartAction->isChecked());
		settings.setValue(ShowContourLineTruth, viewShowContourLineTruthAction->isChecked());
		settings.setValue(ShowContourLine, viewShowContourLineAction->isChecked());
		settings.setValue(ShowContourLineMin, viewShowContourLineMinAction->isChecked());
		settings.setValue(ShowContourLineSDF, viewShowContourLineSDFAction->isChecked());
		settings.setValue(ShowContourLineMax, viewShowContourLineMaxAction->isChecked());
		settings.setValue(ShowContourLineSorted, viewShowContourLineSortedAction->isChecked());
		settings.setValue(ShowContourLineSortedSDF, viewShowContourLineSortedSDFAction->isChecked());
		settings.setValue(ShowContourLineResampled, viewShowContourLineResampledAction->isChecked());
		settings.setValue(ShowContourLineMean, viewShowContourLineMeanAction->isChecked());
		settings.setValue(ShowContourLineMedian, viewShowContourLineMedianAction->isChecked());
		settings.setValue(ShowContourLineOutlier, viewShowContourLineOutlierAction->isChecked());
		settings.setValue(ShowClusterBS, viewShowClusterBSAction->isChecked());
		settings.setValue(ShowClusterBV, viewShowClusterBVAction->isChecked());

//		settings.setValue(ShowBeliefEllipse, _pControlWidget->ui.radioButtonBackgroundCluster->isChecked());

		event->accept();
	}
// 	else
// 		event->ignore();
}

