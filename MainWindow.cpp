#include <QtGui>



#include <QtWidgets/QGridLayout>
#include <QLabel>
#include <QMenu>
#include <QToolBar>
#include <QMenuBar>
#include <QSet>
#include <QMessageBox>
#include <QDockWidget>




#include "MainWindow.h"
#include "MainView.h"

#include "CustomTableModel.h"
#include "DetailWindow.h"
#include "TableWidget.h"

#include "myglwidget.h"
#include "MeteModel.h"
#include "ClusteringWidget.h"


// using the original temperature data
#define TEMPERATURE

#ifdef TEMPERATURE
const char* g_strFileT = "../../data/data5.dat";
const char* g_strFileE = "../../data/data6.dat";
const char* g_strFileB = "../../data/data8/pro-2008-01-01-96h-t2-%1.dat";
#else
const char* g_strFileT = "../../data/data5.dat";
const char* g_strFileE = "../../data/data6.dat";
const char* g_strFileB = "../../data/data9/pro-2016-10-01-120h-h500-%1.dat";
#endif

// 	char *g_strFileE = "../../data/data10/pre-mod-cma-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-cptec-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-eccc-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-ecmwf-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-jma-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-kma-20160802-00-96.txt";
// 	char *g_strFileE = "../../data/data10/pre-mod-ncep-20160802-00-96.txt";

const QString ShowGridLines("ShowGridLines");
const QString ShowBackground("ShowBackground");
const QString ShowIntersection("ShowIntersection");
const QString ShowUnionB("ShowUnionB");
const QString ShowUnionE("ShowUnionE");
const QString ShowGradientE("ShowGradientE");
const QString ShowLineChart("ShowLineChart");
const QString ShowContourLineTruth("ShowContourLineTruth");
const QString ShowContourLine("ShowContourLine");
const QString ShowContourLineMin("ShowContourLineMin");
const QString ShowContourLineMax("ShowContourLineMax");
const QString ShowContourLineMean("ShowContourLineMean");
const QString ShowClusterBS("ShowClusterBS");
const QString ShowClusterBV("ShowClusterBV");
const QString ShowBeliefEllipse("ShowBeliefEllipse");


MainWindow::~MainWindow(){
	delete[] _gridDataMinB;
	delete[] _gridDataMeanB;
	delete[] _gridDataMaxB;
	delete[] _dataB;
	delete[] _gridLabels;
	delete[] _gridLabelsF;
	delete[]_arrLatVarB;
	for (int i = 0; i < g_gradient_l; i++)
	{
		delete[](_gridDataMinB_3[i]);
		delete[](_gridDataMaxB_3[i]);
		delete[](_gridDataMaxGradientB[i]);
		delete[](_gridDataMinGradientB[i]);
	}

	delete[] _gridIntersectionMaxB;
	delete[] _gridIntersectionMinB;

	delete[] _gridDataUnionMaxB;
	delete[] _gridDataUnionMinB;
	delete[] _gridDataMeanMeanB;
	delete[] _gridDataMeanVarB;
	delete[] _gridDataVarMeanB;
	delete[] _normalizedMVB;
	delete[] _normalizedVMB;
	delete[] _dataT;
	delete[] _arrIntersections;

	if (_pModel) delete _pModel;
}

MainWindow::MainWindow()
{
	resize(QSize(800, 600));
	_gridDataMinB = new double[g_focus_l * g_lenEnsembles];
	_gridDataMeanB = new double[g_focus_l * g_lenEnsembles];
	_gridDataMaxB = new double[g_focus_l * g_lenEnsembles];
	for (int i = 0; i < g_gradient_l;i++)
	{

		_gridDataMinB_3[i] = new double[g_focus_l * g_lenEnsembles];
		_gridDataMaxB_3[i] = new double[g_focus_l * g_lenEnsembles];
		_gridDataMaxGradientB[i] = new double[g_focus_l];
		_gridDataMinGradientB[i] = new double[g_focus_l];
	}

	_dataB = new double[g_focus_l * 4 * g_lenEnsembles];
	_gridLabels = new int[g_focus_l];
	_gridLabelsF = new double[g_focus_l];
	_arrLatVarB = new double[g_focus_h];

	
	_gridIntersectionMaxB = new double[g_focus_l];
	_gridIntersectionMinB = new double[g_focus_l];

	_gridDataUnionMaxB = new double[g_focus_l];
	_gridDataMeanMeanB = new double[g_focus_l];
	_gridDataUnionMinB = new double[g_focus_l];
	_gridDataMeanVarB = new double[g_focus_l];
	_gridDataVarMeanB = new double[g_focus_l];

	_normalizedMVB = new double[g_focus_l];
	_normalizedVMB = new double[g_focus_l];

	_dataT = new double[g_focus_l];


	_pModel = MeteModel::CreateModel();

	_arrIntersections = new EnsembleIntersections[g_focus_l];

// 	double fValue = 0273.16 + g_temperature + .0001;

//	double fValue = g_temperature + .0001;


// 	readDataB();
// 	readDataT();

// 	doStatisticE();
// 	doStatisticB();
// 	doStatistic3B();
// 	if (!_bClustering){
// // 		readLabels();
// 		doSegment();
// 	}
// 	else{
// 		doCluster();
// 	}
// 

// 	for (int l = 0; l < 1; l++)
// 	{
// 		// 1.generate contours of this ensemble
// 		const double *data = _dataE + l*g_focus_l;
// 		generateContour(data, _listContourE[0], 2);
// 		generateContour(data, _listContourE[0], 5);
// 		generateContour(data, _listContourE[0], 10);
// 		generateContour(data, _listContourE[0], 30);
// 		generateContour(data, _listContourE[0], 50);
// 		generateContour(data, _listContourE[0], 70);
// 		generateContour(data, _listContourE[0], 90);
// 	}

// 	generateContour(_gridUnionMinE, _listContourUnionMinE[0], 2);
// 	generateContour(_gridUnionMaxE, _listContourUnionMaxE[0], 2);

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
	viewShowContourLineMinAction->setChecked(settings.value(ShowContourLineMin, true).toBool());
	viewShowContourLineMaxAction->setChecked(settings.value(ShowContourLineMax, true).toBool());
	viewShowContourLineMeanAction->setChecked(settings.value(ShowContourLineMean, true).toBool());
	viewShowClusterBSAction->setChecked(settings.value(ShowClusterBS, true).toBool());
	viewShowClusterBVAction->setChecked(settings.value(ShowClusterBV, true).toBool());

	_pControlWidget->ui.radioButtonBackgroundCluster->setChecked(settings.value(ShowClusterBV, true).toBool());




    setWindowTitle(tr("Ensembles"));
}

void MainWindow::createSceneAndView(){
	_view3D = new MyGLWidget;
	_view3D->SetDataT(_dataT);
// 	_view3D->SetDataT(_gridDataMeanMeanB);
// 	_view3D->SetDataT(_gridDataMeanVarB);
// 	_view3D->SetDataT(_gridDataVarMeanB);
// 	_view3D->SetDataT(_varE);
// 	_view3D->SetVar(_gridDataVarMeanB, _gridDataMeanVarB);
	_view3D->SetVar(_normalizedVMB, _normalizedMVB);
	_view3D->SetLabels(_gridLabels);

	_view3D->SetDataB(_dataB);
	_view3D->SetModelE(_pModel);

	// 	_view3D->SetData(_gridDataVarB);
	// 	_view3D->SetData(_varE);
	_view3D->SetContourTruth(_listContourT);
// 	_view3D->SetContourMin(_listContourIntersectionMinB);
// 	_view3D->SetContourMax(_listContourIntersectionMaxB);
	// 	_view3D->SetContourMin(_listContourIntersectionMinE);
	// 	_view3D->SetContourMax(_listContourIntersectionMaxE);
	_view3D->SetContourIntervalB(_listContourIntersectionMinB, _listContourIntersectionMaxB);
// 	_view3D->SetContourIntervalE(_listContourUnionMinE, _listContourUnionMaxE);
	// 	_view3D->SetContourMean(_listContourMeanB);
	_view3D->SetContourMean(_listContourEnsembleMeanB);
//	_view3D->SetContour(_listContourE);
// 	_view3D->SetContourMean(_listContourE1);
// 	_view3D->SetContour(_listContourE2);
	_view3D->SetMultiStatistic(_listContourEnsembleMin_3, _listContourEnsembleMax_3);
	_view3D->SetIntersections(_arrIntersections);


	setCentralWidget(_view3D);
}

void MainWindow::createDockWidgets() {
	setDockOptions(QMainWindow::AnimatedDocks);
	QDockWidget::DockWidgetFeatures features =
		QDockWidget::DockWidgetMovable |
		QDockWidget::DockWidgetFloatable;

	_pControlWidget = new ControlWidget();


	QDockWidget *controlDockWidget = new QDockWidget(
		tr("Control"), this);
	controlDockWidget->setFeatures(features);
	controlDockWidget->setWidget(_pControlWidget);
	addDockWidget(Qt::RightDockWidgetArea, controlDockWidget);


	_pWidgetClustering = new ClusteringWidget();
	QDockWidget *pDockWidgetClustering = new QDockWidget(tr("Clustering"), this);
	pDockWidgetClustering->setFeatures(features);
	pDockWidgetClustering->setWidget(_pWidgetClustering);
	addDockWidget(Qt::LeftDockWidgetArea, pDockWidgetClustering);

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
		<< viewShowContourLineMinAction
		<< viewShowContourLineMaxAction
		<< viewShowContourLineMeanAction
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

	viewShowContourLineMinAction = new QAction(tr("Show Contour of Min"), this);
	viewShowContourLineMinAction->setIcon(QIcon(":/images/min.png"));
	viewShowContourLineMinAction->setCheckable(true);

	viewShowContourLineMaxAction = new QAction(tr("Show Contour of Max"), this);
	viewShowContourLineMaxAction->setIcon(QIcon(":/images/max.png"));
	viewShowContourLineMaxAction->setCheckable(true);

	viewShowContourLineMeanAction = new QAction(tr("Show Contour of Mean"), this);
	viewShowContourLineMeanAction->setIcon(QIcon(":/images/b.png"));
	viewShowContourLineMeanAction->setCheckable(true);


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
	connect(viewShowContourLineMinAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMin(bool)));
	connect(viewShowContourLineMaxAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMax(bool)));
	connect(viewShowContourLineMeanAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowContourLineMean(bool)));
	connect(viewShowClusterBSAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowClusterBS(bool)));
	connect(viewShowClusterBVAction, SIGNAL(toggled(bool)), _view3D, SLOT(viewShowClusterBV(bool)));


	connect(_pControlWidget->ui.radioButtonBackgroundMean, SIGNAL(clicked(bool)), this, SLOT(onSelectBackgroundMean(bool)));
	connect(_pControlWidget->ui.radioButtonBackgroundVari, SIGNAL(clicked(bool)), this, SLOT(onSelectBackgroundVari(bool)));
	connect(_pControlWidget->ui.radioButtonBackgroundCluster, SIGNAL(clicked(bool)), this, SLOT(onSelectBackgroundCluster(bool)));
	connect(_pControlWidget->ui.radioButtonBackgroundSDF, SIGNAL(clicked(bool)), this, SLOT(onSelectBackgroundSDF(bool)));
	connect(_pControlWidget->ui.radioButtonBackgroundVarSmooth, SIGNAL(clicked(bool)), this, SLOT(onSelectBackgroundVarSmooth(bool)));


	connect(_pControlWidget->ui.spinBoxVarSmooth, SIGNAL(valueChanged(int)), this, SLOT(updateVarSmooth(int)));


	connect(_pControlWidget->ui.checkBoxShowEllipse, SIGNAL(toggled(bool)), _view3D, SLOT(onCheckShowBeliefEllipse(bool)));


	connect(_pWidgetClustering, SIGNAL(minPtsChanged(int)), _view3D, SLOT(updateMinPts(int)));
	connect(_pWidgetClustering, SIGNAL(epsChanged(double)), _view3D, SLOT(updateEps(double)));


}

void MainWindow::readDataB(){
	for (int l = 0; l < g_lenEnsembles; l++)
	{
// 		QString fileName = QString("../../data/data4/pro-2008-01-01-96h-t2-%1.dat").arg(QString::number(l + 1));
// 		QString fileName = QString("../../data/data7/pro-2008-01-01-96h-t2-%1.dat").arg(QString::number(l + 1));
// 		QString fileName = QString("../../data/data8/pro-2008-01-01-96h-t2-%1.dat").arg(QString::number(l + 1));
		QString fileName = QString(g_strFileB).arg(QString::number(l + 1));


		QFile file(fileName);
		if (!file.open(QIODevice::ReadOnly)) {
			return;
		}
		char* temp = new char[file.size()];
		file.read(temp, file.size());
		double* f = (double*)temp;


		for (int iii = 0; iii < g_grid_height; iii++)
		{
// 			int i = g_grid_height-1 - iii;
			int i = iii;
			int ii = i - g_focus_y;
			for (int j = 0; j < g_grid_width; j++)
			{
				int jj = j - g_focus_x;
				if (i >= g_focus_y&&i < g_focus_y + g_focus_h_raw&&j >= g_focus_x&&j < g_focus_x + g_focus_w_raw){
					// in the scope
					for (int s = 0; s < 4; s++){
						if (ii%g_nSpace == 0 && jj%g_nSpace == 0){
							// on the grid
							int ij = (ii/g_nSpace)*g_focus_w + jj/g_nSpace;
							double fValue = *f;
							_dataB[l*g_focus_l * 4 + (ij)* 4 + s] = *f;
						}
						f++;
					}
				}
				else{
					f += 4;
				}
			}
		}
		file.close();
		delete[] temp;
	}

}

void MainWindow::readDataT(){
// 	char* fileName = "../../data/data5.dat";
	QFile file(g_strFileT);
	if (!file.open(QIODevice::ReadOnly)) {
		return;
	}

	int nLen = file.size();
	char* temp = new char[nLen];
	file.read(temp, file.size());
	double* f = (double*)temp;

	for (int i = 0; i < g_grid_height; i++)
	{
		int ii = i - g_focus_y;
		for (int j = 0; j < g_grid_width; j++)
		{
			int jj = j - g_focus_x;
			if (i >= g_focus_y&&i < g_focus_y + g_focus_h_raw&&j >= g_focus_x&&j < g_focus_x + g_focus_w_raw){
				// in the scope
				if (ii%g_nSpace==0&&jj%g_nSpace==0){
					// on the grid
					_dataT[(ii / g_nSpace)*g_focus_w + jj/g_nSpace] = *f;
				}
			}
			f++;
		}
	}

	file.close();
	delete[] temp;
}

void MainWindow::doStatisticB(){
	// 1.calculate min,mean, and max for each grid point, for each ensemble result
	for (int l = 0; l < g_lenEnsembles; l++){
		double *data = _dataB + l * 4 * g_focus_l;
		for (int g = 0; g < g_focus_l; g++){
			// for each grid point
			double mb = data[0];
			double mv = data[1];
			double m = data[2];
			double v = data[3];
			double fMin = m - v;
			double fMax = m + v;
			double fMean = m;
			fMin = fMin*mv + mb;
			fMax = fMax*mv + mb;
			fMean = fMean*mv + mb;

// 			if (g == 20 * g_focus_w + 43){
// 				qDebug() << mb << "\t" << mb << "\t" << m << "\t" << v;
// 			}
// 			qDebug() << fMean;
			_gridDataMinB[l*g_focus_l + g] = fMin;
			_gridDataMeanB[l*g_focus_l + g] = fMean;
			_gridDataMaxB[l*g_focus_l + g] = fMax;
			data += 4;
		}
	}
	// 2.for each point, calculate the union
	for (int g = 0; g < g_focus_l; g++){
		double fMaxG = -1000;
		double fMinG = 1000;
		double fMeanG = 0;
		double fVarG = 0;
		for (int l = 0; l < g_lenEnsembles; l++)
		{
			
			double fMin = _gridDataMinB[l*g_focus_l + g];
			double fMean = _gridDataMeanB[l*g_focus_l + g];
			double fMax = _gridDataMaxB[l*g_focus_l + g];
			if (fMax > fMaxG) fMaxG = fMax;
			if (fMin < fMinG) fMinG = fMin;
			fMeanG += fMean;
			fVarG += fMax - fMean;
		}

		fMeanG /= g_lenEnsembles;
		fVarG /= g_lenEnsembles;
		_gridDataUnionMinB[g] = fMinG;
		_gridDataUnionMaxB[g] = fMaxG;
		_gridDataMeanMeanB[g] = fMeanG;
		_gridDataMeanVarB[g] = fVarG;
		// calculate variance
		double fVariance = 0;
		for (int l = 0; l < g_lenEnsembles; l++){
			double fMean = _gridDataMeanB[l*g_focus_l + g];
			double fBias = fMeanG - fMean;
// 			fVariance += fBias*fBias;
			fVariance += abs(fBias);
		}
// 		_gridDataVarMeanB[g] = sqrt(fVariance / g_lenEnsembles);	
		_gridDataVarMeanB[g] = fVariance / g_lenEnsembles;

// 		qDebug() << _gridDataMeanVarB[g] << endl;
	}
	// 3.for each point, calculate the intersection
	for (int g = 0; g < g_focus_l; g++){
		EnsembleIntersections intersection;
		for (int l = 0; l < g_lenEnsembles; l++)
		{
			double fMin = _gridDataMinB[l*g_focus_l + g];
			double fMax = _gridDataMaxB[l*g_focus_l + g];
			intersection.Intersect(fMin, fMax);
		}
		// remove invalid intersection
		int nLen = intersection._listIntersection.length();
		for (int i = nLen - 1; i >= 0; i--)
		{
			if (intersection._listIntersection[i]._fMax == intersection._listIntersection[i]._fMin)
			{
				intersection._listIntersection.removeAt(i);
			}
		}
		_arrIntersections[g] = intersection;
	}

	// 4.for each point, calculate the most overlapped intersection
	for (int g = 0; g < g_focus_l; g++)
	{
		EnsembleIntersections intersections = _arrIntersections[g];
// 		intersections.GetMostOverlapped(_gridIntersectionMin[g], _gridIntersectionMax[g]);
		intersections.GetMostOverlapped(_gridIntersectionMinB[g], _gridIntersectionMaxB[g], g_nThreshold);
	}

	// 5.calculate the mean variance of each latitude
	for (int i = 0; i < g_focus_h;i++)
	{
		double fVar = 0;
		for (int j = 0; j < g_focus_w;j++)
		{
			int index = i*g_focus_w + j;
			for (int l = 0; l < g_lenEnsembles;l++)
			{
				fVar += _gridDataMaxB[l*g_focus_l + index] - _gridDataMeanB[l*g_focus_l + index];
			}
		}
		_arrLatVarB[i] = fVar;
	}
// 	qDebug() << "Raw";
// 	for (int i = 0; i < g_focus_h;i++)
// 	{
// 		qDebug() << _arrLatVarB[i];
// 	}
// 	// calculate regression coefficient
// 	double ex = 0;
// 	double ey = 0;
// 	for (int i = 0; i < g_focus_h;i++)
// 	{
// 		ex += i;
// 		ey += arrVar[i];
// 	}
// 	ex /= g_focus_h;
// 	ey /= g_focus_h;
// 	double fNumerator = 0;
// 	double fDenominator = 0;
// 	for (int i = 0; i < g_focus_h;i++)
// 	{
// 		double fXBias = (i - ex);
// 		fNumerator += fXBias*(arrVar[i] - ey);
// 		fDenominator += fXBias*fXBias;
// 	}
// 	double b = fNumerator / fDenominator;
// 	double a = ey - b*ex;
// 	qDebug() << "ex:\t" << ex;
// 	qDebug() << "ey:\t" << ey;
// 	qDebug() << "b:\t" << b;
// 	qDebug() << "a:\t" << a;


	// 6.normalize vm and mv
	double fVMMin = 1000;
	double fMVMin = 1000;
	double fVMMax = -1000;
	double fMVMax = -1000;
	for (int i = 0; i < g_focus_l; i++){
		if (_gridDataMeanVarB[i] > fMVMax) fMVMax = _gridDataMeanVarB[i];
		if (_gridDataMeanVarB[i] < fMVMin) fMVMin = _gridDataMeanVarB[i];
		if (_gridDataVarMeanB[i] > fVMMax) fVMMax = _gridDataVarMeanB[i];
		if (_gridDataVarMeanB[i] < fVMMin) fVMMin = _gridDataVarMeanB[i];
	}
	double fVMRange = fVMMax - fVMMin;
	double fMVRange = fMVMax - fMVMin;
	for (int i = 0; i < g_focus_l;i++)
	{
		_normalizedVMB[i] = (_gridDataVarMeanB[i] - fVMMin) / fVMRange;
		_normalizedMVB[i] = (_gridDataMeanVarB[i] - fMVMin) / fMVRange;
	}
}

void MainWindow::doStatistic3B(){
	for (int s = 0; s < g_gradient_l; s++)
	{
		double fBelief = s / 20.0;
		// 1.calculate min,mean, and max for each grid point, for each ensemble result
		for (int l = 0; l < g_lenEnsembles; l++){
			double *data = _dataB + l * 4 * g_focus_l;
			for (int g = 0; g < g_focus_l; g++){
				// for each grid point
				double mb = data[0];
				double mv = data[1];
				double m = data[2];
				double v = data[3];
				double fMin = m - v*fBelief;
				double fMax = m + v*fBelief;
				double fMean = m;
				fMin = fMin*mv + mb;
				fMax = fMax*mv + mb;
				fMean = fMean*mv + mb;
				_gridDataMinB_3[s][l*g_focus_l + g] = fMin;
				_gridDataMaxB_3[s][l*g_focus_l + g] = fMax;

				data += 4;
			}
		}
		// 2.for each point, calculate the union
		for (int g = 0; g < g_focus_l; g++){
			double fMaxG = -1000;
			double fMinG = 1000;
			for (int l = 0; l < g_lenEnsembles; l++)
			{
				double fMin = _gridDataMinB_3[s][l*g_focus_l + g];
				double fMax = _gridDataMaxB_3[s][l*g_focus_l + g];
				if (fMax>fMaxG) fMaxG = fMax;
				if (fMin < fMinG) fMinG = fMin;
			}
			_gridDataMinGradientB[s][g] = fMinG;
			_gridDataMaxGradientB[s][g] = fMaxG;
		}
	}
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
		settings.setValue(ShowContourLineMax, viewShowContourLineMaxAction->isChecked());
		settings.setValue(ShowContourLineMean, viewShowContourLineMeanAction->isChecked());
		settings.setValue(ShowClusterBS, viewShowClusterBSAction->isChecked());
		settings.setValue(ShowClusterBV, viewShowClusterBVAction->isChecked());

		settings.setValue(ShowBeliefEllipse, _pControlWidget->ui.radioButtonBackgroundCluster->isChecked());

		event->accept();
	}
// 	else
// 		event->ignore();
}

void MainWindow::buildDistanceMatrix_ensemble(double** disMatrix){
	for (int i = 0; i < g_focus_l; i++)
	{
		disMatrix[i] = new double[g_focus_l];
	}
	// vector at point i and j
	double dataI[g_lenEnsembles];
	double dataJ[g_lenEnsembles];
	// 1.calculate the distance
	for (int i = 0; i < g_focus_l; i++)
	{
		int latI = i / g_focus_w;			// latitude of index i
		int lonI = i % g_focus_w;
		for (int l = 0; l < g_lenEnsembles; l++)
		{
			// variance
			dataI[l] = (_gridDataMaxB[l * g_focus_l + i] - _gridDataMeanB[l * g_focus_l + i]) / _arrLatVarB[latI];
		}
		for (int j = 0; j < i; j++)
		{
			int latJ = j / g_focus_w;			// latitude of index j
			int lonJ = j % g_focus_w;
			for (int l = 0; l < g_lenEnsembles; l++)
			{
				// variance
				dataJ[l] = (_gridDataMaxB[l * g_focus_l + j] - _gridDataMeanB[l * g_focus_l + j]) / _arrLatVarB[latJ];
			}
			double fResult = 0;
			for (int l = 0; l < g_lenEnsembles; l++)
			{
				double fBias = dataI[l] - dataJ[l];
				fResult += fBias*fBias;
			}
			double biasLon = lonI - lonJ;
			double biasLat = latI - latJ;
			// 			disMatrix[i][j] = sqrt(fResult) + sqrt(biasLon*biasLon + biasLat*biasLat)*0.0001;
			disMatrix[i][j] = sqrt(fResult);
			// 			qDebug() << sqrt(fResult);
			// 			qDebug() << disMatrix[i][j];
		}
	}
}

void MainWindow::buildDistanceMatrix_vb(double** disMatrix){
	for (int i = 0; i < g_focus_l; i++)
	{
		disMatrix[i] = new double[g_focus_l];
	}
	// 1.calculate the distance
	for (int i = 0; i < g_focus_l; i++)
	{
		int latI = i / g_focus_w;			// latitude of index i
		int lonI = i % g_focus_w;
// 		double fVMI = _gridDataVarMeanB[i];
// 		double fMVI = _gridDataMeanVarB[i];
		double fVMI = _normalizedVMB[i];
		double fMVI = _normalizedMVB[i];
		for (int j = 0; j < i; j++)
		{
			int latJ = j / g_focus_w;			// latitude of index j
			int lonJ = j % g_focus_w;
// 			double fVMJ = _gridDataVarMeanB[j];
// 			double fMVJ = _gridDataMeanVarB[j];
			double fVMJ = _normalizedVMB[j];
			double fMVJ = _normalizedMVB[j];
			double fBias1 = fVMJ - fVMI;
			double fBias2 = fMVJ - fMVI;
			double biasLon = abs(lonI - lonJ);
			double biasLat = abs(latI - latJ);

			// four kinds of distance calculation methods
// 			disMatrix[i][j] = sqrt(fBias1*fBias1 + fBias2*fBias2) + sqrt(biasLat*biasLat + biasLon*biasLon)*0.5;
// 			disMatrix[i][j] = qMax(biasLon, biasLat);
// 			disMatrix[i][j] = sqrt(biasLon*biasLon + biasLat*biasLat);
			disMatrix[i][j] = sqrt(fBias1*fBias1 + fBias2*fBias2);
		}
	}
}

void MainWindow::doSegment(){
	// 1.calculate the distance
	double fThreshold = 0.5;
	for (int i = 0; i < g_focus_l; i++)
	{
		double fVMI = _normalizedVMB[i];
		double fMVI = _normalizedMVB[i];
		_gridLabels[i] = 0;
		if (fVMI > fThreshold) _gridLabels[i] += 1;
		if (fMVI > fThreshold) _gridLabels[i] += 2;
	}
}

void MainWindow::readLabels(){
	char* fileName = "Label.txt";
	QFile caFile(fileName);
	caFile.open(QIODevice::ReadOnly | QIODevice::Text);

	if (!caFile.isOpen()){
		qDebug() << "- Error, unable to open" << "outputFilename" << "for output";
	}
	QTextStream inStream(&caFile);
	int nLabelLen = 0;
	int nTemp = 0;
	for (int i = 0; i < g_focus_l; i++)
	{
		inStream >> _gridLabels[i];
	}
	caFile.close();
	// find the top 20 clusters from the 1000
	int arrCount[1000];
	for (int i = 0; i < 1000;i++)
	{
		arrCount[i] = 0;
	}

	for (int i = 0; i < g_focus_l; i++)
	{
		arrCount[_gridLabels[i]]++;
	}
	QList<int> listTop20;
	for (int i = 0; i < 20;i++)
	{
		int nMax = 0;
		int nIndex;
		for (int j = 0; j < 1000; j++)
		{
			if (arrCount[j]>nMax){
				nMax = arrCount[j];
				nIndex = j;
			}
		}
		qDebug() << nIndex << "\t" << nMax;
		listTop20.push_back(nIndex);
		arrCount[nIndex] = 0;
	}
// 	for each (int n in listTop20)
// 	{
// 		qDebug() << n;
// 	}
	// refine the labels
	// bias between the label and its real label
// 	int nBias = 0;
// 	for (int i = 0; i < g_focus_l;i++)
// 	{
// 		// if i is a valid label
// 		bool bValidLabel = false;
// 		for (int j = 0; j < g_focus_l; j++)
// 		{
// 			if (_gridLabels[j] == i){
// 				_gridLabels[j] -= nBias;
// 				bValidLabel = true;
// 			}
// 		}
// 		if (!bValidLabel) nBias++;
// 	}
// 	for (int i = 0; i < g_focus_l;i++)
// 	{
// 		_gridLabelsF[i] = _gridLabels[i];
// 	}

	{
		// print them
		QFile caFile("LabelRefine.txt");
		caFile.open(QIODevice::WriteOnly | QIODevice::Text);

		if (!caFile.isOpen()){
			qDebug() << "- Error, unable to open" << "outputFilename" << "for output";
		}
		QTextStream outStream(&caFile);

		QList<int> list;
		for (int i = 0; i < g_focus_l; i++)
		{
			list.push_back(_gridLabels[i]);
		}
		qSort(list);
		for each (int n in list)
		{
			outStream << n << endl;
		}

		caFile.close();
	}
}

void MainWindow::onSelectBackgroundMean(bool bChecked) {
	_pModel->SetBgFunctionMean(MeteModel::bg_mean);
	_view3D->ReloadTexture();
}

void MainWindow::onSelectBackgroundVari(bool bChecked) {
	_pModel->SetBgFunctionMean(MeteModel::bg_vari);
	_view3D->ReloadTexture();

}

void MainWindow::onSelectBackgroundCluster(bool bChecked) {
	_pModel->SetBgFunctionMean(MeteModel::bg_cluster);
	_view3D->ReloadTexture();
}

void MainWindow::onSelectBackgroundSDF(bool bChecked) {
	_pModel->SetBgFunctionMean(MeteModel::bg_sdf);
	_view3D->ReloadTexture();
}



void MainWindow::onSelectBackgroundVarSmooth(bool bChecked) {
	_pModel->SetBgFunctionMean(MeteModel::bg_vari_smooth);
	_view3D->ReloadTexture();
}

void MainWindow::updateVarSmooth(int nSmooth) {
	_pModel->SetSmooth(nSmooth);
	_view3D->ReloadTexture();
}