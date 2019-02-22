#include "DisplayCtrlWidget.h"

#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QComboBox>

#include <QDebug>
#include "def.h"

DisplayCtrlWidget::DisplayCtrlWidget(QWidget *parent)
	: QWidget(parent)
{
	createWidgets();
	createLayout();
	createConnections();
}

DisplayCtrlWidget::~DisplayCtrlWidget()
{
}

void DisplayCtrlWidget::createWidgets() {

	_pCombBgFunction = new QComboBox;
	_pCombBgFunction->addItem("Mean", 0);
	_pCombBgFunction->addItem("Variance", 1);
	_pCombBgFunction->addItem("Cluster", 2);
	_pCombBgFunction->addItem("VarThreshold", 3);
	_pCombBgFunction->addItem("Smoothed", 4);
	_pCombBgFunction->addItem("DipValue", 5);
	_pCombBgFunction->addItem("DipValueThreshold", 6);
	_pCombBgFunction->addItem("EOF", 7);
	_pCombBgFunction->addItem("Obs", 8);
	_pCombBgFunction->addItem("Error", 9);
	_pCombBgFunction->addItem("SDF", 10);
	_pCombBgFunction->addItem("Line Kernel Density", 11);
	_pCombBgFunction->addItem("Line Kernel Density X", 12);
	_pCombBgFunction->addItem("Line Kernel Density Y", 13);
	_pCombBgFunction->addItem("Line Kernel Density Z", 14);
	_pCombBgFunction->addItem("ContourDensity X", 15);
	_pCombBgFunction->addItem("ContourDensity Y", 16);
	_pCombBgFunction->addItem("ContourDensity Z", 17);
	_pCombBgFunction->addItem("ContourDensity", 18);

	_pSpinBoxSmooth = new QSpinBox;
	_pSpinBoxSmooth->setRange(1, 10);
	_pSpinBoxSmooth->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxSmooth->setValue(1);
	_pSpinBoxSmooth->setSingleStep(1);

	_pSpinBoxAreas = new QSpinBox;
	_pSpinBoxAreas->setRange(1, 6);
	_pSpinBoxAreas->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxAreas->setValue(4);
	_pSpinBoxAreas->setSingleStep(1);

	_pSpinBoxFocusedCluster = new QSpinBox;
	_pSpinBoxFocusedCluster->setRange(0, 5);
	_pSpinBoxFocusedCluster->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxFocusedCluster->setValue(0);
	_pSpinBoxFocusedCluster->setSingleStep(1);

	_pSpinBoxFocusedRegion = new QSpinBox;
	_pSpinBoxFocusedRegion->setRange(0, 5);
	_pSpinBoxFocusedRegion->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxFocusedRegion->setValue(0);
	_pSpinBoxFocusedRegion->setSingleStep(1);


	_pSpinBoxEOF = new QSpinBox;
	_pSpinBoxEOF->setRange(1, g_nEOFLen);
	_pSpinBoxEOF->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxEOF->setValue(1);
	_pSpinBoxEOF->setSingleStep(1);

	_pSpinBoxMember = new QSpinBox;
	_pSpinBoxMember->setRange(0, g_nEnsembles);
	_pSpinBoxMember->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxMember->setValue(0);
	_pSpinBoxMember->setSingleStep(1);

	_pSpinBoxEnsCluster = new QSpinBox;
	_pSpinBoxEnsCluster->setRange(0, g_nEnsClusterLen);
	_pSpinBoxEnsCluster->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxEnsCluster->setValue(0);
	_pSpinBoxEnsCluster->setSingleStep(1);

	_pSpinBoxEnsClusterLen = new QSpinBox;
	_pSpinBoxEnsClusterLen->setRange(1, g_nEnsembles);
	_pSpinBoxEnsClusterLen->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxEnsClusterLen->setValue(1);
	_pSpinBoxEnsClusterLen->setSingleStep(1);


	_pSpinBoxContourLevel = new QSpinBox;
	_pSpinBoxContourLevel->setRange(0, 5);
	_pSpinBoxContourLevel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxContourLevel->setValue(0);
	_pSpinBoxContourLevel->setSingleStep(1);


	_pSpinBoxTimeStep = new QSpinBox;
	_pSpinBoxTimeStep->setRange(g_nTimeStart, g_nTimeEnd);
	_pSpinBoxTimeStep->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxTimeStep->setValue(g_nTimeStart);
	_pSpinBoxTimeStep->setSingleStep(6);

}

void DisplayCtrlWidget::createLayout() {

	QFormLayout *layout = new QFormLayout;
	layout->addRow(tr("Background:"), _pCombBgFunction);
	layout->addRow(tr("Smooth:"), _pSpinBoxSmooth);
	layout->addRow(tr("Areas:"), _pSpinBoxAreas);
	layout->addRow(tr("Focused Cluster:"), _pSpinBoxFocusedCluster);
	layout->addRow(tr("Focused Region:"), _pSpinBoxFocusedRegion);
	layout->addRow(tr("EOF:"), _pSpinBoxEOF);
	layout->addRow(tr("Member:"), _pSpinBoxMember);
	layout->addRow(tr("Ens Cluster:"), _pSpinBoxEnsCluster);
	layout->addRow(tr("Number of Ens Clusters:"), _pSpinBoxEnsClusterLen);
	layout->addRow(tr("Contour Level:"), _pSpinBoxContourLevel);
	layout->addRow(tr("Time Step:"), _pSpinBoxTimeStep);
	setLayout(layout);
}

void DisplayCtrlWidget::createConnections() {
	connect(_pCombBgFunction, SIGNAL(currentIndexChanged(int)), this, SLOT(updateBgFunction(int)));
	connect(_pSpinBoxSmooth, SIGNAL(valueChanged(int)), this, SLOT(updateSmooth(int)));
	connect(_pSpinBoxAreas, SIGNAL(valueChanged(int)), this, SLOT(updateUncertaintyAreas(int)));
	connect(_pSpinBoxFocusedCluster, SIGNAL(valueChanged(int)), this, SLOT(updateFocusedCluster(int)));
	connect(_pSpinBoxEOF, SIGNAL(valueChanged(int)), this, SLOT(updateEOF(int)));
	connect(_pSpinBoxMember, SIGNAL(valueChanged(int)), this, SLOT(updateMember(int)));
	connect(_pSpinBoxEnsCluster, SIGNAL(valueChanged(int)), this, SLOT(updateEnsCluster(int)));
	connect(_pSpinBoxEnsClusterLen, SIGNAL(valueChanged(int)), this, SLOT(updateEnsClusterLen(int)));
	connect(_pSpinBoxContourLevel, SIGNAL(valueChanged(int)), this, SLOT(updateContourLevel(int)));
	connect(_pSpinBoxTimeStep, SIGNAL(valueChanged(int)), this, SLOT(updateTimeStep(int)));
}

void DisplayCtrlWidget::updateBgFunction(int nBgFunction)
{
	emit bgFunctionChanged(nBgFunction);
}

void DisplayCtrlWidget::updateSmooth(int nSmooth)
{
	emit smoothChanged(nSmooth);
}

void DisplayCtrlWidget::updateUncertaintyAreas(int nAreas)
{
	_pSpinBoxFocusedCluster->setRange(0, nAreas-1);
	emit areasChanged(nAreas);
}

void DisplayCtrlWidget::updateFocusedCluster(int nFocusedCluster) {
	emit focusedClusterChanged(nFocusedCluster);
}


void DisplayCtrlWidget::updateEOF(int nEOF) {

	emit EOFChanged(nEOF);
}
void DisplayCtrlWidget::updateMember(int nMember) {

	emit MemberChanged(nMember);
}
void DisplayCtrlWidget::updateEnsCluster(int nEnsCluster) {

	emit EnsClusterChanged(nEnsCluster);
}
void DisplayCtrlWidget::updateEnsClusterLen(int nEnsClusterLen) {

	emit EnsClusterLenChanged(nEnsClusterLen);
}
void DisplayCtrlWidget::updateContourLevel(int nLevel) {

	emit ContourLevelChanged(nLevel);
}
void DisplayCtrlWidget::updateTimeStep(int nTS) {

	emit TimeStepChanged(nTS);
}