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

}

void DisplayCtrlWidget::createLayout() {

	QFormLayout *layout = new QFormLayout;
	layout->addRow(tr("Background:"), _pCombBgFunction);
	layout->addRow(tr("Smooth:"), _pSpinBoxSmooth);
	layout->addRow(tr("Areas:"), _pSpinBoxAreas);
	layout->addRow(tr("Focused Cluster:"), _pSpinBoxFocusedCluster);
	layout->addRow(tr("Focused Region:"), _pSpinBoxFocusedRegion);
	layout->addRow(tr("EOF:"), _pSpinBoxEOF);
	setLayout(layout);
}

void DisplayCtrlWidget::createConnections() {
	connect(_pCombBgFunction, SIGNAL(currentIndexChanged(int)), this, SLOT(updateBgFunction(int)));
	connect(_pSpinBoxSmooth, SIGNAL(valueChanged(int)), this, SLOT(updateSmooth(int)));
	connect(_pSpinBoxAreas, SIGNAL(valueChanged(int)), this, SLOT(updateUncertaintyAreas(int)));
	connect(_pSpinBoxFocusedCluster, SIGNAL(valueChanged(int)), this, SLOT(updateFocusedCluster(int)));
	connect(_pSpinBoxFocusedRegion, SIGNAL(valueChanged(int)), this, SLOT(updateFocusedRegion(int)));
	connect(_pSpinBoxEOF, SIGNAL(valueChanged(int)), this, SLOT(updateEOF(int)));
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

void DisplayCtrlWidget::updateFocusedRegion(int nFocusedRegion) {
//	qDebug() << "updateFocusedRegion";

	emit focusedRegionChanged(nFocusedRegion);
}

void DisplayCtrlWidget::updateEOF(int nEOF) {
	//	qDebug() << "updateFocusedRegion";

	emit EOFChanged(nEOF);
}