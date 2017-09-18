#include "DisplayCtrlWidget.h"

#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QComboBox>

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
	_pCombBgFunction->addItem("Threshold", 3);
	_pCombBgFunction->addItem("Smoothed", 4);

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

}

void DisplayCtrlWidget::createLayout() {

	QFormLayout *layout = new QFormLayout;
	layout->addRow(tr("Background:"), _pCombBgFunction);
	layout->addRow(tr("Smooth:"), _pSpinBoxSmooth);
	layout->addRow(tr("Areas:"), _pSpinBoxAreas);
	layout->addRow(tr("Focused Cluster:"), _pSpinBoxFocusedCluster);
	setLayout(layout);
}

void DisplayCtrlWidget::createConnections() {
	connect(_pCombBgFunction, SIGNAL(currentIndexChanged(int)), this, SLOT(updateBgFunction(int)));
	connect(_pSpinBoxSmooth, SIGNAL(valueChanged(int)), this, SLOT(updateSmooth(int)));
	connect(_pSpinBoxAreas, SIGNAL(valueChanged(int)), this, SLOT(updateUncertaintyAreas(int)));
	connect(_pSpinBoxFocusedCluster, SIGNAL(valueChanged(int)), this, SLOT(updateFocusedCluster(int)));
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

