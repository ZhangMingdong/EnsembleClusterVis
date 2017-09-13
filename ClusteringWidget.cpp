#include "ClusteringWidget.h"

#include <QFormLayout>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QComboBox>

ClusteringWidget::ClusteringWidget(QWidget *parent)
	: QWidget(parent)
{
	createWidgets();
	createLayout();
	createConnections();
}

ClusteringWidget::~ClusteringWidget()
{
}

void ClusteringWidget::createWidgets() {

	_pCombMethod = new QComboBox;
	_pCombMethod->addItem("AHC", 0);
	_pCombMethod->addItem("K-Means", 1);
	_pCombMethod->addItem("DBSCAN", 2);

	_pSpinBoxMinPts = new QSpinBox;
	_pSpinBoxMinPts->setRange(1,1000);
	_pSpinBoxMinPts->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxMinPts->setValue(10);
	_pSpinBoxMinPts->setSingleStep(1);

	_pSpinBoxEps = new QDoubleSpinBox;
	_pSpinBoxEps->setRange(1.0, 20.0);
	_pSpinBoxEps->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
	_pSpinBoxEps->setValue(1.0);
	_pSpinBoxEps->setSingleStep(1.0);
	_pSpinBoxEps->setDecimals(1);
}

void ClusteringWidget::createLayout() {

	QFormLayout *layout = new QFormLayout;
	layout->addRow(tr("Method:"), _pCombMethod);
	layout->addRow(tr("MinPts:"), _pSpinBoxMinPts);
	layout->addRow(tr("Eps:"), _pSpinBoxEps);
	setLayout(layout);
}

void ClusteringWidget::createConnections() {
	connect(_pCombMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(updateMethod(int)));
	connect(_pSpinBoxMinPts, SIGNAL(valueChanged(int)), this, SLOT(updateMinPts(int)));
	connect(_pSpinBoxEps, SIGNAL(valueChanged(double)), this, SLOT(updateEps(double)));
}

void ClusteringWidget::updateMethod(int method)
{
	emit methodChanged(method);
}

void ClusteringWidget::updateMinPts(int minPts)
{
	emit minPtsChanged(minPts);
}

void ClusteringWidget::updateEps(double eps)
{
	emit epsChanged(eps);
}

