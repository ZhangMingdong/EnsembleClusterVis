#pragma once

#include <QWidget>

class QDoubleSpinBox;
class QSpinBox;
class QComboBox;
/*
	widget to set the parameter of the clustering
*/
class ClusteringWidget : public QWidget
{
	Q_OBJECT

public:
	ClusteringWidget(QWidget *parent=0);
	~ClusteringWidget();

private:
	// clustering method
	QComboBox *_pCombMethod;
	// spinbox of the eps
	QDoubleSpinBox *_pSpinBoxEps;
	QSpinBox *_pSpinBoxMinPts;

private:
	void createWidgets();
	void createLayout();
	void createConnections();

private slots:
	void updateMethod(int method);
	void updateMinPts(int minPts);
	void updateEps(double eps);
signals:
	void methodChanged(int minPts);
	void minPtsChanged(int minPts);
	void epsChanged(double eps);
};
