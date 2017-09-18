#pragma once

#include <QWidget>
class QComboBox;
class QSpinBox;

class DisplayCtrlWidget : public QWidget
{
	Q_OBJECT

public:
	DisplayCtrlWidget(QWidget *parent=NULL);
	~DisplayCtrlWidget();
private:
	// clustering method
	QComboBox *_pCombBgFunction;
	QSpinBox *_pSpinBoxSmooth;
	QSpinBox *_pSpinBoxAreas;
	QSpinBox *_pSpinBoxFocusedCluster;

private:
	void createWidgets();
	void createLayout();
	void createConnections();

private slots:
	void updateBgFunction(int nFunction);
	void updateSmooth(int nSmooth);
	void updateUncertaintyAreas(int nArea);
	void updateFocusedCluster(int nFocusedCluster);
signals:
	void bgFunctionChanged(int nFunction);
	void smoothChanged(int nSmooth);
	void areasChanged(int nArea);
	void focusedClusterChanged(int nFocusedCluster);
};
