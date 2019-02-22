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
	QSpinBox *_pSpinBoxFocusedRegion;
	QSpinBox *_pSpinBoxEOF;
	QSpinBox *_pSpinBoxMember;
	QSpinBox *_pSpinBoxEnsCluster;
	QSpinBox *_pSpinBoxEnsClusterLen;

	/*
		displayed level of contour
		0-1;
		1-3;
		2-7;
		3-15;
		4-31;
		5-50;
	*/
	QSpinBox *_pSpinBoxContourLevel;

	/*
		time step
		0-360 by 6
	*/
	QSpinBox *_pSpinBoxTimeStep;

private:
	void createWidgets();
	void createLayout();
	void createConnections();

private slots:
	void updateBgFunction(int nFunction);
	void updateSmooth(int nSmooth);
	void updateUncertaintyAreas(int nArea);
	void updateFocusedCluster(int nFocusedCluster);
	void updateEOF(int nEOF);
	void updateMember(int nMember);
	void updateEnsCluster(int nEnsCluster);
	void updateEnsClusterLen(int nEnsCluster);
	void updateContourLevel(int nLevel);
	void updateTimeStep(int nTS);
signals:
	void bgFunctionChanged(int nFunction);
	void smoothChanged(int nSmooth);
	void areasChanged(int nArea);
	void focusedClusterChanged(int nFocusedCluster);
	void EOFChanged(int nEOF);
	void MemberChanged(int nMember);
	void EnsClusterChanged(int nCluster);
	void EnsClusterLenChanged(int nClusterLen);
	void ContourLevelChanged(int nLevel);
	void TimeStepChanged(int nTS);
};
