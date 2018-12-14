#pragma once
#include <QList>
#include "BasicStruct.h"
#include "GridFrame.h"

class DataField;
class UnCertaintyArea;
/*
	This class used to capsulate the features generated by an isovalue
	Mingdong
	2018/12/04
*/
class FeatureSet:public GridFrame
{
public:
	FeatureSet(DataField* pData,double _dbIsoValue,int nWidth,int nHeight,int _nEnsembleLen,int nFocusX,int nFocusY,int nFocusW,int nFocusH);
	~FeatureSet();
private:
	DataField* _pData;									// reference to the data set
	QList<QList<ContourLine>> _listContour;				// list of contours of ensemble members
	QList<ContourLine> _listContourMinE;				// list of contours of minimum of E
	QList<ContourLine> _listContourMaxE;				// list of contours of maximum of E
	QList<ContourLine> _listContourMeanE;				// list of contours of mean of E
	double _dbIsoValue = 0;								// isovalue of this feature set;

	QList<ContourLine> _listContourMedianE;				// list of contours of median of E
	QList<ContourLine> _listContourMinValid;			// list of contours of minimum of valid members
	QList<ContourLine> _listContourMaxValid;			// list of contours of maximum of valid members
	QList<ContourLine> _listContourMinHalf;				// list of contours of minimum of half valid members
	QList<ContourLine> _listContourMaxHalf;				// list of contours of maximum of half valid members

	double* _gridHalfMax;			// maximum of half valid ensemble results
	double* _gridHalfMin;			// minimum of half valid ensemble results
	double* _gridValidMax;			// maximum of valid ensemble results
	double* _gridValidMin;			// minimum of valid ensemble results
	int _nMedianIndex = -1;			// index of median of the contour
	double* _pSDF;					// data of signed distance function
	double* _pICD;					// data of iso-contour density
	double* _pSortedSDF;			// data of signed distance function
	bool * _pSet;					// set state of the grid point given iso-value
	bool* _pGridDiverse;			// diversed grid points
	int *_pSetBandDepth;			// sBandDepth
	int *_pMemberType;				// type of each member.0:outlier,1-100%,2-50%.
	int _nOutlierThreshold = 1;		// threshold for outliars
	int _nDiverseCount = 0;			// count of diverse grids

	QList<UnCertaintyArea*> _listAreaValid;			// list of the uncertainty area of union of valid members
	QList<UnCertaintyArea*> _listAreaHalf;			// list of the uncertainty area of union of half valid members

	QList<UnCertaintyArea*> _listUnionAreaE;			// list of the uncertainty area of union of E
	QList<QList<ContourLine>> _listContourSorted;		// list of contours of sorted ensemble members
	QList<QList<ContourLine>> _listContourSDF;			// list of contours generated form SDF
	QList<QList<ContourLine>> _listContourSortedSDF;	// list of contours generated form sorted SDF
public:
	QList<QList<ContourLine>>& GetContours() { return _listContour; }
	QList<ContourLine>& GetContourMin() { return _listContourMinE; }
	QList<ContourLine>& GetContourMax() { return _listContourMaxE; }
	QList<ContourLine>& GetContourMean(){ return _listContourMeanE; }
	QList<ContourLine>& GetContourMedian() { return _listContourMedianE; }
	QList<QList<ContourLine>>& GetContourSDF() { return _listContourSDF; }
	QList<QList<ContourLine>>& GetContourSorted() { return _listContourSorted; }
	QList<QList<ContourLine>>& GetContourSortedSDF() { return _listContourSortedSDF; }

	const double* GetValidMax() { return _gridValidMax; };
	const double* GetValidMin() { return _gridValidMin; };
	const double* GetHalfMax() { return _gridHalfMax; };
	const double* GetHalfMin() { return _gridHalfMin; };
	const double* GetMedian();
	const double* GetSDF() { return _pSDF; }
	const double* GetSDF(int l) const { return _pSDF + l * _nGrids; }
	const double* GetSortedSDF(int l) { return _pSortedSDF + l * _nGrids; }
	bool* GetSet(int l) { return _pSet + l * _nGrids; }
	int GetMemberType(int l) { return _pMemberType[l]; }
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaValid() { return _listAreaValid; }
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaHalf() { return _listAreaHalf; }
	virtual QList<UnCertaintyArea*> GetUncertaintyArea() { return _listUnionAreaE; }
	double* GetICD() { return _pICD; }


private:
	void sortBuf(const double* pS, double* pD);	// space segmentation
	void generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas);

	double* getSDF(int l) { return _pSDF + l * _nGrids; }
// basic operation
private:
	void calculateSet();
	void calculateBandDepth();		// new version, use the method in contour boxplot
	void calculateBandDepth_v1();	// old version, use a threshold
	void calculateMemberType();
	void doStatistics();			// calculate meadian, valid, and 50%
	void generateContours();
	void buildSortedSDF();
	void calculateICD();
public:
	void CalculateSDF(double dbIsovalue);	// used in constructor and Artificial Model

// dimension-reduction and clustering
private:
	void calculateDiverse();			// calculate diverse grids and diverse count
	void calculatePCA();				// calculate pca
	void doClustering();				// clustering
	void doPCAClustering();				// clustering based on pca result
	double* createDiverseArray();		// create diverse grids array, the result should be manually freed.
	void calculatePCARecovery();		// test, recovering,calculate the box of pca
	void calculatePCABox();				// test, recovering,calculate the box of pca
};

