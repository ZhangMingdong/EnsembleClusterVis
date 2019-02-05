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
	FeatureSet(DataField* pData,double _dbIsoValue,int nWidth,int nHeight,int _nEnsembleLen);
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
	double* _pSortedSDF;			// data of signed distance function
	bool * _pSet;					// set state of the grid point given iso-value
	bool* _pGridDiverse;			// diversed grid points
	int *_pSetBandDepth;			// sBandDepth
	int *_pMemberType;				// type of each member.0:outlier,1-100%,2-50%.
	int _nOutlierThreshold = 1;		// threshold for outliars
	int _nDiverseCount = 0;			// count of diverse grids
	double* _pResampledSDF;			// resampled SDF using kernel density function
	int _nResampleLen = 63;			// length of resample 2^7

	int _nPCLen=10;					// length of PC
	double* _arrPC;					// array of PCs

	int _nClusters = 5;			// number of clusters
	int* _arrLabels;				// labels of each member

	QList<UnCertaintyArea*> _listAreaValid;			// list of the uncertainty area of union of valid members
	QList<UnCertaintyArea*> _listAreaHalf;			// list of the uncertainty area of union of half valid members

	QList<UnCertaintyArea*> _listUnionAreaE;			// list of the uncertainty area of union of E
	QList<QList<ContourLine>> _listContourSorted;		// list of contours of sorted ensemble members
	QList<QList<ContourLine>> _listContourSDF;			// list of contours generated form SDF
	QList<QList<ContourLine>> _listContourSortedSDF;	// list of contours generated form sorted SDF
	QList<QList<ContourLine>> _listContourResampled;	// list of resampled contours
	QList<QList<ContourLine>> _listContourSmooth;		// list of smoothed contours

	int _nDetailScale = 100;								// scale of detail texture



public:
	QList<QList<ContourLine>>& GetContours() { return _listContour; }
	QList<QList<ContourLine>>& GetContoursSmooth() { return _listContourSmooth; }
	QList<ContourLine>& GetContourMin() { return _listContourMinE; }
	QList<ContourLine>& GetContourMax() { return _listContourMaxE; }
	QList<ContourLine>& GetContourMean(){ return _listContourMeanE; }
	QList<ContourLine>& GetContourMedian() { return _listContourMedianE; }
	QList<QList<ContourLine>>& GetContourSDF() { return _listContourSDF; }
	QList<QList<ContourLine>>& GetContourSorted() { return _listContourSorted; }
	QList<QList<ContourLine>>& GetContourSortedSDF() { return _listContourSortedSDF; }
	QList<QList<ContourLine>>& GetContourResampled() { return _listContourResampled; }

	const double* GetValidMax() { return _gridValidMax; };
	const double* GetValidMin() { return _gridValidMin; };
	const double* GetHalfMax() { return _gridHalfMax; };
	const double* GetHalfMin() { return _gridHalfMin; };
	const double* GetMedian();
	const double* GetSDF() { return _pSDF; }


	const double* GetSortedSDF(int l) { return _pSortedSDF + l * _nGrids; }
	bool* GetSet(int l) { return _pSet + l * _nGrids; }
	int GetMemberType(int l) { return _pMemberType[l]; }
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaValid() { return _listAreaValid; }
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaHalf() { return _listAreaHalf; }
	virtual QList<UnCertaintyArea*> GetUncertaintyArea() { return _listUnionAreaE; }
	int GetDetailScale() { return _nDetailScale; }
	int GetLabel(int l) { return _arrLabels[l]; }
	double GetPC(int l,int nIndex) { return _arrPC[l * _nPCLen+nIndex]; }


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
	void smoothContours();
	void buildSortedSDF();
	void resampleContours();
	void calculateSimilarityMatrix();		// calculate similarity matrix of the ensemble members
	/*
		calculate similarity between two members;
		paras:
			l1,l2: index of the two member
			min,max: min and max of the sdf
	*/
	double calculateSimilarity(int l1, int l2,double dbMin, double dbMax);	
public:
	static void CalculateSDF(double dbIsoValue, int nEnsembleLen, int nWidth, int nHeight, DataField* pData, double* pSDF
		, QList<QList<ContourLine>> listContour);

// dimension-reduction and clustering
private:
	void calculateDiverse();			// calculate diverse grids and diverse count
	void calculatePCA();				// calculate pca
	void calculatePCA_MDS();			// calculate pca using MDS
	void calculatePCA_MDS_Dis();		// calculate pca using MDS by distance
	void calculatePCA_MDS_Whole();		// calculate pca using MDS of the whole grids
	void calculatePCA_MDS_Whole_Density();		// calculate pca using MDS of the whole grids by density
	void calculatePCA_MutualInformation();		// calculate pca using mutual information
	void doClustering();				// clustering
	void doPCAClustering();				// clustering based on pca result
	double* createDiverseArray();		// create diverse grids array, the result should be manually freed.
	void calculatePCARecovery();		// test, recovering,calculate the box of pca
	void calculatePCABox();				// test, recovering,calculate the box of pca


//=============================density================================
private:
	double* _pICD;					// data of iso-contour density, using kde

	double* _pICDVX;				// data of iso-contour density vector
	double* _pICDVY;				// data of iso-contour density vector
	double* _pICDVW;				// data of iso-contour density vector
	double* _pICD_LineKernel;		// data of iso-contour density using line kernel
	double* _pICDX;					// data of iso-contour density using line kernel X
	double* _pICDY;					// data of iso-contour density using line kernel Y
	double* _pICDZ;					// data of iso-contour density using line kernel Z

private:
	void calculateICD();			// calculate _pICD

	void buildICD_LineKernel();		// build iso-contour density using line kernel
	void buildICD_Vector();			// build iso-contour density using vector kernel
	void buildICDV();
public:
	const double* GetICD() { return _pICD; }

	const double* GetICD_LineKernel() { return _pICD_LineKernel; }
	const double* GetICDX() { return _pICDX; }
	const double* GetICDY() { return _pICDY; }
	const double* GetICDZ() { return _pICDZ; }
	const double* GetSDF(int l) const { return _pSDF + l * _nGrids; }
	const double* GetICDVX(int l) const { return _pICDVX + l * _nGrids; }
	const double* GetICDVY(int l) const { return _pICDVY + l * _nGrids; }
	const double* GetICDVW(int l) const { return _pICDVW + l * _nGrids; }
};

