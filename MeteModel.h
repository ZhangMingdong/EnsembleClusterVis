
#pragma once
#include <QGLWidget>
#include "ContourGenerator.h"
#include "ContourStripGenerator.h"
#include <MathTypes.hpp>
#include "def.h"

class DataField;
class OneSpatialCluster;

// record the cluster result
struct ClusterResult {
	int _nM = 0;								// dimension of each item. also length of _arrLabels;
	int _nK=0;									// number of clusters 1~5;
	std::vector<int> _vecItems[g_nClusterMax];	// items of each cluster, only the first _nK are used
	int _arrLabels[g_nEnsembles];				// labels of each ensemble member
	// push a label
	void PushLabel(int nIndex, int nLabel);
	// match this cluster index with c, minimized the items that changing cluster index
	void Match(ClusterResult& c);
	// align with a map
	void AlighWith(int* arrMap);
	// reset the result
	void Reset(int nM, int nK);
	// sort clusters by number of items in them
	void Sort();
	// generate items by the array of labels
	void generateItemsByLabels();

};


/*
	model of the meteorology data
	Mingdong
	2017/05/05
*/
class MeteModel
{
public:
	MeteModel();
	virtual ~MeteModel();
public:
	// initialize the model
	virtual void InitModel(int nEnsembleLen, int nWidth, int nHeight, int nFocusX, int nFocusY, int nFocusW, int nFocusH
		, QString strFile, bool bBinary = false
		, int nWest = -179, int nEast = 180, int nSouth = -90, int nNorth = 90
		, int nFocusWest = -179, int nFocusEast = 180, int nFocusSouth = -90, int nFocusNorth = 90, bool bFiltered = false);
	// generate color mapping texture
	virtual GLubyte* generateTextureNew();

	virtual int GetW(){ return _nWidth; }
	virtual int GetH(){ return _nHeight; }
	virtual int GetFocusW(){ return _nFocusW; }
	virtual int GetFocusH(){ return _nFocusH; }
	virtual int GetFocusX(){ return _nFocusX; }
	virtual int GetFocusY(){ return _nFocusY; }
	virtual int GetWest() { return _nWest; }
	virtual int GetEast() { return _nEast; }
	virtual int GetSouth() { return _nSouth; }
	virtual int GetNorth() { return _nNorth; }
	virtual int GetFocusWest() { return _nFocusWest; }
	virtual int GetFocusEast() { return _nFocusEast; }
	virtual int GetFocusSouth() { return _nFocusSouth; }
	virtual int GetFocusNorth() { return _nFocusNorth; }
	virtual bool GetFilter() {
		return _bFilter; 
	}
	virtual QList<ContourLine> GetContourMin(){ return _listContourMinE; }
	virtual QList<ContourLine> GetContourMax(){ return _listContourMaxE; }
	virtual QList<ContourLine> GetContourMean(){ return _listContourMeanE; }
	virtual QList<QList<ContourLine>> GetContour();
	virtual QList<QList<ContourLine>> GetContourBrushed();
	virtual QList<QList<ContourLine>> GetContourNotBrushed();
	virtual QList<UnCertaintyArea*> GetUncertaintyArea(){ return _listUnionAreaE; }
	virtual int GetUncertaintyAreas() { return _nUncertaintyRegions; }

	virtual const QList<QList<UnCertaintyArea*> > GetUncertaintyAreaG() {
		return _listUnionAreaEG;
	}
	// get a vector of the sorted variance;
	virtual std::vector<double> GetVariance();
	virtual void SetUncertaintyAreas(int nAreas);
protected:
	// specialized model initialization
	virtual void initializeModel();				
protected:
	// read ensemble data from text file
	virtual void readDataFromText();

	/*
		calculate the signed distance function
		arrData:	the input data
		arrSDF:		the calculated sdf
		isoValue:	the iso value
		contour:	the generated contour
	*/
	void calculateSDF(const double* arrData, double* arrSDF, int nW, int nH, double isoValue, QList<ContourLine> contour);
	
	// space segmentation
	void generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas);

	// read data from binary file
	void readData();
public:
	enum enumBackgroundFunction
	{
		bg_mean,			// mean
		bg_vari,			// variance
		bg_cluster,			// spatial cluster
		bg_sdf,				// signed distance function
		bg_vari_smooth		// smooth variance
	};
protected:
	// using which background function
	enumBackgroundFunction _bgFunction= bg_mean;
public:
	// set background function
	void SetBgFunction(enumBackgroundFunction f);
	// brushing
	virtual void Brush(int nLeft, int nRight, int nTop, int nBottom);

protected:	
	// 1.raw data
	DataField* _pData=0;						// the data		
	int _nEnsembleLen;						// number of ensemble members
	int _nWidth;
	int _nHeight;	
	int _nLen;								// _nWidth*_nHeight
	int _nFocusX;
	int _nFocusY;
	int _nFocusW;
	int _nFocusH;
	int _nFocusLen;							 //_nFocusW*_nFocusH
	int _nWest;
	int _nEast;
	int _nSouth;
	int _nNorth;
	int _nFocusWest;
	int _nFocusEast;
	int _nFocusSouth;
	int _nFocusNorth;

	// 2.basic contours and areas
	QList<UnCertaintyArea*> _listUnionAreaE;			// list of the uncertainty area of union of E
	QList<ContourLine> _listContourMinE;				// list of contours of minimum of E
	QList<ContourLine> _listContourMaxE;				// list of contours of maximum of E
	QList<ContourLine> _listContourMeanE;				// list of contours of mean of E
	QList<QList<ContourLine>> _listContour;				// list of contours of ensemble members
	QList<QList<ContourLine>> _listContourBrushed;		// list of brushed contours of ensemble members
	QList<QList<ContourLine>> _listContourNotBrushed;	// list of not brushed contours of ensemble members


	QList<QList<UnCertaintyArea*>> _listUnionAreaEG;	// list of the uncertainty area of union of E	(for gradient)
	QList<QList<ContourLine>> _listContourMinEG;
	QList<QList<ContourLine>> _listContourMaxEG;


	int _nUncertaintyRegions = 6;				// number of uncertainty regions
	int _nClusters = 5;							// number of clusters
	


	// 0.io related	
	QString _strFile;				// file name of the data	
	bool _bBinaryFile;				// whether read binary file	
	bool _bFilter;					// filtered the data between grids	

	double _dbVarThreshold = 1.5;	// threshold of the variance

	int _nMinPts = 1;
	double _dbEps = 10.0;

	int _nSmooth = 1;	// smooth level 1~5

	// data of the texture
	GLubyte* _dataTexture = NULL;


	// points generated by PCA
	std::vector<DPoint3> _vecPCAPoints[g_nUncertaintyAreaMax];

	// cluster result of the first two spatial clusters
	ClusterResult _arrClusterResult[g_nUncertaintyAreaMax];

	// matrix of the similarity between uncertainty regions
	int _matrixSimilarity[g_nUncertaintyAreaMax][g_nUncertaintyAreaMax];

	// matrix of the similarity between uncertainty regions of each cluster
	int _matrixSimilarityC[g_nUncertaintyAreaMax][g_nUncertaintyAreaMax][g_nClusterMax];

public:
	// interface of the creation of model
	static MeteModel* CreateModel();

	void SetVarThreshold(double dbThreshold);
	double GetVarThreshold() { return _dbVarThreshold; }

	void SetMinPts(int minPts) { _nMinPts = _nMinPts; }
	void SetEps(double dbEps) { _dbEps = dbEps; }

	void SetSmooth(int nSmooth);

	// get the cluster similarity between two uncertainty regions
	int GetRegionSimilarity(int nR1, int nR2) { return _matrixSimilarity[nR1][nR2]; }

	int GetRegionSimilarity(int nR1, int nR2,int nC) { return _matrixSimilarityC[nR1][nR2][nC]; }

private:
	// generate texture of clustered variance
	void buildTextureClusteredVariance();

	// generate texture use threshold var
	void buildTextureThresholdVariance();

	// generate texture of colormap of mean or variance
	void buildTextureColorMap();

	// generate texture of smoothed variance
	void buildTextureSmoothedVariance();

	// generate PCA points from a spatial cluster
	void generatePCAPoint(OneSpatialCluster& cluster, std::vector<DPoint3>& points);

	// cluster according to an area
	void clusterSpatialArea(OneSpatialCluster& cluster, std::vector<DPoint3>& points, ClusterResult&cr);

	// reset labels for pca points
	void setLabelsForPCAPoints(std::vector<DPoint3>& points, ClusterResult&cr);

	// generate new texture
	void regenerateTexture();
public:
	// get the pca points of each uncertainty area
	std::vector<DPoint3>* GetPCAPoints() {
		return _vecPCAPoints;
	}
	// get the cluster result of each uncertainty area
	const ClusterResult* GetClusterResults() { return _arrClusterResult; }	
};

