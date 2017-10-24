#pragma once
#include <QGLWidget>
#include "ContourGenerator.h"
#include "ContourStripGenerator.h"
#include <MathTypes.hpp>
#include "def.h"
#include "ClusterResult.h"
#include "SpatialCluster.h"

class DataField;

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
		, int nWest = -180, int nEast = 180, int nSouth = -90, int nNorth = 90
		, int nFocusWest = -180, int nFocusEast = 180, int nFocusSouth = -90, int nFocusNorth = 90);
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
	// read the dip value
	virtual void readDipValue(char* strFileName);
	virtual void readDipValueG(char* strFileName);

	// read data from text file for global area, set the left line according to the right line.
	virtual void readDataFromTextG();

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
		bg_mean,				// mean
		bg_vari,				// variance
		bg_cluster,				// spatial cluster
		bg_varThreshold,		// thresholded variance
		bg_vari_smooth,			// smooth variance
		bg_dipValue,			// variance of variance
		bg_dipValueThreshold,	// thresholded dip value
		bg_EOF					// EOF
	};
protected:
	// using which background function
	enumBackgroundFunction _bgFunction= bg_mean;
public:
	// set background function
	void SetBgFunction(enumBackgroundFunction f);
	// get the background function
	enumBackgroundFunction GetBgFunction() { return _bgFunction; }
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

	double _dbVarThreshold = 1.5;	// threshold of the variance



	// smooth level 1~5
	int _nSmooth = 1;	
	// EOF: 1~5
	int _nEOF = 1;

	// data of the texture
	GLubyte* _dataTexture = NULL;


	// points generated by PCA
	std::vector<DPoint3> _vecPCAPoints[g_nUncertaintyAreaMax];

	// cluster result matrix
	ClusterResult _mxClusterResult[g_nUncertaintyAreaMax][g_nUncertaintyAreaMax];

	// matrix of the similarity between uncertainty regions of each cluster
	int _mxSimilarity[g_nUncertaintyAreaMax][g_nUncertaintyAreaMax][g_nClusterMax];

	// grid points of thresholded variance
	int _nThresholdedGridPoints = 0;

	// vector of uncertainty regions
	std::vector<UncertaintyRegion> _vecRegions;

	// index of the focused region
	int _nFocusedRegion = 0;

private:
	// generate texture of clustered variance
	void buildTextureClusteredVariance();

	// generate texture use threshold var
	void buildTextureThresholdVariance();


	// generate texture use thresholded dip value
	void buildTextureThresholdDipValue();

	// generate texture of colormap of mean or variance
	void buildTextureColorMap();
	
	// generate PCA points from a spatial cluster
	void generatePCAPoint(UncertaintyRegion& cluster, std::vector<DPoint3>& points);

	// cluster according to an area
	void clusterSpatialArea(UncertaintyRegion& cluster, std::vector<DPoint3>& points, ClusterResult&cr);

	// reset labels for pca points
	void setLabelsForPCAPoints(std::vector<DPoint3>& points, ClusterResult&cr);

	// generate new texture
	void regenerateTexture();

	// cluster in each region
	void regionCluster();

	// calculate similarities between regions
	void calculateSimilarity();

	// align the cluster results
	void alignClusters();
public:
	// interface of the creation of model
	static MeteModel* CreateModel();
public:
	// wrappers
	void SetVarThreshold(double dbThreshold);
	double GetVarThreshold() { return _dbVarThreshold; }



	void SetSmooth(int nSmooth);
	void SetEOF(int nEOF);

	// get the cluster similarity between two uncertainty regions

	int GetRegionSimilarity(int nR1, int nR2, int nC) { return _mxSimilarity[nR1][nR2][nC]; }

	// get the pca points of each uncertainty area
	std::vector<DPoint3>* GetPCAPoints() {
		return _vecPCAPoints;
	}
	// get the cluster result of each uncertainty area
	const ClusterResult* GetClusterResults() { return _mxClusterResult[0]; }	


	const ClusterResult* GetClusterResultOfFocusedRegion() { return _mxClusterResult[0]+_nFocusedRegion; }

	int GetFocusedRegion() { return _nFocusedRegion; }


	int GetGridLength() { return _nLen; }
	int GetThresholdedGridLength() { return _nThresholdedGridPoints; }

	const std::vector<UncertaintyRegion> GetUncertaintyRegions() { return _vecRegions; }

	void SetFocusedRegion(int nRegion);
};

