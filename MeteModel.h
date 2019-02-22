#pragma once
#include <QGLWidget>
#include "ContourStripGenerator.h"
#include <MathTypes.hpp>
#include "def.h"
#include "SpatialCluster.h"

#include "GridFrame.h"

#include "TimeStep.h"

class DataField;
class FeatureSet;

class TimeStep;

/*
	model of the meteorology data
	Mingdong
	2017/05/05
*/
class MeteModel:public GridFrame
{
	Q_OBJECT
public:
	MeteModel();
	virtual ~MeteModel();
protected:
	enum ModelType {
		MT_T2				// T2 data using t2 color map
		, MT_AF				// artificial model, using colormap centered at 0
	};
	ModelType _type = MT_T2;

public:
	// initialize the model
	virtual void InitModel(int nEnsembleLen, int nWidth, int nHeight
		, QString strFile, bool bBinary = false
		, int nWest = -180, int nEast = 180, int nSouth = -90, int nNorth = 90);
	// generate color mapping texture
	virtual GLubyte* GenerateTexture();


	virtual DataField* GetData() { return _pTimeStep->_pData; }

	virtual QList<ContourLine> GetContourMin(int isoIndex = 0);
	virtual QList<ContourLine> GetContourMax(int isoIndex = 0);
	virtual QList<ContourLine> GetContourMean(int isoIndex = 0);
	virtual QList<ContourLine> GetContourMedian(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContour(int isoIndex=0);
	virtual QList<QList<ContourLine>> GetContourSmooth(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContourOutlier(int isoIndex =0);
	virtual QList<QList<ContourLine>> GetContourSorted(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContourSDF(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContourSortedSDF(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContourResampled(int isoIndex = 0);
	virtual QList<QList<ContourLine>> GetContourBrushed();
	virtual QList<QList<ContourLine>> GetContourNotBrushed();
	virtual QList<UnCertaintyArea*> GetUncertaintyArea(int isoIndex = 0);
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaValid(int isoIndex = 0);
	virtual QList<UnCertaintyArea*> GetUncertaintyAreaHalf(int isoIndex = 0) ;
	virtual int GetUncertaintyAreas() { return _nUncertaintyRegions; }
	virtual int GetLabel(int l);
	virtual void GetMerge(int l, int& nSource, int& nTarget);
	virtual double GetPC(int l, int nIndex);
	virtual int GetClusters();

	virtual const QList<QList<UnCertaintyArea*> > GetUncertaintyAreaG() {
		return _listUnionAreaEG;
	}
	// get a vector of the sorted variance;
	virtual std::vector<double> GetVariance();
	virtual void SetUncertaintyAreas(int nAreas);
protected:
	// specialized model initialization
	virtual void initializeModel();		
	void updateTimeStep();				// update according to current time step
protected:
	// read ensemble data from text file
	virtual void readDataFromText(QString filename);
	
	// read the dip value
	virtual void readDipValue(char* strFileName);

	virtual void readDipValueG(char* strFileName);

	// read data from text file for global area, set the left line according to the right line.
	virtual void readDataFromTextG();


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
		bg_EOF,					// EOF
		bg_Obs,					// obs
		bg_err,					// error
		bg_SDF,					// SDF
		bg_LineKernel,			// SDF
		bg_LineKernelX,			// SDF
		bg_LineKernelY,			// SDF
		bg_LineKernelZ,			// SDF
		bg_ICDVX,				// SDF
		bg_ICDVY,				// SDF
		bg_ICDVW,				// SDF
		bg_IsoContourDensity	// density of the iso-contour
	};
protected:
	// using which background function
	enumBackgroundFunction _bgFunction= bg_mean;

protected:	
	// 0.io related	
	QString _strFile;				// file name of the data	
	bool _bBinaryFile;				// whether read binary file	

	TimeStep* _pTimeStep;			// current time step

	TimeStep* _arrTimeSteps[61];	// time step list

	double* _bufObs = 0;					// observation data


	// 2.basic contours and areas
	QList<QList<ContourLine>> _listContourBrushed;		// list of brushed contours of ensemble members
	QList<QList<ContourLine>> _listContourNotBrushed;	// list of not brushed contours of ensemble members
	QList<QList<ContourLine>> _listContourEOF;			// list of contours of EOF
	QList<QList<ContourLine>> _listMemberContour[g_nEnsembles];			// list of contours of each ensemble member
	QList<QList<ContourLine>> _listEnsClusterContour[g_nEnsClusterLen];			// list of contours of each ensemble cluster
	QList<QList<UnCertaintyArea*>> _listUnionAreaEG;	// list of the uncertainty area of union of E	(for gradient)
	QList<QList<ContourLine>> _listContourMinEG;
	QList<QList<ContourLine>> _listContourMaxEG;

	QList<double> _listIsoValues;						// list of iso values

	// 3.cluster related
	int _nUncertaintyRegions = 6;				// number of uncertainty regions
	double _dbVarThreshold = 1.5;				// threshold of the variance	
	int _nSmooth = 1;							// smooth level 1~5	
	int _nEOF = 1;								// EOF: 1~5
	int _nMember = 0;							// current focused member
	int _nEnsCluster = 0;						// current ensemble cluster
	
	GLubyte* _dataTexture = NULL;				// data of the texture
	int _nTexW;
	int _nTexH;
	
	std::vector<DPoint3> _vecPCAPoints[g_nUncertaintyAreaMax];	// points generated by PCA

	// matrix of the similarity between uncertainty regions of each cluster
	int _mxSimilarity[g_nUncertaintyAreaMax][g_nUncertaintyAreaMax][g_nClusterMax];

	// grid points of thresholded variance
	int _nThresholdedGridPoints = 0;

	// vector of uncertainty regions
	std::vector<UncertaintyRegion> _vecRegions;
	
	// 4.cluster the ensemble member directly
	QList<DataField*> _arrEnsClusterData;	// data of each cluster
	int _arrLabels[g_nEnsembles];	// labels of each member

	double* _gridErr;				// error field
private:
	// generate texture of clustered variance
	void buildTextureClusteredVariance();

	// generate regions for clustering in each region
	void generateRegions();

	// generate texture use threshold var
	void buildTextureThresholdVariance();


	// generate texture use thresholded dip value
	void buildTextureThresholdDipValue();

	// generate texture of colormap of mean or variance
	void buildTextureColorMap();


	// generate texture of signed distance function
	void buildTextureSDF();
	void buildTextureICDVX();
	void buildTextureICDVY();
	void buildTextureICDVW();

	void buildTextureICD_LineKernel();
	void buildTextureICDX();
	void buildTextureICDY();
	void buildTextureICDZ();

	// generate texture of iso-contour density
	void buildTextureICD();
	
	// generate PCA points from a spatial cluster
	void generatePCAPoint(UncertaintyRegion& cluster, std::vector<DPoint3>& points);

	// generate new texture
	void regenerateTexture();

	/*
		ensemble clustering
		2017/11/07
	*/
	void doEnsCluster();


	// read observation data
	void readObsData();
public:
	// interface of the creation of model
	static MeteModel* CreateModel(bool bA=false);
public:
	// wrappers

	// set background function
	void SetBgFunction(enumBackgroundFunction f);
	// get the background function
	enumBackgroundFunction GetBgFunction() { return _bgFunction; }
	// brushing
	virtual void Brush(int nLeft, int nRight, int nTop, int nBottom);

	void SetVarThreshold(double dbThreshold);
	double GetVarThreshold() { return _dbVarThreshold; }



	void SetSmooth(int nSmooth);
	void SetEOF(int nEOF);
	void SetMember(int nMember);
	void SetEnsCluster(int nEnsCluster);
	void SetEnsClusterLen(int nEnsClusterLen);
	void SetContourLevel(int nLevel);

	// get the cluster similarity between two uncertainty regions

	int GetRegionSimilarity(int nR1, int nR2, int nC) { return _mxSimilarity[nR1][nR2][nC]; }

	// get the pca points of each uncertainty area
	std::vector<DPoint3>* GetPCAPoints() {
		return _vecPCAPoints;
	}

	int GetGridLength() { return _nGrids; }
	int GetThresholdedGridLength() { return _nThresholdedGridPoints; }

	const std::vector<UncertaintyRegion> GetUncertaintyRegions() { return _vecRegions; }


	void SetIsoValues(QList<double> listIsoValues);
	QList<double> GetIsoValues() { return _listIsoValues; }


	int GetTexW() { return _nTexW; }
	int GetTexH() { return _nTexH; }

private:
	// detail level of the contours, 0-1;1-3;2-7...
	int _nContourLevel = 0;

public slots:


	void updateTimeStep(int nTS);
protected:
	int _nTime = 0;		// time step
signals:
	void UpdateView();
};

