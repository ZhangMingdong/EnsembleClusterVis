#include "MeteModel.h"

#include <iostream>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QMessageBox>

#include "ContourGenerator.h"
#include "ColorMap.h"

#include "AHCClustering.h"
#include "KMeansClustering.h"
#include "MyPCA.h"
#include "DataField.h"

#include <fstream>
#include <iomanip>
#include <assert.h>
#include <sstream>

#include "def.h"

#include "DBSCANClustering.h"


#include "EnsembleModel.h"
#include "ReanalysisModel.h"
#include "ContourBandDepthModel.h"

using namespace std;



/*
double PointToSegDist(double x, double y, double x1, double y1, double x2, double y2)
{
	double cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
	if (cross <= 0) return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));

	double d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
	if (cross >= d2) return sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));

	double r = cross / d2;
	double px = x1 + (x2 - x1) * r;
	double py = y1 + (y2 - y1) * r;
	return sqrt((x - px) * (x - px) + (py - y) * (py - y));
}
*/
double PointToSegDist(double x, double y, double x1, double y1, double x2, double y2)
{
	double cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
	if (cross <= 0)
		return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));

	double d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
	if (cross >= d2)
		return sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));

	double r = cross / d2;
	double px = x1 + (x2 - x1) * r;
	double py = y1 + (y2 - y1) * r;
	return sqrt((x - px) * (x - px) + (y - py) * (y - py));
}

MeteModel::MeteModel() :_pData(0)
, _pSDF(0)
, _arrLabels(0)
, _arrGridLabels(0)
, _pClusterCenter(0)
, _bgFunction(bg_mean)
{
}

MeteModel::~MeteModel()
{
	if (_pData)
	{
		delete _pData;
	}
	if (_pSDF) {
		delete _pSDF;
	}

	for each (UnCertaintyArea* pArea in _listUnionAreaE)
		delete pArea;

	for each (DataField* pData in _listClusterData)
	{
		delete pData;
	}

	for (size_t i = 0; i < _listUnionAreaEC.length(); i++)
	{
		for each (UnCertaintyArea* pArea in _listUnionAreaEC[i])
			delete pArea;
	}
	if (_arrLabels)
	{
		delete[] _arrLabels;
	}
	if (_arrGridLabels)
	{
		delete[] _arrGridLabels;
	}
	if (_pClusterCenter)
	{
		delete[] _pClusterCenter;
	}

	for (size_t i = 0; i < _listUnionAreaEG.length(); i++)
	{
		for each (UnCertaintyArea* pArea in _listUnionAreaEG[i])
			delete pArea;
	}
}

void MeteModel::InitModel(int nEnsembleLen, int nWidth, int nHeight, int nFocusX, int nFocusY, int nFocusW, int nFocusH
	, QString strFile, bool bBinary, int nWest, int nEast, int nSouth, int nNorth
	, int nFocusWest, int nFocusEast, int nFocusSouth, int nFocusNorth,bool bFilter) {
	// 0.record states variables
	_nClusters = g_nClusters;

	_nEnsembleLen = nEnsembleLen;
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nLen = nHeight*nWidth;

	_nFocusX = nFocusX;
	_nFocusY = nFocusY;
	_nFocusW = nFocusW;
	_nFocusH = nFocusH;

	_nFocusLen = _nFocusW*_nFocusH;

	_nWest = nWest;
	_nEast = nEast;
	_nSouth = nSouth;
	_nNorth = nNorth;
	_nFocusWest = nFocusWest;
	_nFocusEast = nFocusEast;
	_nFocusSouth = nFocusSouth;
	_nFocusNorth = nFocusNorth;

	_strFile = strFile;
	_bBinaryFile = bBinary;
	_bFilter = bFilter;

	// 2.allocate resource
	_pData = new DataField(_nWidth, _nHeight, _nEnsembleLen);
	_pSDF = new DataField(_nWidth, _nHeight, _nEnsembleLen);
	_arrLabels = new int[_nEnsembleLen];
	_arrGridLabels = new int[_nLen];

	// 3.read data
	if (_bBinaryFile)
	{
		readData();
	}
	else {
		readDataFromText();
	}

	// 4.statistic
	_pData->DoStatistic();

	if (g_bSpatialClustering) {
		for (size_t i = 0; i < _nLen; i++)
		{
			_arrGridLabels[i] = -3;
		}
		// DBSCAN
		doSpatialClustering();
	}

	// 5.generate features
	ContourGenerator generator;
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		QList<ContourLine> contour;
		generator.Generate(_pData->GetLayer(i), contour, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		_listContour.push_back(contour);
	}
	if (g_bSDFBand||g_bClustering)
	{
		// calculate SDF
		for (size_t i = 0; i < _nEnsembleLen; i++)
		{
			calculateSDF(_pData->GetLayer(i), _pSDF->GetEditableLayer(i), _nWidth, _nHeight, g_fThreshold, _listContour[i]);
		}
		_pSDF->DoStatistic();
	}


	if (g_bSDFBand)
	{
		generator.Generate(_pSDF->GetUMin(), _listContourMinE, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pSDF->GetUMax(), _listContourMaxE, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pSDF->GetMean(), _listContourMeanE, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		// generate gradient feature
		double* arrTempMin = new double[_nLen];
		double* arrTempMax = new double[_nLen];
		const double* pMean = _pSDF->GetMean();
		const double* pVari = _pSDF->GetVari();
		double dbScale = 0.1;
		for (size_t i = 0; i < g_gradient_l; i++)
		{
			for (size_t j = 0; j < _nLen; j++)
			{
				arrTempMin[j] = pMean[j] - dbScale*pVari[j];
				arrTempMax[j] = pMean[j] + dbScale*pVari[j];
			}
			dbScale += 0.1;
			QList<ContourLine> contourMin;
			QList<ContourLine> contourMax;
			generator.Generate(arrTempMin, contourMin, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			generator.Generate(arrTempMax, contourMax, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourMinEG.push_back(contourMin);
			_listContourMaxEG.push_back(contourMax);
		}
		delete[] arrTempMin;
		delete[] arrTempMax;

		for (size_t i = 0; i < g_gradient_l; i++)
		{
			QList<UnCertaintyArea*> area;
			generateContourImp(_listContourMinEG[i], _listContourMaxEG[i], area);
			_listUnionAreaEG.push_back(area);
		}
	}
	else {
		generator.Generate(_pData->GetUMin(), _listContourMinE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pData->GetUMax(), _listContourMaxE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pData->GetMean(), _listContourMeanE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	}

//	generateContourImp(_listContourMinE, _listContourMaxE, _listUnionAreaE);

	if (g_bClustering)
	{
		// 1.PCA
		doPCA();

		bool bUsingGlobalData = false;	// not just focus the iso value
		double dbIsoValue;
		// 2.generate clustered data
		if (bUsingGlobalData)
		{
			_pData->GenerateClusteredData(_listClusterLen, _arrLabels, _listClusterData);
			dbIsoValue = g_fThreshold;
		}
		else {
			dbIsoValue = 0;
			_pSDF->GenerateClusteredData(_listClusterLen, _arrLabels, _listClusterData);
		}
		// 3.generate features for clusters
		for (size_t i = 0; i < _nClusters; i++)
		{
			QList<ContourLine> listContourMin;
			QList<ContourLine> listContourMax;
			QList<ContourLine> listContourMean;
			QList<QList<ContourLine>> listContour;
			QList<UnCertaintyArea*> listUnionArea;
			generator.Generate(_listClusterData[i]->GetMean(), listContourMean, dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			generator.Generate(_listClusterData[i]->GetUMin(), listContourMin, dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			generator.Generate(_listClusterData[i]->GetUMax(), listContourMax, dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourMinEC.push_back(listContourMin);
			_listContourMaxEC.push_back(listContourMax);
			_listContourMeanEC.push_back(listContourMean);
			_listContourC.push_back(listContour);
			_listUnionAreaEC.push_back(listUnionArea);

			// pca
			QList<ContourLine> listContourMeanPCA;
			generator.Generate(_pClusterCenter + i*_nLen, listContourMeanPCA, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourMeanPCA.push_back(listContourMeanPCA);
		}
		// 4.put the spaghetti into different list according to their clusters
		for (size_t i = 0; i < _nEnsembleLen; i++)
		{
			_listContourC[_arrLabels[i]].push_back(_listContour[i]);
		}
		// 5.generate uncertainty areas
		for (size_t i = 0; i < _nClusters; i++)
		{
			if (_listClusterLen[i]>1)
			{
				generateContourImp(_listContourMinEC[i], _listContourMaxEC[i], _listUnionAreaEC[i]);
			}
		}
	}
	// specializaed initialization
	initializeModel();
}

void MeteModel::readDataFromText() {

	int nTimeStep = g_nTimeStep;

	QFile file(_strFile);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	int nCount = 0;
	QTextStream in(&file);
	int tt = 0;
	while (!in.atEnd()) {
		// every time step
		for (size_t t = 0; t <= nTimeStep; t++)
		{
			// every ensemble member
			for (int i = 0; i < _nEnsembleLen; i++)
			{
				QString line = in.readLine();
				// every grid
				for (int j = 0; j < _nLen; j++)
				{
					QString line = in.readLine();
					nCount++;
					if (t == nTimeStep)
						_pData->SetData(i, j, line.toFloat());
					int r = j / _nWidth;
					int c = j%_nWidth;
					
					if (_bFilter&&c < _nWidth - 1) {

						in.readLine();
						nCount++;
					}
					else if (_bFilter&&r < _nHeight - 1) {
						for (size_t ii = 0, length = 2 * (_nWidth - 1) + 1; ii < length; ii++)
						{
							in.readLine();
							nCount++;
						}
					}
				}
			}
		}
		tt++;
		if(tt ==1) break;		// use the second data

	}

	file.close();
	/*
	int nTimeStep = g_nTimeStep;

	QFile file(_strFile);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	QTextStream in(&file);
	while (!in.atEnd()) {
		for (size_t t = 0; t <= nTimeStep; t++)
		{
			for (int i = 0; i < _nEnsembleLen; i++)
			{
				QString line = in.readLine();
				for (int j = 0; j < _nLen; j++)
				{
					QString line = in.readLine();
					if (t == nTimeStep)
						_pData->SetData(i, j, line.toFloat());
				}
			}
		}
		break;

	}

	file.close();
	*/
}

void MeteModel::doPCA() {
	int nPCALen = _nEnsembleLen - 1;

	// 1.pca
	MyPCA pca;
	double* arrData = new double[_nEnsembleLen*nPCALen];
	pca.DoPCA(_pSDF->GetData(), arrData, _nEnsembleLen, _nLen, g_bNewData);

	// try to recover the data, used to test DoPCA and Recover
//	pca.Recover(arrData, _gridSDF, _nEnsembleLen);
//	pca.Recover(arrData, _gridSDF, 1);

	// 2.clustering
//	AHCClustering clustering;
	KMeansClustering clustering;
	clustering.DoCluster(_nEnsembleLen, nPCALen, _nClusters, arrData, _arrLabels);

	// 3.counting length of each cluster
	for (size_t i = 0; i < _nClusters; i++)
	{
		_listClusterLen.push_back(0);
	}
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_listClusterLen[_arrLabels[i]]++;
	}

	// 4.calculate mean of clusters
	// mean of each cluster
	double* _arrMeanC = new double[nPCALen*_nClusters];
	for (size_t i = 0; i < nPCALen*_nClusters; i++)
	{
		_arrMeanC[i] = 0;
	}
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		int nCluster = _arrLabels[i];
		for (size_t j = 0; j < nPCALen; j++)
		{
			_arrMeanC[nCluster*nPCALen + j] += arrData[i*nPCALen + j];
		}
	}
	for (size_t i = 0; i < _nClusters; i++)
	{
		for (size_t j = 0; j < nPCALen; j++)
		{
			_arrMeanC[i*nPCALen + j] /= _listClusterLen[i];
		}
	}

	_pClusterCenter = new double[_nClusters*_nLen];
	// 5.project back of the cluster means
	// transfer back
	pca.Recover(_arrMeanC, _pClusterCenter, _nClusters);

	delete[] _arrMeanC;


	// 4.output labels
	if (false)
	{
		cout << "Labels:" << endl;
		for (size_t i = 0; i < _nEnsembleLen; i++)
		{
			cout << _arrLabels[i] << endl;
		}
		cout << "~Labels:" << endl;

	}
	
	delete[] arrData;
}

void MeteModel::calculateSDF(const double* arrData, double* arrSDF, int nW, int nH, double isoValue, QList<ContourLine> contour) {
	for (size_t i = 0; i < nH; i++)
	{
		for (size_t j = 0; j < nW; j++)
		{
			double dbMinDis = nW*nH;
			for (size_t k = 0,lenK=contour.length(); k < lenK; k++)
			{
				for (size_t l = 0,lenL=contour[k]._listPt.length()-1; l < lenL; l++)
				{
					double dbDis = PointToSegDist(j,i
						, contour[k]._listPt[l].x(), contour[k]._listPt[l].y()
						, contour[k]._listPt[l+1].x(), contour[k]._listPt[l+1].y());
					if (dbDis < dbMinDis)
					{
						dbMinDis = dbDis;
					}

				}
			}
			if (arrData[i*_nWidth + j] < isoValue) dbMinDis = -dbMinDis;	// sign
			arrSDF[i*_nWidth + j] = dbMinDis;
		}
	}
}

GLubyte* MeteModel::generateTexture() {
	if (_bgFunction==bg_cluster)
	{
		// using cluster lable
		return generateTextureGridCluster();
	}
	if (_bgFunction == bg_sdf)
	{
		// using cluster lable
		return generateTextureSDF();
	}
	// using mean or variance
	const double* pData = _bgFunction == bg_mean ? _pData->GetMean() : _pData->GetVari();
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
			dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
			dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
		}
	}
	return dataTexture;
}

GLubyte* MeteModel::generateTextureMean() {
	const double* pData = _bgFunction==bg_mean? _pData->GetMean(): _pData->GetVari();
	// old version, whole area
	if (false)
	{
		GLubyte* dataTexture = new GLubyte[4 * _nLen];

		// color map
		ColorMap* colormap = ColorMap::GetInstance();

		for (int i = 0; i < _nLen; i++)
		{
			MYGLColor color = colormap->GetColor(pData[i]);
			// using transparency and the blue tunnel
			dataTexture[4 * i + 0] = color._rgb[0];
			dataTexture[4 * i + 1] = color._rgb[1];
			dataTexture[4 * i + 2] = color._rgb[2];
			dataTexture[4 * i + 3] = (GLubyte)255;
		}

		return dataTexture;
	}
	// new version, focused area
	if (true)
	{
		GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

		// color map
		ColorMap* colormap = ColorMap::GetInstance();
		for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
			for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

				int nIndex = i*_nWidth + j;
				int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
				MYGLColor color = colormap->GetColor(pData[nIndex]);
				// using transparency and the blue tunnel
				dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
				dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
				dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

			}
		}
		return dataTexture;
	}
}

GLubyte* MeteModel::generateTextureVar() {
	const double* pData = _bgFunction == bg_mean ? _pData->GetMean() : _pData->GetVari();
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]>_dbVarThreshold)
			{
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
			}
			else {
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)128;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)100;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

			}
			/*
			MYGLColor color = colormap->GetColor((int)pData[nIndex]);
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
			dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
			dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
*/
		}
	}
	return dataTexture;
}


GLubyte* MeteModel::generateTextureDiscreteSummary() {
	const double* pData = _pData->GetVari();
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];


	ColorMap* colormap = ColorMap::GetInstance(ColorMap::CP_RainBow);

	// clustering
	bool bDBSCAN = true;
	if (bDBSCAN) {
		// 1.create input data for dbscan
		int* arrState = new int[_nLen];
		int* arrGridLabels= new int[_nLen];
		// generate big variance points
		for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
			for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
				int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
				int nIndex = i*_nWidth + j;
				arrState[nIndex] = 1;
				if (pData[nIndex]>_dbVarThreshold)
					arrState[nIndex] = 1;
				else
					arrState[nIndex] = 0;

			}
		}

		// 2.dbscan
		DBSCANClustering dbscan;
		dbscan.DoCluster(_nFocusH, _nFocusW, g_nMinPts, g_dbEps, arrState, arrGridLabels);
		// 3.count the clusters
		int nClusters = 0;
		for (size_t i = 0; i < _nLen; i++)
		{
			if (arrGridLabels[i] > nClusters) nClusters = arrGridLabels[i];
		}

		// generate big variance points
		for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
			for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
				int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
				int nIndex = i*_nWidth + j;
				if (arrGridLabels[nIndex]>0)
				{
					MYGLColor color = colormap->GetColor(arrGridLabels[nIndex]-7);

					dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
					dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
					dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
					dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
				}
				else {
					dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
					dataTexture[4 * nIndexFocus + 1] = (GLubyte)255;
					dataTexture[4 * nIndexFocus + 2] = (GLubyte)255;
					dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

				}

			}
		}

		delete[]arrState;
		delete[]arrGridLabels;
	}
	else {
		std::vector<int> vecX;
		std::vector<int> vecY;

		// generate big variance points
		for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
			for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
				int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
				int nIndex = i*_nWidth + j;
				if (pData[nIndex]>_dbVarThreshold)
				{
					vecX.push_back(i);
					vecY.push_back(j);
				}

				dataTexture[4 * nIndexFocus + 0] = (GLubyte)100;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)100;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)100;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

			}
		}
		//	AHCClustering clustering;
		KMeansClustering clustering;
		int number = vecX.size();
		int nClusters = 4;
		double* arrBuf = new double[number * 2];
		for (size_t i = 0; i < number; i++)
		{
			arrBuf[i * 2] = vecX[i];
			arrBuf[i * 2 + 1] = vecY[i];
		}
		int* arrLabel = new int[number];
		clustering.DoCluster(number, 2, nClusters, arrBuf, arrLabel);

		for (size_t i = 0; i < number; i++)
		{

			int nIndexFocus = (vecX[i] - _nFocusY)*_nFocusW + vecY[i] - _nFocusX;

			switch (arrLabel[i]) {
			case 0:
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
				break;
			case 1:
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)255;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
				break;
			case 2:
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)255;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
				break;
			case 3:
				dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
				dataTexture[4 * nIndexFocus + 1] = (GLubyte)128;
				dataTexture[4 * nIndexFocus + 2] = (GLubyte)128;
				dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
				break;
			}
		}

		delete arrBuf;
		delete arrLabel;
	}

	return dataTexture;
}

vector<double> MeteModel::GetVariance() {

	const double* pData = _pData->GetVari();
	vector<double> vecVar;
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			vecVar.push_back(pData[nIndex]);
		}
	}
	std::sort(vecVar.begin(), vecVar.end());
	return vecVar;
}

GLubyte* MeteModel::generateTextureGridCluster() {
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	GLubyte colors[40][3] = {
		{ 31	,119	,180 },
		{ 174	,199	,232 },
		{ 255	,127	,14 },
		{ 255	,187	,120 },
		{ 44	,160	,44 },
		{ 152	,223	,138 },
		{ 214	,39		,40 },
		{ 255	,152	,150 },
		{ 148	,103	,189 },
		{ 197	,176	,213 },
		{ 140	,86		,75 },
		{ 196	,156	,148 },
		{ 227	,119	,194 },
		{ 247	,182	,210 },
		{ 127	,127	,127 },
		{ 188	,189	,34 },
		{ 219	,219	,141 },
		{ 23	,190	,207 },
		{ 158	,218	,229 },
		{ 57	,59		,121 },
		{ 82	,84		,163 },
		{ 107	,110	,207 },
		{156	,158	,222 },
		{99		,121	,57	 },
		{140	,162	,82	 },
		{181	,207	,107 },
		{206	,219	,156 },
		{140	,109	,49	 },
		{189	,158	,57	 },
		{231	,186	,82	 },
		{231	,203	,148 },
		{132	,60		,57	 },
		{173	,73		,74	 },
		{214	,97		,107 },
		{231	,150	,156 },
		{123	,65		,115 },
		{165	,81		,148 },
		{206	,109	,189 },
		{222	,158	,214 }
	};

	double fMin = 100000;
	double fMax = -100000;

	for (int i = 0; i < _nLen; i++)
	{
		double f = _arrGridLabels[i];
		if (f > fMax) fMax = f;
		if (f < fMin) fMin = f;
	}
	double fRange = fMax - fMin;

	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;

			GLbyte color[3];
			switch (_arrGridLabels[nIndex])
			{
			case -3:
				color[0] = 255;
				color[1] = 255;
				color[2] = 255;
				break;
			case -2:
				color[0] = 220;
				color[1] = 220;
				color[2] = 220;
				break;
			case -1:
				color[0] = 200;
				color[1] = 200;
				color[2] = 200;
				break;
			default:
				if (_arrGridLabels[nIndex]<40)
				{
					color[0] = colors[_arrGridLabels[nIndex]][0];
					color[1] = colors[_arrGridLabels[nIndex]][1];
					color[2] = colors[_arrGridLabels[nIndex]][2];
				}
				else {
					color[0] = 0;
					color[1] = 0;
					color[2] = 0;

				}
				break;
			}
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = (GLubyte)color[0];
			dataTexture[4 * nIndexFocus + 1] = (GLubyte)color[1];
			dataTexture[4 * nIndexFocus + 2] = (GLubyte)color[2];
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
		}
	}

	return dataTexture;
}

GLubyte* MeteModel::generateTextureRange() {
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	double* pMeanBias = new double[_nLen];
	for (size_t i = 0; i < _nLen; i++)
	{
		pMeanBias[i] = _pData->GetMean()[i];// -_gridMean[2][i];
	}

	double *pMean = pMeanBias;

	double fMin = 100000;
	double fMax = -100000;

	for (int i = 0; i < _nLen; i++)
	{
		double f = _pData->GetData(0, i);
		if (pMean[i] > fMax) fMax = pMean[i];
		if (pMean[i] < fMin) fMin = pMean[i];
	}
	double fRange = fMax - fMin;

	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 2] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
		}
	}

	delete[] pMeanBias;
	return dataTexture;
}

GLubyte* MeteModel::generateTextureRange(int nIndex) {
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	const double *pMean = _pData->GetLayer(nIndex);

	double fMin = 100000;
	double fMax = -100000;

	for (int i = 0; i < _nLen; i++)
	{
		double f = _pData->GetData(0, i);
		if (pMean[i] > fMax) fMax = pMean[i];
		if (pMean[i] < fMin) fMin = pMean[i];
	}
	double fRange = fMax - fMin;

	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 2] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
		}
	}

	return dataTexture;
}

GLubyte* MeteModel::generateTextureSDF() {
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

	double fMin = 100000;
	double fMax = -100000;
	const double* pData = _pSDF->GetMean();
//	double* pData = _pClusterCenter + 2 * _nLen;	// using the third cluster

	for (int i = 0; i < _nLen; i++)
	{
		if (pData[i] > fMax) fMax = pData[i];
		if (pData[i] < fMin) fMin = pData[i];
	}
	double fRange = fMax - fMin;

	for (int i = 0; i < _nFocusH; i++) {
		for (int j = 0; j < _nFocusW; j++) {

			int nIndex = (i + _nFocusY)*_nWidth + j + _nFocusX;
			int nIndexFocus = i*_nFocusW + j;
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 2] = (GLubyte)((pData[nIndex] - fMin) * 254 / fRange);
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)((pData[nIndex] - fMin) * 254 / fRange);
		}
	}

	return dataTexture;
}

void MeteModel::generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas) {
	ContourStripGenerator generator;
	generator.Generate(areas, contourMin, contourMax, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
}

void MeteModel::readData() {
	QFile file(_strFile);
	if (!file.open(QIODevice::ReadOnly)) {
		return;
	}

	char* temp = new char[file.size()];
	file.read(temp, file.size());
	float* f = (float*)temp;
	int x = file.size();

	for (int l = 0; l < _nEnsembleLen; l++)
	{
		// 第一个集合之前多一个数据，之后的集合之前多两个数据
		if (l == 0) f++;
		else f += 2;

		for (int i = 0; i < _nHeight; i++)
		{
			for (int j = 0; j < _nWidth; j++)
			{
				double fT = *f;
				_pData->SetData(l, i, j, *f++);
				// 只取整度，过滤0.5度
				if (_bFilter&&j < _nWidth - 1) f++;
			}
			// 只取整度，过滤0.5度
			if (_bFilter&&i < _nHeight - 1) f += (_nWidth * 2 - 1);
		}
	}
}

void MeteModel::doSpatialClustering() {
	// 1.create input data for dbscan
	int* arrState = new int[_nLen];
	for (size_t i = 0; i < _nLen; i++)
	{
		arrState[i] = _pData->GetMean()[i]>g_fThreshold ? 1 : 0;
	}

	// 2.dbscan
	DBSCANClustering dbscan;
	dbscan.DoCluster(_nFocusH, _nFocusW, g_nMinPts, g_dbEps, arrState, _arrGridLabels);
	// 3.count the clusters
	int nClusters = 0;
	for (size_t i = 0; i < _nLen; i++)
	{
		if (_arrGridLabels[i] > nClusters) nClusters = _arrGridLabels[i];
	}
	bool bRaw = false;	// weather using raw data or belief eclipse
	nClusters++;
	_points.clear();
	if (bRaw) {
		for (size_t i = 0; i < _nHeight; i++)
		{
			for (size_t j = 0; j < _nWidth; j++)
			{
				if (_arrGridLabels[i*_nWidth + j] == 1) {

					_points.push_back(Point(DPoint3(j, i, 0)));
				}
			}
		}
	}
	else {
		if (g_nConfidenceEllipseIndex>0)
		{
			if (g_nConfidenceEllipseIndex < nClusters) {
				MyPCA::generateEllipse(_points, _arrGridLabels, g_nConfidenceEllipseIndex, _nWidth, _nHeight,g_dbMDis);
			}
		}
		else {
			for (size_t i = 0; i < nClusters; i++)
			{
				MyPCA::generateEllipse(_points, _arrGridLabels, i, _nWidth, _nHeight, g_dbMDis);
			}
		}
	}
	delete[] arrState;
}

void MeteModel::initializeModel() {
}

void MeteModel::Brush(int nLeft, int nRight, int nTop, int nBottom) {
	_listContourBrushed.clear();
	_listContourNotBrushed.clear();

	int nFocusL = _nWest + 180;
	int nFocusB = _nSouth + 90;

	int nFocusW = _nEast - _nWest;
	int nFocusH = _nNorth - _nSouth;

	nLeft -= nFocusL;
	nRight -= nFocusL;
	nTop -= nFocusB;
	nBottom -= nFocusB;

	if (nLeft < 0) nLeft = 0;
	if (nBottom < 0)nBottom = 0;
	if (nRight > nFocusW) nRight = nFocusW;
	if (nTop > nFocusH) nTop = nFocusH;

	for (size_t l = 0; l < _nEnsembleLen; l++)
	{
		bool bCover = false;
		bool bHigher = false;
		bool bLower = false;

		for (int i = nBottom; i <= nTop; i++)
		{
			for (int j = nLeft; j <= nRight; j++)
			{
//				std::cout << _pData->GetData(l, i, j) << endl;
				_pData->GetData(l, i, j)>g_fThreshold ? (bHigher = true) : (bLower = true);
				if (bHigher&&bLower) {
					bCover = true;
					break;
				}
			}
			if (bCover)
			{
				break;
			}
		}
		if (bCover)
		{
			_listContourBrushed.push_back(_listContour[l]);
		}
		else _listContourNotBrushed.push_back(_listContour[l]);
	}
}

QList<QList<ContourLine>> MeteModel::GetContourBrushed()
{
	return _listContourBrushed;
}

QList<QList<ContourLine>> MeteModel::GetContourNotBrushed()
{
	return _listContourNotBrushed;
}

QList<QList<ContourLine>> MeteModel::GetContour()
{
	return _listContour;
}


MeteModel* MeteModel::CreateModel() {
	if (g_bEnsembleModel) {

		MeteModel* pModel = NULL;
		int nWidth = g_globalW;
		int nHeight = g_globalH;

		int nFocusX = 0;
		int nFocusY = 0;
		int nFocusW = nWidth;
		int nFocusH = nHeight;

		int nWest;
		int nEast;
		int nNorth;
		int nSouth;
		int nFocusWest;
		int nFocusEast;
		int nFocusNorth;
		int nFocusSouth;



		if (g_bGlobalArea)
		{
			nFocusX = 0;
			nFocusY = 0;
			nFocusW = nWidth;
			nFocusH = nHeight;
		}
		else {

			nWidth = 31;
			nHeight = 31;

			if (g_bSubArea)
			{
				nFocusX = 0;
				nFocusY = 40;
				nFocusW = 31;
				nFocusH = 11;
				nWest = 60;
				nEast = 150;
				nSouth = 10;
				nNorth = 60;


				nFocusWest = 60;
				nFocusEast = 90;
				nFocusSouth = 50;
				nFocusNorth = 60;
			}
			else {
				nFocusX = 0;
				nFocusY = 0;
				nFocusW = nWidth;
				nFocusH = nHeight;
				nWest = 100;
				nEast = 130;
				nSouth = 20;
				nNorth = 50;
				nFocusWest = 100;
				nFocusEast = 130;
				nFocusSouth = 20;
				nFocusNorth = 50;
			}

		}



		pModel = new EnsembleModel();
		pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
			//				, "../../data/t2-mod-ecmwf-20160105-00-72-216.txt", false
			//				, "../../data/t2-2007-2017-jan-120h-50.txt", false
			//, "../../data/t2-2007-2017-jan-144 and 240h-50.txt", false
			, "../../data/t2-mod-ecmwf-200701-00-360.txt", false
			, nWest, nEast, nSouth, nNorth
			, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth, g_bFilter);

		return pModel;
	}
	else {
		// old codes
		MeteModel* pModel = NULL;
		int nWidth = g_globalW;
		int nHeight = g_globalH;

		int nFocusX = 0;
		int nFocusY = 0;
		int nFocusW = nWidth;
		int nFocusH = nHeight;

		int nWest;
		int nEast;
		int nNorth;
		int nSouth;
		int nFocusWest;
		int nFocusEast;
		int nFocusNorth;
		int nFocusSouth;



		if (g_bGlobalArea)
		{
			nFocusX = 0;
			nFocusY = 0;
			nFocusW = nWidth;
			nFocusH = nHeight;
		}
		else {
			if (g_bFilter) {
				nWidth = 91;
				nHeight = 51;
			}
			else {
				nWidth = 181;
				nHeight = 101;
			}
			if (g_bSubArea)
			{
				nFocusX = 0;
				nFocusY = 40;
				nFocusW = 31;
				nFocusH = 11;
				nWest = 60;
				nEast = 150;
				nSouth = 10;
				nNorth = 60;


				nFocusWest = 60;
				nFocusEast = 90;
				nFocusSouth = 50;
				nFocusNorth = 60;
			}
			else {
				nFocusX = 0;
				nFocusY = 0;
				nFocusW = nWidth;
				nFocusH = nHeight;
				nWest = 60;
				nEast = 150;
				nSouth = 10;
				nNorth = 60;
				nFocusWest = 60;
				nFocusEast = 150;
				nFocusSouth = 10;
				nFocusNorth = 60;
			}

		}




		bool bNewData = true;

		switch (g_usedModel)
		{
		case PRE_CMA:
			pModel = new EnsembleModel();
			pModel->InitModel(14, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-cma-20160802-00-96.txt"); break;
		case PRE_CPTEC:
			pModel = new EnsembleModel();
			pModel->InitModel(14, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-cptec-20160802-00-96.txt"); break;
		case PRE_ECCC:
			pModel = new EnsembleModel();
			pModel->InitModel(20, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-eccc-20160802-00-96.txt"); break;
		case PRE_ECMWF:
			pModel = new MeteModel();
			pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-ecmwf-20160802-00-96.txt"); break;
		case PRE_JMA:
			pModel = new EnsembleModel();
			pModel->InitModel(26, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-jma-20160802-00-96.txt"); break;
		case PRE_KMA:
			pModel = new EnsembleModel();
			pModel->InitModel(24, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-kma-20160802-00-96.txt"); break;
		case PRE_NCEP:
			pModel = new EnsembleModel();
			pModel->InitModel(20, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH, "../../data/data10/pre-mod-ncep-20160802-00-96.txt"); break;
		case T2_ECMWF:
			pModel = new ContourBandDepthModel();
			if (bNewData) {
				pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
					//				, "../../data/t2-mod-ecmwf-20160105-00-72-216.txt", false
					//				, "../../data/t2-2007-2017-jan-120h-50.txt", false
					, "../../data/t2-2007-2017-jan-144 and 240h-50.txt", false
					, nWest, nEast, nSouth, nNorth
					, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth, g_bFilter);
			}
			else {

				pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
					, "../../data/t2_mod_20080101-96h.dat", true
					, nWest, nEast, nSouth, nNorth
					, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth, g_bFilter);
			}
			break;
		case T2_Reanalysis:
			pModel = new ReanalysisModel();
			pModel->InitModel(1209, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
				, "../../data/t2_obs_1979-2017_1_china.txt", false
				, nWest, nEast, nSouth, nNorth
				, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth, false);
			break;
		defaut:
			break;
		}

		return pModel;
	}
/*

	*/
}