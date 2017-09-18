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

#include "SpatialCluster.h"


using namespace std;

/* 
	allign these result
	nK: number of clusters
	arrMap: map of the ensemble members to its new sequence
*/
void Align(ClusterResult* arrResult, int nK,int* arrMap) {	

	// initialize the map to natral sequence
	for (size_t i = 0; i < g_nEnsembles; i++) arrMap[i] = i;

	// swap sorting
	for (size_t i = 0; i < g_nEnsembles-1; i++)
	{
		for (size_t j = i + 1; j < g_nEnsembles; j++) {
			bool bSwap = false;
			for (size_t k = 0; k < nK; k++)
			{
				if (arrResult[k]._arrLabels[arrMap[i]]>arrResult[k]._arrLabels[arrMap[j]]) {
					bSwap = true;
					break;
				}
				else if (arrResult[k]._arrLabels[arrMap[i]] < arrResult[k]._arrLabels[arrMap[j]]) {
					break;
				}
			}
			if (bSwap)
			{
				int nTemp = arrMap[i];
				arrMap[i] = arrMap[j];
				arrMap[j] = nTemp;
			}
		}
	}
}

void ClusterResult::PushLabel(int nIndex, int nLabel) {
	_vecItems[nLabel].push_back(nIndex);
	_arrLabels[nIndex] = nLabel;
}


void ClusterResult::Match(ClusterResult& mc) {
	int mtxWeight[5][5];					// weight matrix
	int arrMatch[5];						// the result of the perfect match
	int nClusters = 5;						// number of clusters
	// 0.initialize matrix
	for (size_t i = 0; i < nClusters; i++)
	{
		for (size_t j = 0; j < nClusters; j++) {
			mtxWeight[i][j] = 0;
		}
	}
	// 1.calculate weight
	for (size_t i = 0; i < _nM; i++)
	{
		mtxWeight[_arrLabels[i]][mc._arrLabels[i]]++;
	}

	for (size_t i = 0; i < nClusters; i++)
	{
		for (size_t j = 0; j < nClusters; j++) {
			cout<< mtxWeight[i][j]<<"\t";
		}
		cout << endl;
	}

	// 2.calculate max weight perfect matching using Hungarian algorithm
	// use approximate algorithm, always searching for the biggest value
	for (size_t l = 0; l < nClusters; l++) {
		// find biggest value
		int iMax = 0; 
		int jMax = 0;
		int vMax = -1;
		for (size_t i = 0; i < nClusters; i++)
		{
			for (size_t j = 0; j < nClusters; j++) {
				if (mtxWeight[i][j]>vMax) {
					vMax = mtxWeight[i][j];
					iMax = i;
					jMax = j;
				}
			}
		}
		arrMatch[iMax] = jMax;
		for (size_t i = 0; i < nClusters; i++)
		{
			mtxWeight[i][jMax] = mtxWeight[iMax][i] = -1;
		}
	}
	/*
	qDebug() << "match:";
	for (size_t i = 0; i < nClusters; i++)
	{
		qDebug() << arrMatch[i];
	}
	*/
	// 3.reset labels
	for (size_t i = 0; i < _nM; i++)
	{
		_arrLabels[i] = arrMatch[_arrLabels[i]];
	}
	generateItemsByLabels();

}


void ClusterResult::AlighWith(int* arrMap) {
	for (size_t i = 0; i < _nK; i++)
	{
		int nLen = _vecItems[i].size();
		for (size_t j = 0; j < nLen-1; j++)
		{
			for (size_t k = j + 1; k < nLen; k++) {
				if (arrMap[_vecItems[i][j]]>arrMap[_vecItems[i][k]])
				{
					int nTemp = _vecItems[i][j];
					_vecItems[i][j] = _vecItems[i][k];
					_vecItems[i][k] = nTemp;
				}
			}
		}
	}
}

void ClusterResult::generateItemsByLabels() {
	for (size_t i = 0; i < 5; i++)
	{
		_vecItems[i].clear();
	}
	for (size_t i = 0; i < _nM; i++)
	{
		_vecItems[_arrLabels[i]].push_back(i);
	}
}

void ClusterResult::Sort() {
	// array recording count of each cluster
	int arrCount[g_nClusterMax];
	for (size_t i = 0; i < g_nClusterMax; i++)
	{
		arrCount[i] = _vecItems[i].size();
	}
	// get the projection
	int arrMap[g_nClusterMax];
	for (size_t i = 0; i < 5; i++)
	{
		int nMax = -1;
		int nMaxIndex = -1;
		for (size_t j = 0; j < g_nClusterMax; j++)
		{
			if (arrCount[j] > nMax) {
				nMax = arrCount[j];
				nMaxIndex = j;
			}
		}
		arrMap[nMaxIndex] = i;
		arrCount[nMaxIndex] = -1;
	}
	// reset label
	for (size_t i = 0; i < _nM; i++)
	{
		_arrLabels[i] = arrMap[_arrLabels[i]];
	}
	generateItemsByLabels();
}

void ClusterResult::Reset(int nM, int nK) {
	_nM = nM;
	_nK = nK;
	for (int i = 0; i < 5; i++) {
		_vecItems[i].clear();
	}
}

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

MeteModel::MeteModel()
{
}

MeteModel::~MeteModel()
{
	if (_pData)
	{
		delete _pData;
	}


	for each (UnCertaintyArea* pArea in _listUnionAreaE)
		delete pArea;



	for (size_t i = 0; i < _listUnionAreaEG.length(); i++)
	{
		for each (UnCertaintyArea* pArea in _listUnionAreaEG[i])
			delete pArea;
	}
	if (_dataTexture) delete[] _dataTexture;
}

void MeteModel::InitModel(int nEnsembleLen, int nWidth, int nHeight, int nFocusX, int nFocusY, int nFocusW, int nFocusH
	, QString strFile, bool bBinary, int nWest, int nEast, int nSouth, int nNorth
	, int nFocusWest, int nFocusEast, int nFocusSouth, int nFocusNorth,bool bFilter) {
	// 0.record states variables
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


	// 5.generate features
	ContourGenerator generator;
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		QList<ContourLine> contour;
		generator.Generate(_pData->GetLayer(i), contour, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		_listContour.push_back(contour);
	}




		generator.Generate(_pData->GetUMin(), _listContourMinE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pData->GetUMax(), _listContourMaxE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		generator.Generate(_pData->GetMean(), _listContourMeanE, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);


	generateContourImp(_listContourMinE, _listContourMaxE, _listUnionAreaE);


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


void MeteModel::buildTextureThresholdVariance() {
	const double* pData = _bgFunction == bg_mean ? _pData->GetMean() : _pData->GetVari(_nSmooth);

	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]>_dbVarThreshold)
			{
				_dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
				_dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;
			}
			else {
				_dataTexture[4 * nIndexFocus + 0] = (GLubyte)128;
				_dataTexture[4 * nIndexFocus + 1] = (GLubyte)100;
				_dataTexture[4 * nIndexFocus + 2] = (GLubyte)0;
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

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
}

// cluster comparison function, used for sort
bool ClusterComparison(OneSpatialCluster c1, OneSpatialCluster c2) {
	return c1._nArea > c2._nArea;
}

// generate texture of clustered variance
void MeteModel::buildTextureClusteredVariance() {
	// 1.find uncertainty regions and sort them according to the areas
	const double* pData = _pData->GetVari(_nSmooth);
	SpatialCluster cluster;
	std::vector<OneSpatialCluster> result;
	cluster.DoCluster(pData, _nWidth, _nHeight, _dbVarThreshold,result);
	sort(result.begin(), result.end(), ClusterComparison);


	// 2.cluster in each uncertainty region
	for (size_t i = 0,length=std::min((int)result.size(), _nUncertaintyRegions); i < length; i++)
	{
		clusterSpatialArea(result[i], _vecPCAPoints[i], _arrClusterResult[i]);
		if (i==0)
		{
			// first region, sort the clusters according to the number of their members
			_arrClusterResult[i].Sort();
		}
		else {
			// the following region match the region before
			_arrClusterResult[i].Match(_arrClusterResult[i - 1]);
		}
	}

	// 3.align the results
	int arrMap[g_nEnsembles];
	Align(_arrClusterResult, _nUncertaintyRegions, arrMap);

	int arrMap2[g_nEnsembles];
	for (size_t i = 0; i < g_nEnsembles; i++)
	{
		arrMap2[arrMap[i]] = i;
	}

	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		_arrClusterResult[i].AlighWith(arrMap2);
	}




	// 4.reset labels for pca
	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		setLabelsForPCAPoints(_vecPCAPoints[i], _arrClusterResult[i]);
	}

	// calculate the similarity between different uncertainty regions
	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		for (size_t j = 0; j < _nUncertaintyRegions; j++) {
			_matrixSimilarity[i][j] = 0;
			for (size_t k = 0; k < _nClusters; k++)
			{
				_matrixSimilarityC[i][j][k] = 0;
			}
		}
	}
	for (size_t l = 0; l < _nEnsembleLen; l++)
	{
		for (size_t i = 0; i < _nUncertaintyRegions; i++)
		{
			for (size_t j = 0; j < i; j++) {
				if (_arrClusterResult[i]._arrLabels[l] == _arrClusterResult[j]._arrLabels[l])
				{
					_matrixSimilarity[i][j]++;
					_matrixSimilarityC[i][j][_arrClusterResult[i]._arrLabels[l]]++;
				}
			}
		}
	}


	// 5.initialize the texture
	for (size_t i = 0; i < _nFocusLen; i++)
	{
		_dataTexture[4 * i + 0] = (GLubyte)100;
		_dataTexture[4 * i + 1] = (GLubyte)100;
		_dataTexture[4 * i + 2] = (GLubyte)100;
		_dataTexture[4 * i + 3] = (GLubyte)255;
	}

	// 6.set the cluster color
	for (size_t i = 0, length = result.size(); i < length; i++)
	{
		
		for (size_t j = 0, length = result[i]._vecPoints.size(); j < length; j++) {

			int nIndex = result[i]._vecPoints[j].x*_nWidth + result[i]._vecPoints[j].y;

			_dataTexture[4 * nIndex + 0] = (GLubyte)ColorMap::GetCategory20I(i, 0);
			_dataTexture[4 * nIndex + 1] = (GLubyte)ColorMap::GetCategory20I(i, 1);
			_dataTexture[4 * nIndex + 2] = (GLubyte)ColorMap::GetCategory20I(i, 2);
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
		}
	}
}

// generate texture of colormap of mean or variance
void MeteModel::buildTextureColorMap() {
	const double* pData = _bgFunction == bg_mean ? _pData->GetMean() : _pData->GetVari();
	_dataTexture = new GLubyte[4 * _nFocusLen];

	ofstream output("variance.txt");
	output << (_bgFunction == bg_mean ? "mean" : "variance") << endl;
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			output << j << "\t" << i << "\t" << pData[nIndex] << endl;
		}
	}

	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
			_dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
			_dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
			_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

		}
	}
}

// generate texture of smoothed variance
void MeteModel::buildTextureSmoothedVariance()
{
	const double* pData = _pData->GetVari(_nSmooth);
	_dataTexture = new GLubyte[4 * _nFocusLen];

	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndexFocus + 0] = color._rgb[0];
			_dataTexture[4 * nIndexFocus + 1] = color._rgb[1];
			_dataTexture[4 * nIndexFocus + 2] = color._rgb[2];
			_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

		}
	}
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

void MeteModel::generatePCAPoint(OneSpatialCluster& cluster, std::vector<DPoint3>& points) {
	// 0.clear the points
	points.clear();
	// 1.set parameter
	int mI = cluster._nArea;
	int mO = 2;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = new double[mI*n];
	double* arrOutput = new double[mO*n];
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < mI; j++) {
			int x = cluster._vecPoints[j].x;
			int y = cluster._vecPoints[j].y;
			arrInput[i*mI + j] = _pData->GetData(i, x, y);
		}
	}
	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);
	// 4.generate points from the output
	for (size_t i = 0; i < n; i++)
	{
		points.push_back(DPoint3(arrOutput[i * 2], arrOutput[i * 2 + 1], 0));
	}
	// 5.release the buffer
	delete[] arrInput;
	delete[] arrOutput;
}

void MeteModel::clusterSpatialArea(OneSpatialCluster& cluster, std::vector<DPoint3>& points, ClusterResult& cr) {
	// 1.pca
	generatePCAPoint(cluster, points);

	// 2.cluster
	CLUSTER::Clustering* pClusterer = new CLUSTER::KMeansClustering();
	int nN = _nEnsembleLen;			// number of data items
	int nK = _nClusters;						// clusters
	int* arrLabel = new int[nN];
	bool bClusterUsingPCA = false;	// whether using pca points to cluster
	if (bClusterUsingPCA)
	{
		int nM = 2;					// dimension
		double* arrBuf = new double[nN * nM];
		for (size_t i = 0; i < nN; i++)
		{
			arrBuf[i*nM] = points[i].x;
			arrBuf[i*nM + 1] = points[i].y;
		}
		pClusterer->DoCluster(nN, nM, nK, arrBuf, arrLabel);
		delete arrBuf;
	}
	else {
		int nM = cluster._nArea;		// dimension
		double* arrBuf = new double[nN * nM];
		for (size_t i = 0; i < nN; i++)
		{
			for (size_t j = 0; j < nM; j++)
			{
				int x = cluster._vecPoints[j].x;
				int y = cluster._vecPoints[j].y;
				arrBuf[i * nM + j] = _pData->GetData(i, x, y);
			}
		}
		pClusterer->DoCluster(nN, nM, nK, arrBuf, arrLabel);
		delete arrBuf;
	}


	// 3.set label for pca points
	// no use, the setting is after this function
	for (size_t i = 0; i < nN; i++)
	{
		points[i].z = arrLabel[i] + .5;
	}

	// 4.record label
	cr.Reset(_nEnsembleLen, nK);
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		cr.PushLabel(i, arrLabel[i]);
	}

	// 5.release the resouse
	delete arrLabel;
	delete pClusterer;
}

void MeteModel::setLabelsForPCAPoints(std::vector<DPoint3>& points, ClusterResult&cr) {
	// 3.set label for pca points
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		points[i].z = cr._arrLabels[i] + .5;
	}
}

void MeteModel::SetVarThreshold(double dbThreshold)
{ 
	if (abs(dbThreshold - _dbVarThreshold) < 0.0001) return;

	_dbVarThreshold = dbThreshold; 

	regenerateTexture();
}

void MeteModel::SetSmooth(int nSmooth) {
	if (nSmooth == _nSmooth) return;
	_nSmooth = nSmooth;
	regenerateTexture();
}

void MeteModel::regenerateTexture() {
	switch (_bgFunction)
	{
	case MeteModel::bg_mean:
	case MeteModel::bg_vari:
		buildTextureColorMap();
		break;
	case MeteModel::bg_cluster:
		//		return generateClusteredVarianceTexture_old();
		buildTextureClusteredVariance();
		break;
	case MeteModel::bg_sdf:
		buildTextureThresholdVariance();
		break;
	case MeteModel::bg_vari_smooth:
		buildTextureSmoothedVariance();
		break;
	default:
		break;
	}

}

void MeteModel::SetBgFunction(enumBackgroundFunction f) 
{ 
	_bgFunction = f; 
	regenerateTexture();
}

void MeteModel::SetUncertaintyAreas(int nAreas)
{ 
	_nUncertaintyRegions = nAreas;
	regenerateTexture();
}

GLubyte* MeteModel::generateTextureNew() {
	if (!_dataTexture) {
		_dataTexture = new GLubyte[4 * _nFocusLen];
		regenerateTexture();
	}
	return _dataTexture;
}