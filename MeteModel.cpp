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
	, int nFocusWest, int nFocusEast, int nFocusSouth, int nFocusNorth) {
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

	// 2.allocate resource
	_pData = new DataField(_nWidth, _nHeight, _nEnsembleLen);

	// 3.read data
	if (_bBinaryFile)
	{
		readData();
	}
	else {
		if (g_bGlobalArea)
		{
			readDataFromTextG();
		}
		else {
			readDataFromText();
		}
	}
	// read dip value (temp)
//	readDipValue("../../data/t2-2007-2017-jan-144 and 240h-50(1Degree, single timestep,skip 3)_dipValue.txt");
//	readDipValue("../../data/t2-2007-2017-jan-144 and 240h-50(1Degree, single timestep)_dipValue.txt");
	readDipValue("../../data/t2-2007-2017-jan-144 and 240h-50(1Degree, single timestep)_dipValue_P.txt");

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
				}
			}
		}
		tt++;
		if(tt ==1) break;		// use the second data
	}

	file.close();
}


void MeteModel::readDipValue(char* strFileName) {

	QFile file(strFileName);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}


	QTextStream in(&file);
	// every grid
	for (int j = 0; j < _nLen; j++)
	{
		QString line = in.readLine();
//		double dbDipValue = line.toFloat() * 100;
		double dbDipValue = line.toFloat() * 10;
		_pData->SetDipValue(j, dbDipValue);
		//qDebug() << j<<"\t" << dbDipValue;
	}

	file.close();
}

void MeteModel::readDataFromTextG() {

	int nTimeStep = g_nTimeStep;

	QFile file(_strFile);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	QTextStream in(&file);
	int tt = 0;

	int nWidth = _nWidth - 1;				// width of a scalar field in the file
	int nLen = nWidth*_nHeight;				// length of a scalar field in the file
	while (!in.atEnd()) {
		// every time step
		for (size_t t = 0; t <= nTimeStep; t++)
		{
			// every ensemble member
			for (int i = 0; i < _nEnsembleLen; i++)
			{
				QString line = in.readLine();
				// every grid
				for (int j = 0; j < nLen; j++)
				{
					QString line = in.readLine();
					if (t == nTimeStep) {
						int r = j / nWidth;
						int c = j % nWidth;
						int nIndex = r*_nWidth + c+1;
						_pData->SetData(i, nIndex, line.toFloat());
						if (c==nWidth-1)
						{
							int nIndex = r*_nWidth;
							_pData->SetData(i, nIndex, line.toFloat());
						}
					}
				}
			}
		}
		tt++;
		if (tt == 1) break;		// use the second data
	}

	file.close();
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
	_nThresholdedGridPoints = 0;
	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]>_dbVarThreshold)
			{
				_dataTexture[4 * nIndexFocus + 0] = ColorMap::GetThresholdColorI(0);
				_dataTexture[4 * nIndexFocus + 1] = ColorMap::GetThresholdColorI(1);
				_dataTexture[4 * nIndexFocus + 2] = ColorMap::GetThresholdColorI(2);
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

				_nThresholdedGridPoints++;
			}
			else {
				_dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 1] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 2] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

			}
		}
	}
}


void MeteModel::buildTextureThresholdDipValue() {
	double dbThreshold = 7;
	const double* pData = _pData->GetDipValue();
	_nThresholdedGridPoints = 0;
	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]<dbThreshold)
			{
				_dataTexture[4 * nIndexFocus + 0] = ColorMap::GetThresholdColorI(0);
				_dataTexture[4 * nIndexFocus + 1] = ColorMap::GetThresholdColorI(1);
				_dataTexture[4 * nIndexFocus + 2] = ColorMap::GetThresholdColorI(2);
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

				_nThresholdedGridPoints++;
			}
			else {
				_dataTexture[4 * nIndexFocus + 0] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 1] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 2] = (GLubyte)255;
				_dataTexture[4 * nIndexFocus + 3] = (GLubyte)255;

			}
		}
	}
}

// cluster comparison function, used for sort
bool ClusterComparison(UncertaintyRegion c1, UncertaintyRegion c2) {
	return c1._nArea > c2._nArea;
}

// generate texture of clustered variance
void MeteModel::buildTextureClusteredVariance() {
	// 1.find uncertainty regions and sort them according to the areas
	const double* pData = _pData->GetVari(_nSmooth);
	SpatialCluster cluster;
	_vecRegions.clear();
	cluster.DoCluster(pData, _nWidth, _nHeight, _dbVarThreshold, _vecRegions);
	sort(_vecRegions.begin(), _vecRegions.end(), ClusterComparison);

	// cluster in each region
	regionCluster();

	// align them
	alignClusters();





	// 6.Generate texture
	// 6.1.initialize the texture
	for (size_t i = 0; i < _nFocusLen; i++)
	{
		_dataTexture[4 * i + 0] = (GLubyte)100;
		_dataTexture[4 * i + 1] = (GLubyte)100;
		_dataTexture[4 * i + 2] = (GLubyte)100;
		_dataTexture[4 * i + 3] = (GLubyte)255;
	}

	// 6.set the cluster color
	for (size_t i = 0, length = _vecRegions.size(); i < _nUncertaintyRegions; i++)
		//for (size_t i = 0, length = result.size(); i < length; i++)
	{
		
		for (size_t j = 0, length = _vecRegions[i]._vecPoints.size(); j < length; j++) {

			int nIndex = _vecRegions[i]._vecPoints[j].x*_nWidth + _vecRegions[i]._vecPoints[j].y;

			_dataTexture[4 * nIndex + 0] = (GLubyte)ColorMap::GetCategory10I(i, 0);
			_dataTexture[4 * nIndex + 1] = (GLubyte)ColorMap::GetCategory10I(i, 1);
			_dataTexture[4 * nIndex + 2] = (GLubyte)ColorMap::GetCategory10I(i, 2);
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
		}
	}
}

// generate texture of colormap of mean or variance
void MeteModel::buildTextureColorMap() {
	const double* pData;
	switch (_bgFunction)
	{
	case MeteModel::bg_mean:
		pData = _pData->GetMean();
		break;
	case MeteModel::bg_vari:
		pData = _pData->GetVari();
		break;
	case MeteModel::bg_dipValue:
		pData = _pData->GetDipValue();
		break;
	default:
		break;
	}

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

	const double* pData = _pData->GetVari(_nSmooth);
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
//				if (_bFilter&&j < _nWidth - 1) f++;
			}
			// 只取整度，过滤0.5度
//			if (_bFilter&&i < _nHeight - 1) f += (_nWidth * 2 - 1);
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
			, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth);

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

			nWidth = 91;
			nHeight = 51;
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
			nWidth = 91;
			nHeight = 51;
			nFocusW = 91;
			nFocusH = 51;
			pModel = new MeteModel();
			pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
//				, "../../data/t2-2007-2017-jan-144 and 240h-50(1Degree).txt", false
				, "../../data/t2-2007-2017-jan-144 and 240h-50(1Degree, single timestep).txt", false
//				, "../../data/t2-2007-2017-jan-144 and 240h-50(1Degree, single timestep,skip 3).txt", false
				, nWest, nEast, nSouth, nNorth
				, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth);
			/*
			pModel = new ContourBandDepthModel();
			if (bNewData) {
				pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
					//				, "../../data/t2-mod-ecmwf-20160105-00-72-216.txt", false
					//				, "../../data/t2-2007-2017-jan-120h-50.txt", false				// 这个区域小一点
					, "../../data/t2-2007-2017-jan-144 and 240h-50.txt", false
					, nWest, nEast, nSouth, nNorth
					, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth);
			}
			else {

				pModel->InitModel(50, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
					, "../../data/t2_mod_20080101-96h.dat", true
					, nWest, nEast, nSouth, nNorth
					, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth, g_bFilter);
			}
			*/
			break;
		case T2_Reanalysis:
			pModel = new ReanalysisModel();
			pModel->InitModel(1209, nWidth, nHeight, nFocusX, nFocusY, nFocusW, nFocusH
				, "../../data/t2_obs_1979-2017_1_china.txt", false
				, nWest, nEast, nSouth, nNorth
				, nFocusWest, nFocusEast, nFocusSouth, nFocusNorth);
			break;
		defaut:
			break;
		}

		return pModel;
	}
}

void MeteModel::generatePCAPoint(UncertaintyRegion& cluster, std::vector<DPoint3>& points) {
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

void MeteModel::clusterSpatialArea(UncertaintyRegion& cluster, std::vector<DPoint3>& points, ClusterResult& cr) {
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
	case MeteModel::bg_dipValue:
		buildTextureColorMap();
		break;
	case MeteModel::bg_cluster:
		//		return generateClusteredVarianceTexture_old();
		buildTextureClusteredVariance();
		break;
	case MeteModel::bg_varThreshold:
		buildTextureThresholdVariance();
		break;
	case MeteModel::bg_dipValueThreshold:
		buildTextureThresholdDipValue();
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

// cluster in each region
void MeteModel::regionCluster() {
	for (size_t i = 0, length = std::min((int)_vecRegions.size(), _nUncertaintyRegions); i < length; i++)
	{
		clusterSpatialArea(_vecRegions[i], _vecPCAPoints[i], _mxClusterResult[0][i]);
	}
}

void MeteModel::calculateSimilarity() {
	// 1.initialize the similarities to 0
	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		for (size_t j = 0; j < _nUncertaintyRegions; j++)
		{
			for (size_t k = 0; k <= _nClusters; k++)
			{
				_mxSimilarity[i][j][k] = 0;
			}
		}
	}
	// 2.counting each cluster
	for (size_t l = 0; l < _nEnsembleLen; l++)
	{
		for (size_t i = 0; i < _nUncertaintyRegions; i++)
		{
			for (size_t j = 0; j < i; j++) {
				if (_mxClusterResult[0][i]._arrLabels[l] == _mxClusterResult[0][j]._arrLabels[l])
				{
					_mxSimilarity[i][j][_mxClusterResult[0][i]._arrLabels[l]]++;
				}
			}
		}
	}

	// 3.counting the similarities for all the clusters
	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		for (size_t j = 0; j < _nUncertaintyRegions; j++)
		{
			for (size_t k = 0; k <= _nEnsembleLen; k++)
			{
				if (_mxClusterResult[0][i]._arrLabels[k] == _mxClusterResult[0][j]._arrLabels[k])
				{
					_mxSimilarity[i][j][_nClusters]++;
				}
			}
		}
	}
}

void MeteModel::SetFocusedRegion(int nRegion) {
	_nFocusedRegion = nRegion;

	alignClusters();
}


// align the cluster results
void MeteModel::alignClusters() {
	qDebug() << "Align Clusters: " << _nFocusedRegion;
	// first region, sort the clusters according to the number of their members
	_mxClusterResult[0][_nFocusedRegion].Sort();

	int nRegionSize = std::min((int)_vecRegions.size(), _nUncertaintyRegions);

	// 2.cluster in each uncertainty region
	for (size_t i = _nFocusedRegion+1; i < nRegionSize; i++)
	{
		// the following region match the region before
		_mxClusterResult[0][i].Match(_mxClusterResult[0][i - 1]);
	}
	for (int i = _nFocusedRegion-1; i >=0; i--)
	{
		// the following region match the region before
		_mxClusterResult[0][i].Match(_mxClusterResult[0][i + 1]);
	}
	// 3.align the results
	ClusterResult::Align(_mxClusterResult[0], _nUncertaintyRegions,_nFocusedRegion);

	// 4.reset labels for pca
	for (size_t i = 0; i < _nUncertaintyRegions; i++)
	{
		setLabelsForPCAPoints(_vecPCAPoints[i], _mxClusterResult[0][i]);
	}

	// 5.calculate the similarity between different uncertainty regions
	calculateSimilarity();
}