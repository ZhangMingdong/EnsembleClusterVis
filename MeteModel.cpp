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
#include "FeatureSet.h"

#include <fstream>
#include <iomanip>
#include <assert.h>
#include <sstream>

#include "def.h"

#include "DBSCANClustering.h"

#include "ArtificialModel.h"
#include "SpatialCluster.h"

#include "Switch.h"


using namespace std;

// 判断点在线段左侧还是右侧，网上找的公式，来不及仔细想了
bool checkLeft(double x, double y, double x1, double y1, double x2, double y2) {
	double dx1 = x2 - x1;
	double dy1 = y2 - y1;
	double dx0 = x1 - x;
	double dy0 = y1 - y;
	return dx0 * dy1 - dy0 * dx1 > 0;
}

MeteModel::MeteModel()
{
	_nWidth = g_nWidth;
	_nHeight = g_nHeight;
	_nWest = g_nWest;
	_nEast = g_nEast;
	_nNorth = g_nNorth;
	_nSouth = g_nSouth;

	_nGrids = _nWidth*_nHeight;

	_bufObs = new double[_nGrids];

	_gridErr = new double[_nGrids];
	//readObsData();

	for (size_t i = 0; i < 62; i++)
	{
		_arrTimeSteps[i] = 0;
	}
}

MeteModel::~MeteModel()
{
	if (_gridErr)
	{
		delete[]_gridErr;
	}
	if (_bufObs)
	{
		delete[]_bufObs;
	}
	//if (_pTimeStep)
	//{
	//	delete _pTimeStep;
	//}	
	for (size_t i = 0; i < 62; i++)
	{
		if (_arrTimeSteps[i]) delete _arrTimeSteps[i];
	}

	for (size_t i = 0; i < _listUnionAreaEG.length(); i++)
	{
		for each (UnCertaintyArea* pArea in _listUnionAreaEG[i])
			delete pArea;
	}
	if (_dataTexture) delete[] _dataTexture;




	ContourGenerator::Release();
}

void MeteModel::InitModel(int nEnsembleLen, int nWidth, int nHeight
	, QString strFile, bool bBinary, int nWest, int nEast, int nSouth, int nNorth) {
	// 0.record states variables
	_nEnsembleLen = nEnsembleLen;
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nGrids = nHeight*nWidth;


	_nWest = nWest;
	_nEast = nEast;
	_nSouth = nSouth;
	_nNorth = nNorth;

	_strFile = strFile;
	_bBinaryFile = bBinary;

	// 1.build data
	//_pTimeStep = new TimeStep();


	// 5.specializaed initialization
	initializeModel();
}

void MeteModel::readDipValue(char* strFileName) {

	QFile file(strFileName);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}


	QTextStream in(&file);
	// every grid
	for (int j = 0; j < _nGrids; j++)
	{
		QString line = in.readLine();
//		double dbDipValue = line.toFloat() * 100;
		double dbDipValue = line.toFloat() * 10;
		_pTimeStep->_pData->SetDipValue(j, dbDipValue);
		//qDebug() << j<<"\t" << dbDipValue;
	}

	file.close();
}

void MeteModel::readDipValueG(char* strFileName) {

	QFile file(strFileName);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}


	QTextStream in(&file);
	int nWidth = _nWidth - 1;
	for (size_t i = 0; i < _nHeight; i++)
	{
		for (size_t j = 0; j < nWidth; j++)
		{
			QString line = in.readLine();
			double dbDipValue = line.toFloat() * 10;
			_pTimeStep->_pData->SetDipValue(i*_nWidth+j, dbDipValue);
			if (j==0)
			{
				_pTimeStep->_pData->SetDipValue(i*_nWidth + _nWidth-1, dbDipValue);
			}
		}
	}

	file.close();
}

void MeteModel::readDataFromTextG() {
	QFile file(_strFile);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	QTextStream in(&file);


	int nWidth = _nWidth - 1;				// width of a scalar field in the file
	int nLen = nWidth*_nHeight;				// length of a scalar field in the file

	int nStep = 4;
	for (size_t t = 0; t < nStep; t++)
	{
		// every ensemble member
		for (int i = 0; i < _nEnsembleLen; i++)
		{
			QString line = in.readLine();
			// every grid
			for (int j = 0; j < nLen; j++)
			{
				QString line = in.readLine();
				int r = j / nWidth;
				int c = j % nWidth;
				int nIndex = r*_nWidth + c + 1;
				_pTimeStep->_pData->SetData(i, nIndex, line.toFloat());
				if (c == nWidth - 1)
				{
					int nIndex = r*_nWidth;
					_pTimeStep->_pData->SetData(i, nIndex, line.toFloat());
				}
			}
		}

	}

	file.close();
}

void MeteModel::buildTextureThresholdVariance() {
	const double* pData = _pTimeStep->_pData->GetVari(_nSmooth);
	_nThresholdedGridPoints = 0;
	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]>_dbVarThreshold)
			{
				_dataTexture[4 * nIndex + 0] = ColorMap::GetThresholdColorI(0);
				_dataTexture[4 * nIndex + 1] = ColorMap::GetThresholdColorI(1);
				_dataTexture[4 * nIndex + 2] = ColorMap::GetThresholdColorI(2);
				_dataTexture[4 * nIndex + 3] = (GLubyte)255;

				_nThresholdedGridPoints++;
			}
			else {
				_dataTexture[4 * nIndex + 0] = (GLubyte)255;
				_dataTexture[4 * nIndex + 1] = (GLubyte)255;
				_dataTexture[4 * nIndex + 2] = (GLubyte)255;
				_dataTexture[4 * nIndex + 3] = (GLubyte)255;

			}
		}
	}
}

void MeteModel::buildTextureSDF() {
	const double* pData = _pTimeStep->_listFeature[0]->GetSDF(_nMember? _nMember-1: _nMember);
	ColorMap* colormap = ColorMap::GetInstance(ColorMap::CP_EOF);


	double dbMax = -100000;
	double dbMin = 100000;	
	for (int i = 0; i < _nGrids; i++)
	{
		if (pData[i] > dbMax) dbMax = pData[i];
		if (pData[i] < dbMin) dbMin = pData[i];
	}
	colormap->SetRange(dbMin, dbMax);
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = color._rgb[0];
			_dataTexture[4 * nIndex + 1] = color._rgb[1];
			_dataTexture[4 * nIndex + 2] = color._rgb[2];
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
		}
	}
}

void MeteModel::buildTextureICDVX() {
	const double* pData = _pTimeStep->_listFeature[0]->GetICDVX(_nMember ? _nMember - 1 : _nMember);
	ColorMap* colormap = ColorMap::GetInstance(ColorMap::CP_EOF);


	double dbMax = -1000;
	double dbMin = 1000;	
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = color._rgb[0];
			_dataTexture[4 * nIndex + 1] = color._rgb[1];
			_dataTexture[4 * nIndex + 2] = color._rgb[2];
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
			if (pData[nIndex] > dbMax) dbMax = pData[nIndex];
			if (pData[nIndex] < dbMin) dbMin = pData[nIndex];
		}
	}
}

void MeteModel::buildTextureICDVY() {
	const double* pData = _pTimeStep->_listFeature[0]->GetICDVY(_nMember ? _nMember - 1 : _nMember);
	ColorMap* colormap = ColorMap::GetInstance(ColorMap::CP_EOF);


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = color._rgb[0];
			_dataTexture[4 * nIndex + 1] = color._rgb[1];
			_dataTexture[4 * nIndex + 2] = color._rgb[2];
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
			if (pData[nIndex] > dbMax) dbMax = pData[nIndex];
			if (pData[nIndex] < dbMin) dbMin = pData[nIndex];
		}
	}
}

void MeteModel::buildTextureICDVW() {
	const double* pData = _pTimeStep->_listFeature[0]->GetICDVW(_nMember ? _nMember - 1 : _nMember);
	ColorMap* colormap = ColorMap::GetInstance(ColorMap::CP_EOF);


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = color._rgb[0];
			_dataTexture[4 * nIndex + 1] = color._rgb[1];
			_dataTexture[4 * nIndex + 2] = color._rgb[2];
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
			if (pData[nIndex] > dbMax) dbMax = pData[nIndex];
			if (pData[nIndex] < dbMin) dbMin = pData[nIndex];
		}
	}
}

void MeteModel::buildTextureICD() {
	qDebug() << "buildTextureICD";
	int nDetaiScale = _pTimeStep->_listFeature[0]->GetDetailScale();
	_nTexW = (_nWidth - 1) * nDetaiScale + 1;
	_nTexH = (_nHeight - 1) * nDetaiScale + 1;
	const double* pData= _pTimeStep->_listFeature[0]->GetICD();	// the data


	for (int i = 0; i < _nTexH; i++) {
		for (int j = 0; j < _nTexW; j++) {

			int nIndex = i * _nTexW + j;
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = 100;
			_dataTexture[4 * nIndex + 1] = 100;
			_dataTexture[4 * nIndex + 2] = 200;
			_dataTexture[4 * nIndex + 3] = (GLubyte)(pData[nIndex]*255);
		}
	}
}

void MeteModel::buildTextureICD_LineKernel(){
	qDebug() << "line kernel" << endl;
	const double* pData = _pTimeStep->_listFeature[0]->GetICD_LineKernel();


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = 100;
			_dataTexture[4 * nIndex + 1] = 100;
			_dataTexture[4 * nIndex + 2] = 200;
			_dataTexture[4 * nIndex + 3] = (GLubyte)(pData[nIndex] * 255);
		}
	}
}

void MeteModel::buildTextureICDX() {
	qDebug() << "line kernel X" << endl;
	const double* pData = _pTimeStep->_listFeature[0]->GetICDX();


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = 100;
			_dataTexture[4 * nIndex + 1] = 100;
			_dataTexture[4 * nIndex + 2] = 200;
			_dataTexture[4 * nIndex + 3] = (GLubyte)(pData[nIndex] * 255);
		}
	}
}

void MeteModel::buildTextureICDY() {
	qDebug() << "line kernel Y" << endl;
	const double* pData = _pTimeStep->_listFeature[0]->GetICDY();


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = 100;
			_dataTexture[4 * nIndex + 1] = 100;
			_dataTexture[4 * nIndex + 2] = 200;
			_dataTexture[4 * nIndex + 3] = (GLubyte)(pData[nIndex] * 255);
		}
	}
}

void MeteModel::buildTextureICDZ() {
	qDebug() << "line kernel Z" << endl;
	const double* pData = _pTimeStep->_listFeature[0]->GetICDZ();


	double dbMax = -1000;
	double dbMin = 1000;	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i * _nWidth + j;
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = 100;
			_dataTexture[4 * nIndex + 1] = 100;
			_dataTexture[4 * nIndex + 2] = 200;
			_dataTexture[4 * nIndex + 3] = (GLubyte)(pData[nIndex] * 255);
		}
	}
}

void MeteModel::buildTextureThresholdDipValue() {
	double dbThreshold = 7;
	const double* pData = _pTimeStep->_pData->GetDipValue();
	_nThresholdedGridPoints = 0;
	// color map
	ColorMap* colormap = ColorMap::GetInstance();
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {
			int nIndex = i*_nWidth + j;
			if (pData[nIndex]<dbThreshold)
			{
				_dataTexture[4 * nIndex + 0] = ColorMap::GetThresholdColorI(0);
				_dataTexture[4 * nIndex + 1] = ColorMap::GetThresholdColorI(1);
				_dataTexture[4 * nIndex + 2] = ColorMap::GetThresholdColorI(2);
				_dataTexture[4 * nIndex + 3] = (GLubyte)255;

				_nThresholdedGridPoints++;
			}
			else {
				_dataTexture[4 * nIndex + 0] = (GLubyte)255;
				_dataTexture[4 * nIndex + 1] = (GLubyte)255;
				_dataTexture[4 * nIndex + 2] = (GLubyte)255;
				_dataTexture[4 * nIndex + 3] = (GLubyte)255;

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

	generateRegions();


	// 6.Generate texture
	// 6.1.initialize the texture
	for (size_t i = 0; i < _nGrids; i++)
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

void MeteModel::generateRegions() {
	_vecRegions.clear();
	/*
	// regions according to uncertainty
	const double* pData = _pData->GetVari(_nSmooth);
	SpatialCluster cluster;
	cluster.DoCluster(pData, _nWidth, _nHeight, _dbVarThreshold, _vecRegions);
	sort(_vecRegions.begin(), _vecRegions.end(), ClusterComparison);
	if (_vecRegions.size()<_nUncertaintyRegions)
		_nUncertaintyRegions = _vecRegions.size();

		*/

	/*
	// generate six region left up corner
	for (size_t i = 0; i < 6; i++)
	{
		UncertaintyRegion region;
		for (size_t j = 0; j < 11; j++)
		{
			for (size_t k = 0; k < 11; k++) {
				region._vecPoints.push_back(IPoint2(41 + j, 10 * i + k));
			}
		}
		region._nArea = 121;
		_vecRegions.push_back(region);
	}
	*/
	// hierarchical clustering

//	CLUSTER::Clustering* pClusterer = new CLUSTER::KMeansClustering();
	CLUSTER::Clustering* pClusterer = new CLUSTER::AHCClustering();
	// 1.parameters setting
	int nN = _nGrids;							// number of data items
	int nK = 6;							// number of clusters
	int nM = 1;								// dimension
	// 2.input data
	int* arrLabels = new int[nN];
	double* arrBuf = new double[nN * nM];
	for (size_t i = 0; i < nN; i++)
	{
		arrBuf[i] = _pTimeStep->_pData->GetMean()[i];
	}
	// 3.clustering
	pClusterer->DoCluster(nN, nM, nK, arrBuf, arrLabels);
	delete arrBuf;
	// 4.counting each cluster
	for (size_t i = 0; i < nK; i++)
	{
		UncertaintyRegion region;
		_vecRegions.push_back(region);
	}

	for (size_t i = 0; i < nN; i++) {
		int nW = i%_nWidth;
		int nH = i / _nWidth;
		_vecRegions[arrLabels[i]]._vecPoints.push_back(IPoint2(nH, nW));
	}
	/*
	for (size_t i = 0; i < _nHeight; i++)
	{
		for (size_t j = 0; j < _nWidth; j++)
		{
			int nIndex = i*_nWidth + j;
			_vecRegions[arrLabels[nIndex]]._vecPoints.push_back(IPoint2(i, j));
		}
	}
	*/
	delete[] arrLabels;
}

// generate texture of colormap of mean or variance
void MeteModel::buildTextureColorMap() {
	const double* pData;	// the data
	ColorMap* colormap;		// color map function

	switch (_bgFunction)
	{
	case MeteModel::bg_mean:
	case MeteModel::bg_Obs:
		//if (_type==MT_T2)
		//	colormap = ColorMap::GetInstance(ColorMap::CP_T2);
		//else if (_type == MT_AF)
		//	colormap = ColorMap::GetInstance();




		if (_bgFunction== MeteModel::bg_Obs)
		{
			pData = _bufObs;
		}
		//else if (_nEnsCluster)
		//{
		//	pData = _arrEnsClusterData[_nEnsCluster-1]->GetMean();
		//}
		else if (_nMember)
		{
			pData = _pTimeStep->_pData->GetLayer(_nMember - 1);
		}
		else {
			pData = _pTimeStep->_pData->GetMean();
		}
		{
			colormap = ColorMap::GetInstance(ColorMap::CP_EOF);
			double dbMin = 100000;
			double dbMax = -100000;
			for (size_t i = 0; i < _nGrids; i++)
			{
				if (pData[i] > dbMax) dbMax = pData[i];
				if (pData[i] < dbMin) dbMin = pData[i];
			}

			qDebug() << "min:" << dbMin;
			qDebug() << "max:" << dbMax;
			colormap->SetRange(dbMin, dbMax);
		}
		break;
	case MeteModel::bg_vari:
		pData = _pTimeStep->_pData->GetVari();
		colormap = ColorMap::GetInstance();
		break;
	case MeteModel::bg_vari_smooth:
		pData = _pTimeStep->_pData->GetVari(_nSmooth);
		colormap = ColorMap::GetInstance();
		break;
	case MeteModel::bg_dipValue:
		pData = _pTimeStep->_pData->GetDipValue();
		colormap = ColorMap::GetInstance();
		break;
	case MeteModel::bg_EOF:
		{

			pData = _pTimeStep->_pData->GetEOF(_nEOF - 1);
			colormap = ColorMap::GetInstance(ColorMap::CP_EOF);
			double dbMin = 1000;
			double dbMax = -1000;
			for (size_t i = 0; i < _nGrids; i++)
			{
				if (pData[i] > dbMax) dbMax = pData[i];
				if (pData[i] < dbMin) dbMin = pData[i];
			}
			qDebug() << "min:" << dbMin;
			qDebug() << "max:" << dbMax;
			colormap->SetRange(dbMin, dbMax);
		}

		//colormap = ColorMap::GetInstance(ColorMap::CP_T2);
		
		break;
	case MeteModel::bg_err:
		{
			const double* _pBiasBuf;


			//if (_nEnsCluster)
			//{
			//	_pBiasBuf = _arrEnsClusterData[_nEnsCluster - 1]->GetMean();
			//}
			//else 
			if (_nMember)
			{
				_pBiasBuf = _pTimeStep->_pData->GetLayer(_nMember - 1);
			}
			else {
				_pBiasBuf = _pTimeStep->_pData->GetMean();
			}
			double dbAccum= 0;
			for (size_t i = 0; i < _nGrids; i++)
			{
				int nRow = i / _nWidth;
				int nCol = i% _nWidth;

				_gridErr[i] = _pBiasBuf[i] - _bufObs[i];
			}
			// accumulate errors for the given region
			for (size_t i = 0; i < g_nClusterRegionR; i++)
			{
				for (size_t j = 0; j < g_nClusterRegionR; j++)
				{
					int nRow = g_nClusterRegionY + i;
					int nCol = g_nClusterRegionX + j;
					int nIndex = nRow*+_nWidth + nCol;
					dbAccum += abs(_gridErr[nIndex]);
				}

			}
			qDebug() << dbAccum;
			pData = _gridErr;
			colormap = ColorMap::GetInstance(ColorMap::CP_EOF);
		}
		break;
	default:
		break;
	}


	double dbMax = -1000;
	double dbMin = 1000;
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {

			int nIndex = i*_nWidth + j;
			MYGLColor color = colormap->GetColor(pData[nIndex]);
			// using transparency and the blue tunnel
			_dataTexture[4 * nIndex + 0] = color._rgb[0];
			_dataTexture[4 * nIndex + 1] = color._rgb[1];
			_dataTexture[4 * nIndex + 2] = color._rgb[2];
			_dataTexture[4 * nIndex + 3] = (GLubyte)255;
			if (pData[nIndex] > dbMax) dbMax = pData[nIndex];
			if (pData[nIndex] < dbMin) dbMin = pData[nIndex];
		}
	}
//	qDebug() << "Max: " << dbMax;
//	qDebug() << "Min: " << dbMin;
}

vector<double> MeteModel::GetVariance() {

	const double* pData = _pTimeStep->_pData->GetVari(_nSmooth);
	vector<double> vecVar;
	for (int i = 0; i < _nHeight; i++) {
		for (int j = 0; j < _nWidth; j++) {
			int nIndex = i*_nWidth + j;
			vecVar.push_back(pData[nIndex]);
		}
	}
	std::sort(vecVar.begin(), vecVar.end());
	return vecVar;
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
				_pTimeStep->_pData->SetData(l, i, j, *f++);
				// 只取整度，过滤0.5度
//				if (_bFilter&&j < _nWidth - 1) f++;
			}
			// 只取整度，过滤0.5度
//			if (_bFilter&&i < _nHeight - 1) f += (_nWidth * 2 - 1);
		}
	}
}

void MeteModel::initializeModel() {
	_nTime = g_nTimeStart;
	while (_nTime <= g_nTimeEnd) {
		updateTimeStep();
		_nTime += 6;
		//break;
	}
}

void MeteModel::updateTimeStep() {
	int nTimeIndex = _nTime / 6;
	qDebug() << "nTimeIndex" << nTimeIndex;
	if (_arrTimeSteps[nTimeIndex]) {
		_pTimeStep = _arrTimeSteps[nTimeIndex];
	}
	else {
		_pTimeStep = _arrTimeSteps[nTimeIndex] = new TimeStep();
		_pTimeStep->Init(_nWidth, _nHeight, _nEnsembleLen);

		// 1.Read data
		_strFile = g_strPath + QString("-mb") + QString::number(_nTime) + QString(".txt");
		_pTimeStep->_pData->ReadDataFromText(_strFile);

		// 2.statistic
		_pTimeStep->_pData->DoStatistic();

		// 2.Set isovalues
		{
			QList<double> listIsoValue;

			//listIsoValue.append(273.16 - 20);
			//listIsoValue.append(273.16 - 15);
			//listIsoValue.append(273.16 - 10);
			//listIsoValue.append(273.16 - 5);
			//listIsoValue.append(273.16);
			//listIsoValue.append(273.16 + 5);
			//listIsoValue.append(273.16 + 10);
			//listIsoValue.append(273.16 + 15);
			//listIsoValue.append(273.16 + 20);

			//listIsoValue.append(5300);
			//listIsoValue.append(5350);
			//listIsoValue.append(5400);
			//listIsoValue.append(5450);
			//listIsoValue.append(5500);
			//listIsoValue.append(5550);
			//listIsoValue.append(5600);

			//listIsoValue.append(5300);
			//listIsoValue.append(5400);
			//listIsoValue.append(5500);
			//listIsoValue.append(5600);
			//listIsoValue.append(5700);
			//listIsoValue.append(5800);
			//listIsoValue.append(5900);


			//listIsoValue.append(5300);
			//listIsoValue.append(5500);
			//listIsoValue.append(5700);
			//listIsoValue.append(5900);


			//listIsoValue.append(5500);
			//listIsoValue.append(5600);
			//listIsoValue.append(5700);


			/*
			for (int i=-5;i<3;i++)
			{

				listIsoValue.append(5880+i*80);
			}
			*/

			//listIsoValue.append(5580);
			//listIsoValue.append(5600);
			//listIsoValue.append(5620);


			//listIsoValue.append(5500);
			//listIsoValue.append(5600);
			//listIsoValue.append(5700);

			//listIsoValue.append(5880);


			//listIsoValue.append(5500);
			//listIsoValue.append(5540);
			//listIsoValue.append(5580);
			//listIsoValue.append(5620);
			//listIsoValue.append(5660);
			//listIsoValue.append(5700);

			int nBase = 5620;
			int nStep = 30;
			for (int  i = -3; i < 3; i++)
			{
				listIsoValue.append(nBase + i * nStep);
			}
			//listIsoValue.append(nBase);

			SetIsoValues(listIsoValue);
		}

		// 4.generate feature;
		for each (double isoValue in _listIsoValues)
		{
			_pTimeStep->_listFeature.append(new FeatureSet(_pTimeStep->_pData, isoValue, _nWidth, _nHeight, _nEnsembleLen));
		}

	}
}

void MeteModel::Brush(int nLeft, int nRight, int nTop, int nBottom) {
	/*
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
	*/
}

QList<QList<ContourLine>> MeteModel::GetContourBrushed()
{
	return _listContourBrushed;
}

QList<QList<ContourLine>> MeteModel::GetContourNotBrushed()
{
	return _listContourNotBrushed;
}

QList<QList<ContourLine>> MeteModel::GetContour(int isoIndex)
{
	QList<QList<ContourLine>> result;

	QList<QList<ContourLine>>& listContour = _pTimeStep->_listFeature[isoIndex]->GetContours();
	if (_bgFunction==bg_EOF)
	{
		return _listContourEOF;
	}
	else 
	{
		if (_nEnsCluster)
		{
			for (size_t i = 0; i < _nEnsembleLen; i++)
			{
				QList<ContourLine> emptyContour;
				if (GetLabel(i) == _nEnsCluster - 1)
					result.push_back(listContour[i]);
				else
					result.push_back(emptyContour);
			}
			return result;
		}
		else if (_nMember)
		{
			//return _listMemberContour[_nMember - 1];
			result.push_back(listContour[_nMember - 1]);
			return result;
		}
		else return listContour;
	}
}

void MeteModel::GetMerge(int l, int& nSource, int& nTarget) {
	nSource = _pTimeStep->_listFeature[0]->GetMergeSource(l);
	nTarget = _pTimeStep->_listFeature[0]->GetMergeTarget(l);
}

int MeteModel::GetLabel(int l) {
	return _pTimeStep->_listFeature[0]->GetLabel(l);
}

double MeteModel::GetPC(int l,int nIndex) {
	return _pTimeStep->_listFeature[0]->GetPC(l,nIndex);
}

int MeteModel::GetClusters() { 
	return _pTimeStep->_listFeature[0]->GetClusters();
}

QList<QList<ContourLine>> MeteModel::GetContourSmooth(int isoIndex)
{
	QList<QList<ContourLine>>& listContour = _pTimeStep->_listFeature[isoIndex]->GetContoursSmooth();
	if (_nMember)
	{
		//return _listMemberContour[_nMember - 1];
		QList<QList<ContourLine>> result;
		result.push_back(listContour[_nMember - 1]);
		return result;

	}
	else return listContour;
}

QList<QList<ContourLine>> MeteModel::GetContourOutlier(int isoIndex)
{
	QList<QList<ContourLine>> listContour = _pTimeStep->_listFeature[isoIndex]->GetContours();
	QList<QList<ContourLine>> result;
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		if (_pTimeStep->_listFeature[isoIndex]->GetMemberType(i)==0)
		{
			result.push_back(listContour[i]);
		}

	}
	return result;
}

// add item from source to target
void addContour(const QList<QList<ContourLine>>& source, QList<QList<ContourLine>>& target, int nIndex0, int nIndex1, int nLevel) {
	if (nIndex1 > nIndex0) {
		int nMedian = (nIndex1 + nIndex0) / 2;
		target.push_back(source[nMedian]);
		if (nLevel > 0) {
			addContour(source, target, nIndex0, nMedian, nLevel - 1);
			addContour(source, target, nMedian+1, nIndex1, nLevel - 1);
		}
	}
}

QList<QList<ContourLine>> MeteModel::GetContourSorted(int isoIndex)
{
	QList<QList<ContourLine>> listResult;
	QList<QList<ContourLine>> contours= _pTimeStep->_listFeature[isoIndex]->GetContourSorted();
	addContour(contours, listResult, 0, contours.size(), _arrContourLevel[isoIndex]);
	return listResult;
}

QList<QList<ContourLine>> MeteModel::GetContourSortedSDF(int isoIndex)
{
	QList<QList<ContourLine>> listResult;
	QList<QList<ContourLine>> contours = _pTimeStep->_listFeature[isoIndex]->GetContourSortedSDF();
	addContour(contours, listResult, 0, contours.size(), _arrContourLevel[isoIndex]);
	return listResult;
}

QList<QList<ContourLine>> MeteModel::GetContourResampled(int isoIndex)
{
	QList<QList<ContourLine>> listResult;
	QList<QList<ContourLine>> contours = _pTimeStep->_listFeature[isoIndex]->GetContourResampled();
	addContour(contours, listResult, 0, contours.size(), _arrContourLevel[isoIndex]);
	return listResult;
}

QList<QList<ContourLine>> MeteModel::GetContourResampled_C()
{
	QList<QList<ContourLine>> listResult;
	QList<QList<ContourLine>> contours = _pTimeStep->_listFeature[0]->GetContourResampled_C();


	int nClusters = _pTimeStep->_listFeature[0]->GetClusters();
	if (_arrContourLevel[0] ==0)
	{
		for (int i = 0; i < nClusters; i++) {
			listResult.push_back(contours[i * 11 + 6]);
		}
	}
	else if (_arrContourLevel[0] == 1)
	{
		for (int i = 0; i < nClusters; i++) {
			listResult.push_back(contours[i * 11 + 4]);
			listResult.push_back(contours[i * 11 + 8]);
		}
	}
	else if (_arrContourLevel[0] == 2)
	{
		for (int i = 0; i < nClusters; i++) {
			listResult.push_back(contours[i * 11 + 3]);
			listResult.push_back(contours[i * 11 + 6]);
			listResult.push_back(contours[i * 11 + 9]);
		}
	}
	return listResult;
}
QList<QList<ContourLine>> MeteModel::GetContourDomainResampled(int isoIndex)
{
	QList<QList<ContourLine>> listResult;
	QList<QList<ContourLine>> contours = _pTimeStep->_listFeature[isoIndex]->GetContourDomainResampled();
	addContour(contours, listResult, 0, contours.size(), _arrContourLevel[isoIndex]);
	return listResult;
}

QList<QList<ContourLine>> MeteModel::GetContourSDF(int isoIndex)
{
	QList<QList<ContourLine>> contours = _pTimeStep->_listFeature[isoIndex]->GetContourSDF();
	if (_nMember)
	{
		//return _listMemberContour[_nMember - 1];
		QList<QList<ContourLine>> result;
		result.push_back(contours[_nMember - 1]);
		return result;
	}
	else return contours;
}

MeteModel* MeteModel::CreateModel(bool bA) {
	// 1.Create model instance
	//MeteModel* pModel = new MeteModel();
	MeteModel* pModel = bA ? new ArtificialModel() : new MeteModel();



	// 3.Initialize model
	int nWidth = g_nWidth;
	int nHeight		= g_nHeight;
	int nFocusX		= g_nFocusX;
	int nFocusY		= g_nFocusY;
	int nFocusW		= g_nFocusW;
	int nFocusH		= g_nFocusH;
	int nWest		= g_nWest;
	int nEast		= g_nEast;
	int nNorth		= g_nNorth;
	int nSouth		= g_nSouth;
	int nFocusWest	= g_nWest;
	int nFocusEast	= g_nEast;
	int nFocusNorth	= g_nNorth;
	int nFocusSouth	= g_nSouth;

	switch (g_usedModel)
	{
	case PRE_CMA:		
		pModel->InitModel(14, nWidth, nHeight, "../../data/data10/pre-mod-cma-20160802-00-96.txt"); break;
	case PRE_CPTEC:
		pModel->InitModel(14, nWidth, nHeight, "../../data/data10/pre-mod-cptec-20160802-00-96.txt"); break;
	case PRE_ECCC:
		pModel->InitModel(20, nWidth, nHeight, "../../data/data10/pre-mod-eccc-20160802-00-96.txt"); break;
	case PRE_ECMWF:
		pModel->InitModel(50, nWidth, nHeight, "../../data/data10/pre-mod-ecmwf-20160802-00-96.txt"); break;
	case PRE_JMA:
		pModel->InitModel(26, nWidth, nHeight, "../../data/data10/pre-mod-jma-20160802-00-96.txt"); break;
	case PRE_KMA:
		pModel->InitModel(24, nWidth, nHeight, "../../data/data10/pre-mod-kma-20160802-00-96.txt"); break;
	case PRE_NCEP:
		pModel->InitModel(20, nWidth, nHeight, "../../data/data10/pre-mod-ncep-20160802-00-96.txt"); break;
	case T2_ECMWF:
		pModel->InitModel(50, nWidth, nHeight, g_strFileName,false, nWest, nEast, nSouth, nNorth);
		break;
	case PRE_ECMWF_2017:
		pModel->InitModel(50, nWidth, nHeight, g_strFileName);		
		break;
	defaut:
		break;
	}

	return pModel;
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
			arrInput[i*mI + j] = _pTimeStep->_pData->GetData(i, x, y);
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

void MeteModel::SetEOF(int nEOF) {
	if (nEOF == _nEOF) return;
	_nEOF = nEOF;

	_listContourEOF.clear();
	for (int i = -10; i <= 10; i += 20)
	{
		QList<ContourLine> contour;
		ContourGenerator::GetInstance()->Generate(_pTimeStep->_pData->GetEOF(_nEOF-1), contour, i, _nWidth, _nHeight);
		_listContourEOF.push_back(contour);
	}
	regenerateTexture();
}

void MeteModel::SetMember(int nMember) {
	_nMember = nMember;
	regenerateTexture();
}

void MeteModel::SetEnsCluster(int nEnsCluster) {
	_nEnsCluster = nEnsCluster;
	regenerateTexture();
}

void MeteModel::SetEnsClusterLen(int nEnsClusterLen) {
	_pTimeStep->_listFeature[0]->SetClustersLen(nEnsClusterLen);
	regenerateTexture();
}

void MeteModel::regenerateTexture() {
	_nTexW = _nWidth;
	_nTexH = _nHeight;
	if (!_dataTexture)
	{
		double dbScale = 100;
		if (_pTimeStep->_listFeature.size())
			dbScale = _pTimeStep->_listFeature[0]->GetDetailScale();
		_dataTexture = new GLubyte[4 * _nGrids*dbScale*dbScale];
	}

	switch (_bgFunction)
	{
	case MeteModel::bg_mean:
	case MeteModel::bg_Obs:
	case MeteModel::bg_vari:
	case MeteModel::bg_vari_smooth:
	case MeteModel::bg_dipValue:
	case MeteModel::bg_EOF:
	case MeteModel::bg_err:
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
	case MeteModel::bg_SDF:
		buildTextureSDF();
		break;
	case MeteModel::bg_ICDVX:
		buildTextureICDVX();
		break;
	case MeteModel::bg_ICDVY:
		buildTextureICDVY();
		break;
	case MeteModel::bg_ICDVW:
		buildTextureICDVW();
		break;
	case MeteModel::bg_IsoContourDensity:
		buildTextureICD();
		break;
	case MeteModel::bg_LineKernel:
		buildTextureICD_LineKernel();
		break;
	case MeteModel::bg_LineKernelX:
		buildTextureICDX();
		break;
	case MeteModel::bg_LineKernelY:
		buildTextureICDY();
		break;
	case MeteModel::bg_LineKernelZ:
		buildTextureICDZ();
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

GLubyte* MeteModel::GenerateTexture() {
	if (!_dataTexture) {
		regenerateTexture();
	}
	return _dataTexture;
}

void MeteModel::doEnsCluster() {
	CLUSTER::Clustering* pClusterer = new CLUSTER::KMeansClustering();
	// 1.parameters setting
	int nN = _nEnsembleLen;					// number of data items
	int nK = g_nEnsClusterLen;				// number of clusters
	/*
	int nM = _nLen;							// dimension
	// 2.input data
	double* arrBuf = new double[nN * nM];
	for (size_t i = 0; i < nN; i++)
	{
		for (size_t j = 0; j < nM; j++)
		{
			arrBuf[i * nM + j] = _pData->GetData(i, j);
		}
	}
	*/
	int nM = g_nClusterRegionR*g_nClusterRegionR;							// dimension
	// 2.input data
	double* arrBuf = new double[nN * nM];
	for (size_t i = 0; i < nN; i++)
	{
		for (size_t j = 0; j < g_nClusterRegionR; j++)
		{
			for (size_t k = 0; k < g_nClusterRegionR; k++) {
				arrBuf[i * nM + j * g_nClusterRegionR + k] = _pTimeStep->_pData->GetData(i,g_nClusterRegionY+j,g_nClusterRegionX+k);
			}
		}
	}
	// 3.clustering
	pClusterer->DoCluster(nN, nM, nK, arrBuf, _arrLabels);
	delete arrBuf;
	// 4.counting each cluster
	int arrLens[g_nEnsClusterLen];
	for (size_t i = 0; i < g_nEnsClusterLen; i++)
	{
		arrLens[i] = 0;
	}

	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		arrLens[_arrLabels[i]]++;
	}


	QList<int> listLen;
	for (size_t i = 0; i < g_nEnsClusterLen; i++)
	{
		listLen.append(arrLens[i]);
	}
	// 5.postprocess the clustered data
	_pTimeStep->_pData->GenerateClusteredData(listLen, _arrLabels, _arrEnsClusterData);

	// 6.output labels
	for (size_t i = 0; i < g_nEnsClusterLen; i++)
	{
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (_arrLabels[j] == i) qDebug() << j;
		}
		qDebug()<<"\n";
	}
}

void MeteModel::readObsData() {
	int nBias = g_nBiasY * 31 + g_nBiasD;
	QFile file(g_strObsFileName);

	if (!file.open(QIODevice::ReadOnly)) {
		for (int j = 0; j < _nGrids; j++)
		{
			_bufObs[j] = 0;
		}
		QMessageBox::information(0, "error", file.errorString());
		return;
	}	

	QTextStream in(&file);

	// every ensemble member
	for (int i = 0; i < nBias; i++)
	{
		// skip first line
		QString line = in.readLine();
		// every grid
		for (int j = 0; j < _nGrids; j++)
		{
			QString line = in.readLine();
		}
	}
	// skip first line
	QString line = in.readLine();
	/*
	// every grid
	for (int j = 0; j < _nLen; j++)
	{
		QString line = in.readLine();
		_bufObs[j] = line.toFloat();
	}
	*/
	// observation is upside down
	// every grid
	for (size_t i = 0; i < _nHeight; i++)
	{
		for (size_t j = 0; j < _nWidth; j++) {
			QString line = in.readLine();
			int index = (_nHeight - i - 1)*_nWidth + j;
			_bufObs[index] = line.toFloat();
		}
	}
	file.close();
}

QList<ContourLine> MeteModel::GetContourMin(int isoIndex)  { return _pTimeStep->_listFeature[isoIndex]->GetContourMin() ; }
QList<ContourLine> MeteModel::GetContourMax(int isoIndex)  { return _pTimeStep->_listFeature[isoIndex]->GetContourMax() ; }
QList<ContourLine> MeteModel::GetContourMean(int isoIndex) { return _pTimeStep->_listFeature[isoIndex]->GetContourMean(); }
QList<ContourLine> MeteModel::GetContourMedian(int isoIndex) { return _pTimeStep->_listFeature[isoIndex]->GetContourMedian(); }

QList<UnCertaintyArea*> MeteModel::GetUncertaintyAreaValid(int nTimeIndex,int isoIndex) { return _arrTimeSteps[nTimeIndex]->_listFeature[isoIndex]->GetUncertaintyAreaValid(); }
QList<UnCertaintyArea*> MeteModel::GetUncertaintyAreaHalf(int nTimeIndex, int isoIndex) { return _arrTimeSteps[nTimeIndex]->_listFeature[isoIndex]->GetUncertaintyAreaHalf(); }
QList<UnCertaintyArea*> MeteModel::GetUncertaintyArea(int nTimeIndex, int isoIndex) { return _arrTimeSteps[nTimeIndex]->_listFeature[isoIndex]->GetUncertaintyArea(); }

void MeteModel::SetIsoValues(QList<double> listIsoValues) {
	_listIsoValues = listIsoValues;
}

void MeteModel::updateTimeStep(int nTS) {
	qDebug() << "MeteModel::updateTimeStep" << nTS;
	_nTime = nTS;
	updateTimeStep();
	UpdateView();
}

int MeteModel::GetCurrentTimeIndex() {
	return (_nTime - g_nTimeStart) / 6;
}

void MeteModel::updateContourLevel0(int nLevel) { _arrContourLevel[0] = nLevel; UpdateView();
qDebug() << "updateContourLevel0";
}
void MeteModel::updateContourLevel1(int nLevel) { _arrContourLevel[1] = nLevel; UpdateView();}
void MeteModel::updateContourLevel2(int nLevel) { _arrContourLevel[2] = nLevel; UpdateView();}
void MeteModel::updateContourLevel3(int nLevel) { _arrContourLevel[3] = nLevel; UpdateView();}
void MeteModel::updateContourLevel4(int nLevel) { _arrContourLevel[4] = nLevel; UpdateView();}
void MeteModel::updateContourLevel5(int nLevel) { _arrContourLevel[5] = nLevel; UpdateView();}
void MeteModel::updateContourLevel6(int nLevel) { _arrContourLevel[6] = nLevel; UpdateView();}
void MeteModel::updateContourLevel7(int nLevel) { _arrContourLevel[7] = nLevel; UpdateView();}
void MeteModel::updateContourLevel8(int nLevel) { _arrContourLevel[8] = nLevel; UpdateView();}
void MeteModel::updateContourLevel9(int nLevel) { _arrContourLevel[9] = nLevel; UpdateView();}