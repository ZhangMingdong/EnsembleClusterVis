#include "ReanalysisModel.h"
#include "def.h"
#include "DataField.h"

#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QMessageBox>

#include <iostream>

using namespace std;

ReanalysisModel::ReanalysisModel():_pRelationH(NULL), _pRelationV(NULL)
{
}


ReanalysisModel::~ReanalysisModel()
{
	if (_pRelationH) delete _pRelationH;
	if (_pRelationV) delete _pRelationV;
}

void ReanalysisModel::InitModel(int nEnsembleLen, int nWidth, int nHeight, int nFocusX, int nFocusY, int nFocusW, int nFocusH
	, QString strFile, bool bBinary, int nWest, int nEast, int nSouth, int nNorth
	, int nFocusWest, int nFocusEast, int nFocusSouth, int nFocusNorth, bool bFilter) {
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


	_pRelationH = new DataField(_nWidth-1, _nHeight, _nEnsembleLen);
	_pRelationV = new DataField(_nWidth, _nHeight-1, _nEnsembleLen);

	readDataFromText();

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

	// 6.calculate relation
	calculateRelation();

	return;
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
	}

	generateContourImp(_listContourMinE, _listContourMaxE, _listUnionAreaE);

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



void ReanalysisModel::readDataFromText() {

	int nTimeStep = g_nTimeStep;

	QFile file(_strFile);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	int nCount = 0;
	QTextStream in(&file);
	int tt = 0;
	while (!in.atEnd()) {
		for (size_t t = 0; t <= nTimeStep; t++)
		{
			for (int i = 0; i < _nEnsembleLen; i++)
			{
				QString line = in.readLine();
				for (int j = 0; j < _nLen; j++)
				{
					QString line = in.readLine();
					nCount++;
					int r = j / _nWidth;
					int c = j%_nWidth;
					r = _nHeight - r - 1;
					if (t == nTimeStep)
						_pData->SetData(i, r*_nWidth+c, line.toFloat());					
				}
			}
		}
		tt++;
		if (tt == 2) break;		// use the second data
	}

	file.close();
}

void ReanalysisModel::calculateRelation() {
	for (size_t l = 0; l < _nEnsembleLen; l++)
	{
		for (size_t i = 0; i < _nHeight-1; i++)
		{
			for (size_t j = 0; j < _nWidth-1; j++)
			{
				int nIndex = i*_nWidth + j;
				int nIndexR = nIndex + 1;
				int nIndexB = nIndex + _nWidth;
				double dbValue = _pData->GetData(l, nIndex);
				double dbValueR = _pData->GetData(l, nIndexR);
				double dbValueB = _pData->GetData(l, nIndexB);
				double dbRelationH = dbValueR - dbValue;
				double dbRelationV = dbValueB - dbValue;
				_pRelationH->SetData(l, nIndex, dbRelationH);
				_pRelationV->SetData(l, nIndex, dbRelationV);
			}
		}
	}
}