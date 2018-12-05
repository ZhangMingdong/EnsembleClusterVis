#include "FeatureSet.h"
#include "DataField.h"
#include "ContourGenerator.h"
#include "UnCertaintyArea.h"
#include "ContourStripGenerator.h"

#include <qDebug>


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

FeatureSet::FeatureSet(DataField* pData, double dbIsoValue, int nWidth, int nHeight, int nEnsembleLen, int nFocusX, int nFocusY, int nFocusW, int nFocusH):
	_pData(pData)
	,_dbIsoValue(dbIsoValue)
{
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nGrids = nWidth * nHeight;
	_nFocusX = nFocusX;
	_nFocusY = nFocusY;
	_nFocusW = nFocusW;
	_nFocusH = nFocusH;
	_nEnsembleLen = nEnsembleLen;


	_gridHalfMax = new double[_nGrids];
	_gridHalfMin = new double[_nGrids];
	_gridValidMax = new double[_nGrids];
	_gridValidMin = new double[_nGrids];
	_pSDF = new double[_nGrids*_nEnsembleLen];
	_pSortedSDF = new double[_nGrids*_nEnsembleLen];
	_pSet = new bool[_nGrids*_nEnsembleLen];
	_pSetBandDepth = new int[_nEnsembleLen];
	_pRegionType = new int[_nEnsembleLen];
}

FeatureSet::~FeatureSet()
{
	delete[] _gridHalfMax;
	delete[] _gridHalfMin;
	delete[] _gridValidMax;
	delete[] _gridValidMin;
	delete[] _pSDF;
	delete[] _pSortedSDF;
	delete[] _pSet;
	delete[] _pSetBandDepth;
	delete[] _pRegionType;

	for each (UnCertaintyArea* pArea in _listAreaValid)
		delete pArea;
	for each (UnCertaintyArea* pArea in _listAreaHalf)
		delete pArea;
	for each (UnCertaintyArea* pArea in _listUnionAreaE)
		delete pArea;

}

void FeatureSet::GenerateContours() {
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// spaghetti
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetLayer(i), contour, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContour.push_back(contour);

			// calculate sdf
			calculateSDF(_pData->GetData(i), GetSDF(i), _nWidth, _nHeight, _dbIsoValue, contour);
		}
		// sorted spaghetti
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetSortedLayer(i), contour, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSorted.push_back(contour);
		}


	}
	ContourGenerator::GetInstance()->Generate(_pData->GetUMin(), _listContourMinE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(_pData->GetUMax(), _listContourMaxE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(_pData->GetMean(), _listContourMeanE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetValidMin(), _listContourMinValid, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetValidMax(), _listContourMaxValid, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetHalfMin(), _listContourMinHalf, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetHalfMax(), _listContourMaxHalf, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetMedian(), _listContourMedianE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);


	generateContourImp(_listContourMinValid, _listContourMaxValid, _listAreaValid);
	generateContourImp(_listContourMinHalf, _listContourMaxHalf, _listAreaHalf);

}

void FeatureSet::BuildSortedSDF() {
	sortBuf(_pSDF, _pSortedSDF);

	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// contours from SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(GetSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSDF.push_back(contour);
		}
		// contours from sorted SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(GetSortedSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSortedSDF.push_back(contour);
		}
	}
}

void FeatureSet::sortBuf(const double* pS, double* pD) {
	for (int i = 0; i < _nGrids; i++) {
		for (int j = 0; j < _nEnsembleLen; j++) {
			double dbValue = pS[j*_nGrids + i];
			int k = 0;
			while (k < j && pD[k*_nGrids + i] < dbValue) k++;
			int l = j;
			while (l > k) {
				pD[l*_nGrids + i] = pD[(l - 1)*_nGrids + i];
				l--;
			}
			pD[k*_nGrids + i] = dbValue;
		}
	}
}

void FeatureSet::CalculateSet() {
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		for (int i = 0; i < _nGrids; i++) {
			_pSet[l*_nGrids + i] = (_pData->GetData(l,i) > _dbIsoValue);
		}
	}

	int nThreshold = 30;
	// caclulate sBand depth
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pSetBandDepth[i] = 0;
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (i == j) continue;
			for (size_t k = 0; k < _nEnsembleLen; k++)
			{
				if (i == k || j == k) continue;
				int nFailed = 0;
				for (int l = 0; l < _nGrids; l++)
				{
					if ((_pSet[i*_nGrids + l] && !_pSet[j*_nGrids + l] && !_pSet[k*_nGrids + l])
						|| (!_pSet[i*_nGrids + l] && _pSet[j*_nGrids + l] && _pSet[k*_nGrids + l])) {
						nFailed++;
					}
					if (nFailed > nThreshold) break;
				}
				if (nFailed < nThreshold) _pSetBandDepth[i]++;
			}
		}
	}
	/*
	// print depth value
	qDebug() << "Set Band Depth:";
	for (size_t i = 0; i < _nL; i++)
		qDebug() << _pSetBandDepth[i];

	qDebug() << "Set Band Depth Finished";
	*/

	// calculate region type

	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pRegionType[i] = 1;
	}
	int nOutliers = 0;
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		if (_pSetBandDepth[i] < _nOutlierThreshold) {
			_pRegionType[i] = 0;
			nOutliers++;
		}
	}
	int nValidLen = (_nEnsembleLen - nOutliers) / 2;
	for (int i = 0; i < nValidLen; i++)
	{
		int nMaxDepth = 0;
		int nIndex = -1;
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (_pRegionType[j] == 1 && _pSetBandDepth[j] > nMaxDepth)
			{
				nIndex = j;
				nMaxDepth = _pSetBandDepth[j];
			}
		}
		_pRegionType[nIndex] = 2;
	}

	qDebug() << "Region Types";
	for (size_t i = 0; i < _nEnsembleLen; i++)
		qDebug() << _pRegionType[i];

	qDebug() << "Region Types Finished";



	// migrate from calculate()

	// calculate valid max and min and mean
	for (int i = 0; i < _nGrids; i++) {
		std::vector<double> vecDataValid;
		std::vector<double> vecDataHalf;
		// 1.calculate mean
		for (int j = 0; j < _nEnsembleLen; j++)
		{
			if (_pRegionType[j])
			{
				vecDataValid.push_back(_pData->GetData(j, i));
				if (_pRegionType[j] == 2)
				{
					vecDataHalf.push_back(_pData->GetData(j, i));
				}
			}
		}


		std::sort(vecDataValid.begin(), vecDataValid.end());
		std::sort(vecDataHalf.begin(), vecDataHalf.end());
		_gridHalfMin[i] = vecDataHalf[0];
		_gridHalfMax[i] = vecDataHalf[vecDataHalf.size() - 1];
		_gridValidMin[i] = vecDataValid[0];
		_gridValidMax[i] = vecDataValid[vecDataValid.size() - 1];

	}

	// find median
	int nMaxDepth = 0;
	_nMedianIndex = -1;
	for (int i = 0; i < _nEnsembleLen; i++) {
		if (_pSetBandDepth[i] > nMaxDepth)
		{
			nMaxDepth = _pSetBandDepth[i];
			_nMedianIndex = i;
		}
	}

}

const double* FeatureSet::GetMedian() { return _pData->GetData(_nMedianIndex); }

void FeatureSet::generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas) {
	ContourStripGenerator generator;
	generator.Generate(areas, contourMin, contourMax, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
}

void FeatureSet::calculateSDF(const double* arrData, double* arrSDF, int nW, int nH, double isoValue, QList<ContourLine> contour) {
	for (size_t i = 0; i < nH; i++)
	{
		for (size_t j = 0; j < nW; j++)
		{
			double dbMinDis = nW * nH;
			for (size_t k = 0, lenK = contour.length(); k < lenK; k++)
			{
				for (size_t l = 0, lenL = contour[k]._listPt.length() - 1; l < lenL; l++)
				{
					double dbDis = PointToSegDist(j, i
						, contour[k]._listPt[l].x(), contour[k]._listPt[l].y()
						, contour[k]._listPt[l + 1].x(), contour[k]._listPt[l + 1].y());
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
