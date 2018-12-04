#include "DataField.h"

#include <QList>
#include <QDebug>

#include <math.h>
#include <vector>
#include <algorithm>

#include "MyPCA.h"

DataField::DataField(int w, int h, int l):_nW(w),_nH(h),_nL(l)
{
	_nGrids = _nW * _nH;

	_pBuf = new double[w*h*l];
	_pDipBuf = new double[w*h];
	_gridVari = new double[w*h];
	_gridVV = new double[w*h];
	_gridMean = new double[w*h];
	_gridUMax = new double[w*h];
	_gridUMin = new double[w*h];
	_gridHalfMax = new double[w*h];
	_gridHalfMin = new double[w*h];
	_gridValidMax = new double[w*h];
	_gridValidMin = new double[w*h];
	_pSortedBuf = new double[w*h*l];
	_pSDF = new double[w*h*l];
	_pSortedSDF = new double[w*h*l];
	_pSet = new bool[w*h*l];
	_pSetBandDepth = new int[l];
	_pRegionType = new int[l];

	for (size_t i = 0; i < g_nEOFLen; i++)
	{
		_gridEOF[i]= new double[w*h];
	}
}

DataField::~DataField()
{
	delete[] _pDipBuf;
	delete[] _gridVari;
	delete[] _gridVV;
	delete[] _gridMean;
	delete[] _gridUMax;
	delete[] _gridUMin;
	delete[] _gridHalfMax;
	delete[] _gridHalfMin;
	delete[] _gridValidMax;
	delete[] _gridValidMin;
	delete[] _pBuf;
	delete[] _pSortedBuf;
	delete[] _pSDF;
	delete[] _pSortedSDF;
	delete[] _pSet;
	delete[] _pSetBandDepth;
	delete[] _pRegionType;
	for (size_t i = 0; i < _nSmooth; i++)
	{
		delete[] _gridVarSmooth[i];
	}
	for (size_t i = 0; i < g_nEOFLen; i++)
	{
		delete[] _gridEOF[i];
	}
}

const double* DataField::GetLayer(int l) {
	return _pBuf + l* _nGrids;
}

const double* DataField::GetSortedLayer(int l) {
	return _pSortedBuf + l * _nGrids;
}


double* DataField::GetEditableLayer(int l) {
	return _pBuf + l* _nGrids;
}

const double* DataField::GetVari(int nSmooth) {
	if (nSmooth==0)
	{
		return _gridVari;
	}
	else {
		for (size_t i= _nSmooth; i < nSmooth; i++)
		{
			smoothVar(i);
		}
		_nSmooth = nSmooth;
		return _gridVarSmooth[nSmooth-1];
	}
}

const double* DataField::GetEOF(int nSeq) {
	return _gridEOF[nSeq];
}

const double* DataField::GetDipValue() {
	return _pDipBuf;
}

const double* DataField::GetMean() { return _gridMean; }

const double* DataField::GetMedian() { return _pBuf+_nMedianIndex*_nGrids; }

const double* DataField::GetUMax() { return _gridUMax; }

const double* DataField::GetUMin() { return _gridUMin; }

void DataField::SetData(int l, int bias, double dbValue) {
	_pBuf[l*_nW*_nH + bias] = dbValue;
}

void DataField::SetDipValue(int bias, double dbValue) {
	_pDipBuf[bias] = dbValue;

}

void DataField::SetData(int l, int y, int x, double dbValue) {
	_pBuf[l*_nW*_nH + y*_nW + x] = dbValue;
}

double DataField::GetData(int l, int bias) {
	return _pBuf[l*_nW*_nH + bias];
}

double DataField::GetData(int l, int r, int c) {
	return _pBuf[l*_nW*_nH + r*_nW + c];
}

void DataField::sortBuf(const double* pS, double* pD) {
	for (int i = 0; i < _nGrids; i++) {
		for (int j = 0; j < _nL; j++) {
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

void DataField::buildSortedBuf() {
	sortBuf(_pBuf, _pSortedBuf);
}

void DataField::BuildSortedSDF() {
	sortBuf(_pSDF, _pSortedSDF);
}

void DataField::CalculateSet(double dbIsoValue) {
	for (int l = 0; l < _nL; l++)
	{
		for (int i = 0; i < _nGrids; i++) {
			_pSet[l*_nGrids + i] = (_pBuf[l*_nGrids + i] > dbIsoValue);
		}
	}

	int nThreshold =30;
	// caclulate sBand depth
	for (size_t i = 0; i < _nL; i++)
	{
		_pSetBandDepth[i] = 0;
		for (size_t j = 0; j < _nL; j++)
		{
			if (i == j) continue;
			for (size_t k = 0;  k< _nL; k++)
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
				if(nFailed < nThreshold) _pSetBandDepth[i]++;
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

	for (size_t i = 0; i < _nL; i++)
	{
		_pRegionType[i] = 1;
	}
	int nOutliers = 0;
	for (size_t i = 0; i < _nL; i++)
	{
		if (_pSetBandDepth[i] < _nOutlierThreshold) {
			_pRegionType[i] = 0;
			nOutliers++;
		}
	}
	int nValidLen = (_nL - nOutliers) / 2;
	for (int i = 0; i < nValidLen; i++)
	{
		int nMaxDepth = 0;
		int nIndex = -1;
		for (size_t j = 0; j < _nL; j++)
		{
			if (_pRegionType[j]==1 && _pSetBandDepth[j]>nMaxDepth)
			{
				nIndex = j;
				nMaxDepth = _pSetBandDepth[j];
			}
		}
		_pRegionType[nIndex] = 2;
	}

	qDebug() << "Region Types";
	for (size_t i = 0; i < _nL; i++)
		qDebug() << _pRegionType[i];

	qDebug() << "Region Types Finished";
}

void DataField::DoStatistic() {
	buildSortedBuf();

	// for each grid point
	for (int i = 0; i < _nGrids; i++)
	{
		std::vector<double> vecData;	// store the data for vv calculation
		// 1.calculate mean
		double fMean = 0;
		for (int j = 0; j < _nL; j++)
		{
			vecData.push_back(this->GetData(j, i));
			fMean += this->GetData(j, i);
		}
		fMean /= _nL;
		_gridMean[i] = fMean;
		// 2.calculate variance
		double fVariance = 0;
		for (int j = 0; j < _nL; j++)
		{
			double fBias = this->GetData(j, i) - fMean;
			fVariance += fBias*fBias;
		}
		_gridVari[i] = sqrt(fVariance / _nL);
		// 3.caluclate variance of variance
		// sort the data
		std::sort(vecData.begin(), vecData.end());
		// find the largest space
		const int nQueueLen = 10;
		double arrMax[nQueueLen];
		double arrMin[nQueueLen];
		for (size_t i = 0; i < nQueueLen; i++)
		{
			arrMax[i] = 0;
		}
		for (size_t i = 0; i < nQueueLen; i++)
		{
			arrMin[i] = 100000;
		}
		double dbRange = vecData[vecData.size() - 1] - vecData[0];
		for (size_t i = 1; i < vecData.size(); i++)
		{
			double dbSpace = vecData[i] - vecData[i - 1];
			// max array
			for (size_t j = 0; j < nQueueLen; j++)
			{
				if (dbSpace > arrMax[j]) {
					for (size_t k = nQueueLen-1; k > j; k--)
					{
						arrMax[k] = arrMax[k - 1];
					}
					arrMax[j] = dbSpace;
					break;
				}
			}
			// min array
			for (size_t j = 0; j < nQueueLen; j++)
			{
				if (dbSpace < arrMin[j]) {
					for (size_t k = nQueueLen - 1; k > j; k--)
					{
						arrMin[k] = arrMin[k - 1];
					}
					arrMin[j] = dbSpace;
					break;
				}
			}
		}
		double dbMax = 0;
		double dbMin = 0;
		for (size_t i = 0; i < nQueueLen; i++)
		{
			dbMax += arrMax[i];
			dbMin += arrMin[i];
		}
		double vv = dbMax / dbMin/10.0;
		_gridVV[i] = vv;
		/*
		// the codes for calculating variance of variance. old version
		double dbVV = 0;
		for (int j = 0; j < _nL; j++)
		{
			double fBias = this->GetData(j, i) - fMean;
			double fBB = fBias*fBias - _gridVari[i] * _gridVari[i];
			dbVV += fBB*fBB;
		}
		_gridVV[i] = sqrt(sqrt(dbVV / _nL));
		*/
		// 4.calculate max and min
		if (true)
		{
			// calculate max and min
			QList<double> list;
			for (int j = 0; j < _nL; j++)
			{
				list.append(this->GetData(j, i));
			}
			qSort(list);
			_gridUMin[i] = list[0];
			_gridUMax[i] = list[_nL - 1];
		}
		else {
			// calculate max and min
			double dbFactor = 1.0;
			_gridUMin[i] = _gridMean[i] - dbFactor*_gridVari[i];
			_gridUMax[i] = _gridMean[i] + dbFactor*_gridVari[i];
		}
	}

	// calculate valid max and min and mean
	for (int i = 0; i < _nGrids; i++) {
		std::vector<double> vecDataValid;
		std::vector<double> vecDataHalf;
		// 1.calculate mean
		for (int j = 0; j < _nL; j++)
		{
			if (_pRegionType[j])
			{
				vecDataValid.push_back(this->GetData(j, i));
				if (_pRegionType[j] ==2)
				{
					vecDataHalf.push_back(this->GetData(j, i));
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
	for (int i = 0; i < _nL; i++) {
		if (_pSetBandDepth[i] > nMaxDepth)
		{
			nMaxDepth = _pSetBandDepth[i];
			_nMedianIndex = i;
		}
	}
}

void DataField::GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData) {
	// number of clusters
	int nClusters = listClusterLens.length();
	// generate data field according to the length of element in each cluster
	QList<int> listCurrentIndex;
	for (size_t i = 0; i < nClusters; i++)
	{
		arrData.push_back(new DataField(_nW, _nH, listClusterLens[i]));
		listCurrentIndex.push_back(0);
	}
	// set data for each new field
	for (size_t i = 0; i < _nL; i++)
	{
		int nLabel = arrLabels[i];
		for (size_t j = 0; j < _nGrids; j++)
		{
			arrData[nLabel]->SetData(listCurrentIndex[nLabel], j, this->GetData(i, j));
		}
		listCurrentIndex[nLabel]++;
	}
	// do statistic for each new field
	for (size_t i = 0; i < nClusters; i++)
	{
		arrData[i]->DoStatistic();
	}
}

void DataField::smoothVar(int nSmooth) {
	const double* pVar = (nSmooth == 0) ? _gridVari : _gridVarSmooth[nSmooth - 1];

	double* _pVarNew= _gridVarSmooth[nSmooth] = new double[_nGrids];

	// smooth
	for (size_t i = 1; i < _nH-1; i++)
	{
		for (size_t j = 1; j < _nW-1; j++)
		{
			double dbVar = 0;
			for (int ii = -1; ii < 2; ii++)
			{
				for (int jj = -1; jj < 2; jj++) {
					dbVar += pVar[(i + ii)*_nW + j + jj];
				}
			}
			int nIndex = i*_nW + j;

			_pVarNew[nIndex] = dbVar / 9.0;
		}
	}

	// zero border
	for (size_t i = 0; i < _nH; i++)
	{
		_pVarNew[i*_nW + 0] = 0;
		_pVarNew[i*_nW + _nW-1] = 0;
	}
	for (size_t j = 1; j < _nW - 1; j++)
	{
		_pVarNew[j] = 0;
		_pVarNew[(_nH - 1)*_nW + j] = 0;
	}

}

void DataField::DoEOF() {
	// 1.set parameter
	int mI = _nL;
	int mO = g_nEOFLen;
	int n = _nGrids;
	// 2.allocate input and output buffer
	double* arrInput = new double[mI*n];
	double* arrOutput = new double[mO*n];
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < mI; j++) 
		{
			arrInput[i*mI + j] = GetData(j, i);
		}
	}
	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);
	// 4.generate points from the output
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < mO; j++)
		{
			_gridEOF[j][i] = arrOutput[i * mO + j];
		}
	}
	// 5.release the buffer
	delete[] arrInput;
	delete[] arrOutput;
}

