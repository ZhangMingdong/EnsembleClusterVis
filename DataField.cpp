#include "DataField.h"

#include <QList>
#include <math.h>

DataField::DataField(int w, int h, int l):_nW(w),_nH(h),_nL(l)
{
	_pBuf = new double[w*h*l];
	_gridVari = new double[w*h];
	_gridMean = new double[w*h];
	_gridUMax = new double[w*h];
	_gridUMin = new double[w*h];
}

DataField::~DataField()
{
	delete[] _gridVari;
	delete[] _gridMean;
	delete[] _gridUMax;
	delete[] _gridUMin;
	delete[] _pBuf;
}

const double* DataField::GetLayer(int l) {
	return _pBuf + l*_nW*_nH;
}

double* DataField::GetEditableLayer(int l) {
	return _pBuf + l*_nW*_nH;
}

const double* DataField::GetVari() { return _gridVari;}

const double* DataField::GetMean() { return _gridMean; }

const double* DataField::GetUMax() { return _gridUMax; }

const double* DataField::GetUMin() { return _gridUMin; }

void DataField::SetData(int l, int bias, double dbValue) {
	_pBuf[l*_nW*_nH + bias] = dbValue;
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

void DataField::DoStatistic() {
	int nLen = _nW*_nH;
	for (int i = 0; i < nLen; i++)
	{
		// calculate mean
		double fMean = 0;
		for (int j = 0; j < _nL; j++)
		{
			fMean += this->GetData(j, i);
		}
		fMean /= _nL;
		_gridMean[i] = fMean;
		// calculate variance
		double fVariance = 0;
		for (int j = 0; j < _nL; j++)
		{
			double fBias = this->GetData(j, i) - fMean;
			fVariance += fBias*fBias;
		}
		_gridVari[i] = sqrt(fVariance / _nL);
		if (false)
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
}

void DataField::GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData) {
	int nClusters = listClusterLens.length();
	QList<int> listCurrentIndex;
	for (size_t i = 0; i < nClusters; i++)
	{
		arrData.push_back(new DataField(_nW, _nH, listClusterLens[i]));
		listCurrentIndex.push_back(0);
	}
	int nLen = _nW*_nH;
	for (size_t i = 0; i < _nL; i++)
	{
		int nLabel = arrLabels[i];
		for (size_t j = 0; j < nLen; j++)
		{
			arrData[nLabel]->SetData(listCurrentIndex[nLabel], j, this->GetData(i, j));
		}
		listCurrentIndex[nLabel]++;
	}
	for (size_t i = 0; i < nClusters; i++)
	{
		arrData[i]->DoStatistic();
	}
}