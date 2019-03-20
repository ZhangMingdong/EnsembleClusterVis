#include "DataField.h"

#include <QList>
#include <QDebug>
#include <QFile>
#include <QMessageBox>

#include <math.h>
#include <vector>
#include <algorithm>
#include <Windows.h>

#include "MyPCA.h"
#include "Switch.h"

DataField::DataField(int w, int h, int l)
{
	_nWidth = w;
	_nHeight = h;
	_nEnsembleLen = l;
	_nGrids = w * h;

	_pBuf = new double[w*h*l];
	_pDipBuf = new double[w*h];
	_gridVari = new double[w*h];
	_gridVV = new double[w*h];
	_gridMean = new double[w*h];
	_gridUMax = new double[w*h];
	_gridUMin = new double[w*h];
	_pSortedBuf = new double[w*h*l];
	_pResampledBuf = new double[_nGrids*_nResampleLen];



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
	delete[] _pBuf;
	delete[] _pSortedBuf;

	delete[] _pResampledBuf;


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

void DataField::SetData(int l, int bias, double dbValue) {
	_pBuf[l*_nGrids + bias] = dbValue;
}

void DataField::SetDipValue(int bias, double dbValue) {
	_pDipBuf[bias] = dbValue;

}

void DataField::SetData(int l, int y, int x, double dbValue) {
	_pBuf[l*_nGrids + y*_nWidth + x] = dbValue;
}

double DataField::GetData(int l, int bias) {
	return _pBuf[l*_nGrids + bias];
}

double DataField::GetData(int l, int r, int c) {
	return _pBuf[l*_nGrids + r* _nWidth + c];
}

double DataField::GetSortedData(int l, int bias) {
	return _pSortedBuf[l*_nGrids + bias];
}

void DataField::DoStatistic() {
	long t0,t1;

	if (g_bDebugDataField) { qDebug() << "DoStatistic:"; }
	// 0.calculate min and max
	if (g_bDebugDataField) { qDebug() << "calculate min and max"; t0 = GetTickCount(); }
	double dbMin = 100000;
	double dbMax = -100000;
	for (size_t i = 0; i < _nEnsembleLen*_nGrids; i++)
	{
		if (_pBuf[i] > dbMax) dbMax = _pBuf[i];
		if (_pBuf[i] < dbMin) dbMin = _pBuf[i];
	}
	if (g_bDebugDataField) { 
		qDebug() << "Min:" << dbMin; 
		qDebug() << "Max:" << dbMax; 
		t1 = GetTickCount();
		qDebug() << "============================TimeSpan:" << t1 - t0;
		t0 = t1;
	}

	// 1.sort buffer
	if (g_bDebugDataField) qDebug() << "sort buffer";
	sortBuf(_pBuf, _pSortedBuf);
	if (g_bDebugDataField) {
		t1 = GetTickCount();
		qDebug() << "============================TimeSpan:" << t1 - t0;
		t0 = t1;
	}
	// 2.resample buffer
	if (g_bDebugDataField) qDebug() << "resample buffer";
	if(g_bDomainResampleContours) resampleBuf(_pSortedBuf, _pResampledBuf, _nResampleLen);
	if (g_bDebugDataField) {
		t1 = GetTickCount();
		qDebug() << "============================TimeSpan:" << t1 - t0;
		t0 = t1;
	}

	double dbMaxMean = 0;
	double dbMaxVari = 0;
	double dbMinMean = 100000;
	double dbMinVari = 100000;

	// for each grid point
	for (int i = 0; i < _nGrids; i++)
	{
		// 1.calculate mean
		double fMean = 0;
		for (int j = 0; j < _nEnsembleLen; j++)
		{
			fMean += this->GetData(j, i);
		}
		fMean /= _nEnsembleLen;
		_gridMean[i] = fMean;
		
		// 2.calculate variance
		double fVariance = 0;
		for (int j = 0; j < _nEnsembleLen; j++)
		{
			double fBias = this->GetData(j, i) - fMean;
			fVariance += fBias*fBias;
		}
		_gridVari[i] = sqrt(fVariance / _nEnsembleLen);


		if (_gridMean[i] > dbMaxMean) dbMaxMean = _gridMean[i];
		if (_gridMean[i] < dbMinMean) dbMinMean = _gridMean[i];
		if (_gridVari[i] > dbMaxVari) dbMaxVari = _gridVari[i];
		if (_gridVari[i] < dbMinVari) dbMinVari = _gridVari[i];



		// 4.calculate max and min
		if (true)
		{
			// calculate max and min
			QList<double> list;
			for (int j = 0; j < _nEnsembleLen; j++)
			{
				list.append(this->GetData(j, i));
			}
			qSort(list);
			_gridUMin[i] = list[0];
			_gridUMax[i] = list[_nEnsembleLen - 1];
		}
		else {
			// calculate max and min
			double dbFactor = 1.0;
			_gridUMin[i] = _gridMean[i] - dbFactor*_gridVari[i];
			_gridUMax[i] = _gridMean[i] + dbFactor*_gridVari[i];
		}
	}

	qDebug() << "Mean:(" << dbMinMean << "," << dbMaxMean << ")";
	qDebug() << "Variance:(" << dbMinVari << "," << dbMaxVari << ")";

	if (g_bDebugDataField) {
		t1 = GetTickCount();
		qDebug() << "============================TimeSpan:" << t1 - t0;
		t0 = t1;
	}
	// 5.EOF
	if (g_bEOF)
		doEOF();

	if (g_bDebugDataField) {
		t1 = GetTickCount();
		qDebug() << "============================TimeSpan:" << t1 - t0;
		t0 = t1;
	}
}

void DataField::GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData) {
	// number of clusters
	int nClusters = listClusterLens.length();
	// generate data field according to the length of element in each cluster
	QList<int> listCurrentIndex;
	for (size_t i = 0; i < nClusters; i++)
	{
		arrData.push_back(new DataField(_nWidth, _nHeight, listClusterLens[i]));
		listCurrentIndex.push_back(0);
	}
	// set data for each new field
	for (size_t i = 0; i < _nEnsembleLen; i++)
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
	for (size_t i = 1; i < _nHeight-1; i++)
	{
		for (size_t j = 1; j < _nWidth-1; j++)
		{
			double dbVar = 0;
			for (int ii = -1; ii < 2; ii++)
			{
				for (int jj = -1; jj < 2; jj++) {
					dbVar += pVar[(i + ii)*_nWidth + j + jj];
				}
			}
			int nIndex = i* _nWidth + j;

			_pVarNew[nIndex] = dbVar / 9.0;
		}
	}

	// zero border
	for (size_t i = 0; i < _nHeight; i++)
	{
		_pVarNew[i*_nWidth + 0] = 0;
		_pVarNew[i*_nWidth + _nWidth -1] = 0;
	}
	for (size_t j = 1; j < _nWidth - 1; j++)
	{
		_pVarNew[j] = 0;
		_pVarNew[(_nHeight - 1)*_nWidth + j] = 0;
	}

}

void DataField::DoEOF_old() {
	qDebug() << "DoEOF";
	// 1.set parameter
	int mI = _nEnsembleLen;
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
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < mO; j++)
		{
			// reverse the order
			_gridEOF[mO-1-j][i] = arrOutput[i * mO + j];
		}
	}
	// 5.release the buffer
	delete[] arrInput;
	delete[] arrOutput;
}


void DataField::doEOF() {
	qDebug() << "doEOF";
	// 1.set parameter
	int mI = _nEnsembleLen;
	int mO = g_nEOFLen;
	int n = _nGrids;
	// 2.allocate input and output buffer
	double* arrInput = new double[mI*n];
	double* arrOutput = new double[mO*n];
	double* arrFactors = new double[mO*n];
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < mI; j++)
		{
			arrInput[i*mI + j] = GetData(j, i);
		}
	}

	// 3.pca
	MyPCA pca;
	pca.DoEOF(arrInput, arrOutput,arrFactors, n, mI, mO);
	// 4.generate points from the output
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < mO; j++)
		{
			// reverse the order
			_gridEOF[j][i] = arrOutput[i * mO + j];
		}
	}
	// 5.release the buffer
	delete[] arrInput;
	delete[] arrOutput;
	delete[] arrFactors;
}


void DataField::ReadDataFromText(QString filename) {
	QFile file(filename);

	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}

	QTextStream in(&file);

	// every ensemble member
	for (int i = 0; i < _nEnsembleLen; i++)
	{
		// skip first line
		QString line = in.readLine();
		// every grid
		for (int j = 0; j < _nGrids; j++)
		{
			QString line = in.readLine();
			SetData(i, j, line.toFloat());
		}
	}

	file.close();
}
