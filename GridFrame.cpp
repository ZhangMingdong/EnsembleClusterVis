#include "GridFrame.h"
#include "MathTypes.hpp"

#include "Switch.h"
#include <QDebug>
#include <Windows.h>


GridFrame::GridFrame()
{
}


GridFrame::~GridFrame()
{
}

void GridFrame::sortBuf(const double* pS, double* pD) {
	for (int i = 0; i < _nGrids; i++) {								// i: index of grid points
		for (int j = 0; j < _nEnsembleLen; j++) {					// j: index of ensemble member
			double dbValue = pS[j*_nGrids + i];						// dbValue: get value of ensemble j
			int k = 0;												// k: finds the first value bigger than dbValue at k
			while (k < j && pD[k*_nGrids + i] < dbValue) k++;
			int l = j;
			while (l > k) {											// move the value behind k one step further
				pD[l*_nGrids + i] = pD[(l - 1)*_nGrids + i];
				l--;
			}
			pD[k*_nGrids + i] = dbValue;							// set new value at k
		}
	}
	// vector version
	/*
		for (int i = 0; i < _nGrids; i++) {
		std::vector<double> vecValues;
		for (int j = 0; j < _nEnsembleLen; j++) {
			vecValues.push_back(pS[j*_nGrids + i]);
		}
		std::sort(vecValues.begin(), vecValues.end());
		for (int j = 0; j < _nEnsembleLen; j++) {
			pD[j*_nGrids + i] = vecValues[j];
		}
	}
	*/
}


const double c_dbK = 1.0 / sqrt(PI2d);

inline double KernelFun(double para) {
	return c_dbK * exp(-para * para / 2.0);
}
void resampleBasedOnKDE_FixStep(double* dbInput, double* dbOutput, int nInputLen, int nOutputLen) {
	double dbStep = 0.01;								// resample step
	double dbMin = dbInput[0] - 5;						// resample min border
	double dbMax = dbInput[nInputLen - 1] + 5;			// resample max border
	int resampleGridLen = (dbMax - dbMin) / dbStep;		// resample grid len
	if (resampleGridLen<10)
	{
		double dbSpan = (dbMax - dbMin) / (nOutputLen + 1);
		dbMin = dbSpan;
		for (size_t i = 0, j = 0; i < nOutputLen; i++)
		{
			dbOutput[i] = dbMin;
			dbMin += dbSpan;
		}

	}
	else {

		double* pDensity = new double[resampleGridLen];		// density field


		// 1.calculate density function
		double dbAccum = 0.0;								// accumulated density
		for (size_t i = 0; i < resampleGridLen; i++)
		{
			double x = dbMin + dbStep * i;
			double y = 0;
			double dbH = .1;
			for (size_t j = 0; j < nInputLen; j++)
			{
				y += KernelFun((x - dbInput[j]) / dbH) / dbH;
			}
			pDensity[i] = y;
			//if(bShow) qDebug() << y;
			dbAccum += y;
		}

		// 2.resampling
		double dbCurrentAccum = 0.0;							// current accumulated density
		double dbDensityStep = dbAccum / (nOutputLen + 1);		// steps for the result density
		for (size_t i = 0, j = 0; i < resampleGridLen; i++)
		{
			double x = dbMin + dbStep * i;
			dbCurrentAccum += pDensity[i];


			//if (bShow) qDebug() << dbCurrentAccum << ",j:" << j << ",accum:" << j * dbDensityStep;
			if (dbCurrentAccum > (j + 1)*dbDensityStep) {
				dbOutput[j] = x;
				//if (bShow) qDebug() << x;
				j++;
			}
			if (j == nOutputLen)
			{
				break;
			}
		}
		delete[] pDensity;
	}
}

void resampleBasedOnKDE_FixLen(double* dbInput, double* dbOutput, int nInputLen, int nOutputLen) {
	int resampleGridLen = 1000;
	double dbRange = dbInput[nInputLen - 1] - dbInput[0];
	double dbMin = dbInput[0] - dbRange * 0.005;
	dbRange *= 1.01;
	double dbMax = dbMin + dbRange;
	double dbStep = dbRange / (resampleGridLen + 1);

	double* pDensity = new double[resampleGridLen];		// density field


	// 1.calculate density function
	double dbAccum = 0.0;								// accumulated density
	for (size_t i = 0; i < resampleGridLen; i++)
	{
		double x = dbMin + dbStep * i;
		double y = 0;
		double dbH = .1;
		for (size_t j = 0; j < nInputLen; j++)
		{
			y += KernelFun((x - dbInput[j]) / dbH) / dbH;
		}
		pDensity[i] = y;
		//if(bShow) qDebug() << y;
		dbAccum += y;
	}

	// 2.resampling
	double dbCurrentAccum = 0.0;							// current accumulated density
	double dbDensityStep = dbAccum / (nOutputLen + 1);		// steps for the result density
	for (size_t i = 0, j = 0; i < resampleGridLen; i++)
	{
		double x = dbMin + dbStep * i;
		dbCurrentAccum += pDensity[i];


		//if (bShow) qDebug() << dbCurrentAccum << ",j:" << j << ",accum:" << j * dbDensityStep;
		if (dbCurrentAccum > (j + 1)*dbDensityStep) {
			dbOutput[j] = x;
			//if (bShow) qDebug() << x;
			j++;
		}
		if (j == nOutputLen)
		{
			break;
		}
	}
	delete[] pDensity;
}

// resample buffer
void GridFrame::resampleBuf(const double* pSortedBuf, double* pResampledBuf, int nResampleLen) {
	double* dbSample = new double[_nEnsembleLen];
	double* dbResample = new double[nResampleLen];

	// 1.calculate resampled SDF for each grid point
	for (size_t i = 0; i < _nGrids; i++)
	{
		for (size_t l = 0; l < _nEnsembleLen; l++)
		{
			dbSample[l] = pSortedBuf[l*_nGrids + i];
		}
		if (g_bDebugDataField) {
			//qDebug() << dbSample[_nEnsembleLen - 1] - dbSample[0];
		}

		resampleBasedOnKDE_FixStep(dbSample, dbResample, _nEnsembleLen, nResampleLen);
		//resampleBasedOnKDE_FixLen(dbSample, dbResample, _nEnsembleLen, nResampleLen);
		for (size_t l = 0; l < nResampleLen; l++)
		{
			pResampledBuf[l*_nGrids + i] = dbResample[l];
		}
	}

	delete[] dbSample;
	delete[] dbResample;
}


// resample buffer
void GridFrame::resampleBuf_C(const double* pUnSortedBuf, const int* pLabels, double* pResampledBuf, int nClusters) {
	int nResampleLen = 11;
	double* dbSample = new double[_nEnsembleLen];
	double* dbResample = new double[nResampleLen];


	for (size_t indexC = 0; indexC < nClusters; indexC++){	// foreach Cluster
		// 0.get index list of this cluster
		std::vector<int> vecClusterMemberIndex;
		for (size_t l = 0; l < _nEnsembleLen; l++)
		{
			if (pLabels[l] == indexC) vecClusterMemberIndex.push_back(l);
		}

		// 1.calculate resampled SDF for each grid point
		for (size_t i = 0; i < _nGrids; i++)
		{
			int nClusterMemberLen = vecClusterMemberIndex.size();

			std::vector<double> vecValues;
			for (size_t l = 0; l < nClusterMemberLen; l++)
			{
				vecValues.push_back(pUnSortedBuf[vecClusterMemberIndex[l] * _nGrids + i]);
			}
			std::sort(vecValues.begin(), vecValues.end());

			for (size_t l = 0; l < nClusterMemberLen; l++)
			{
				dbSample[l] = vecValues[l];
			}


			resampleBasedOnKDE_FixStep(dbSample, dbResample, nClusterMemberLen, nResampleLen);
			//resampleBasedOnKDE_FixLen(dbSample, dbResample, _nEnsembleLen, nResampleLen);
			for (size_t l = 0; l < nResampleLen; l++)
			{
				pResampledBuf[l*_nGrids + i] = dbResample[l];
			}
		}
		pResampledBuf += _nGrids * nResampleLen;
	}


	delete[] dbSample;
	delete[] dbResample;
}