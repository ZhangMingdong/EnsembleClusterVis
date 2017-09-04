#include "ContourBandDepthModel.h"
#include "DataField.h"
#include "def.h"

#include <iostream>

ContourBandDepthModel::ContourBandDepthModel()
{
}

ContourBandDepthModel::~ContourBandDepthModel()
{
}

void ContourBandDepthModel::initializeModel() {
	calculateCBD(g_fThreshold);
}

void ContourBandDepthModel::calculateCBD(double dbIso) {
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		int nDepth = 0;
		// calculate the i'th member
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (j == i) continue;
			for (size_t k = 0; k < _nEnsembleLen; k++)
			{
				if (k == i || k == j) continue;
				// test band (j,k);
				const double* arrFieldI = _pData->GetLayer(i);
				const double* arrFieldJ = _pData->GetLayer(j);
				const double* arrFieldK = _pData->GetLayer(k);	
				if (testDepth(dbIso,arrFieldI, arrFieldJ, arrFieldK)) nDepth++;
			}
		}
		_vecContourBandDepth.push_back(nDepth);
	}
}

bool ContourBandDepthModel::testDepth(double dbIso
	, const double* fieldI
	, const double* fieldJ
	, const double* fieldK) {
	int nTotal = 0;					// total iso points
	int nTarget = 0;				// target iso points
	for (size_t i = 0; i < _nHeight; i++)
	{
		for (size_t j = 0; j < _nWidth; j++)
		{
			double dbI0 = fieldI[i*_nWidth + j];
			double dbJ0 = fieldJ[i*_nWidth + j];
			double dbK0 = fieldK[i*_nWidth + j];
			if (i<_nHeight-1)
			{
				double dbI1 = fieldI[(i + 1)*_nWidth + j];
				double dbJ1 = fieldJ[(i + 1)*_nWidth + j];
				double dbK1 = fieldK[(i + 1)*_nWidth + j];
				if (qMin(dbI0,dbI1)<dbIso&&qMax(dbI0,dbI1)>dbIso)
				{
					nTotal++;
					if (qMin(dbJ0, dbK0)<dbI0
						&& qMax(dbJ0, dbK0)>dbI0
						&& qMin(dbJ1, dbK1)<dbI1
						&& qMax(dbJ1, dbK1)>dbI1
						)
					{
						nTarget++;
					}
				}
			}
			if (j < _nWidth - 1) {
				double dbI1 = fieldI[i*_nWidth + j + 1];
				double dbJ1 = fieldJ[i*_nWidth + j + 1];
				double dbK1 = fieldK[i*_nWidth + j + 1];
				if (qMin(dbI0, dbI1)<dbIso&&qMax(dbI0, dbI1)>dbIso)
				{
					nTotal++;
					if (qMin(dbJ0, dbK0)<dbI0
						&& qMax(dbJ0, dbK0)>dbI0
						&& qMin(dbJ1, dbK1)<dbI1
						&& qMax(dbJ1, dbK1)>dbI1
						)
					{
						nTarget++;
					}
				}
			}
		}
	}
	return nTarget/(double)nTotal>.5;
}

