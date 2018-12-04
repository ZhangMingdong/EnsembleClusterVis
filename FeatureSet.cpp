#include "FeatureSet.h"
#include "DataField.h"
#include "ContourGenerator.h"

FeatureSet::FeatureSet(DataField* pData, double dbIsoValue, int nWidth, int nHeight, int nEnsembleLen, int nFocusX, int nFocusY, int nFocusW, int nFocusH):
	_pData(pData)
	,_dbIsoValue(dbIsoValue)
	,_nWidth(nWidth)
	,_nHeight(nHeight)
	,_nLen	 (nWidth*nHeight)
	,_nFocusX(nFocusX)
	,_nFocusY(nFocusY)
	,_nFocusW(nFocusW)
	,_nFocusH(nFocusH)
	,_nEnsembleLen(nEnsembleLen)
{
	
}


FeatureSet::~FeatureSet()
{
}

void FeatureSet::GenerateContours() {
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		QList<ContourLine> contour;
		ContourGenerator::GetInstance()->Generate(_pData->GetLayer(i), contour, g_fThreshold, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
		_listContour.push_back(contour);
	}
}