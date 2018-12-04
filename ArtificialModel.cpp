#include "ArtificialModel.h"
#include "DataField.h"
#include "FeatureSet.h"
#include "ContourGenerator.h"


ArtificialModel::ArtificialModel()
{
}


ArtificialModel::~ArtificialModel()
{
}

void ArtificialModel::initializeModel() {
	QList<QList<ContourLine>>& listContour = _pFeature->GetContours();
	listContour.clear();
	_listContourSDF.clear();
	_listContourSorted.clear();
	_listContourSortedSDF.clear();

	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < _nWidth*_nHeight;j++)
		{
			_pData->SetData(i, j, 1);
		}
	}
	
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// generate contour
		QList<ContourLine> contour;
		ContourLine line;
		for (size_t j = 0; j < _nWidth; j++)
		{
			//double y = _nHeight / 2 + sin(j*i);

			//double y = _nHeight / 2 + sin(j+i)*2;

			double y = _nHeight / 2 + sin((j + rand()/(double)RAND_MAX)/2) * 2* (rand() / (double)RAND_MAX);

			line._listPt.append(QPointF(j, y));
			for (size_t k = 0; k < y; k++)
			{
				_pData->SetData(i, k*_nWidth + j, -1);
			}
		}
		contour.append(line);
		listContour.push_back(contour);

		// calculate sdf
		calculateSDF(_pData->GetData(i), _pData->GetSDF(i), _nWidth, _nHeight, 0, contour);
	}

	// contours from sorted SDF
	_pData->BuildSortedSDF();
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSDF.push_back(contour);
		}
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetSortedSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSortedSDF.push_back(contour);
		}
	}
}