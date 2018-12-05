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
	/*
	// 1.clear the four group of spaghetti
	QList<QList<ContourLine>>& listContour = _pFeature->GetContours();
	QList<QList<ContourLine>>& listContourSDF = _pFeature->GetContours();
	QList<QList<ContourLine>>& listContourSorted = _pFeature->GetContours();
	QList<QList<ContourLine>>& listContourSortedSDF = _pFeature->GetContours();
	listContour.clear();
	listContourSDF.clear();
	listContourSorted.clear();
	listContourSortedSDF.clear();

	// 2.reset the scalar field to 1
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < _nWidth*_nHeight;j++)
		{
			_pData->SetData(i, j, 1);
		}
	}
	// 3.regenerate a contour line, and set -1 of the grid points of the other side, and calculate sdf
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
		_pFeature->calculateSDF(_pData->GetData(i), _pFeature->GetSDF(i), _nWidth, _nHeight, 0, contour);
	}

	// 4.regenerate contours from sorted SDF
	_pFeature->BuildSortedSDF();
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pFeature->GetSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			listContourSDF.push_back(contour);
		}
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pFeature->GetSortedSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			listContourSortedSDF.push_back(contour);
		}
	}
	*/
}