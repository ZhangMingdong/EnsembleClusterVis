#include "ArtificialModel.h"
#include "DataField.h"
#include "FeatureSet.h"
#include "ContourGenerator.h"

#include <QDebug>


ArtificialModel::ArtificialModel()
{
}


ArtificialModel::~ArtificialModel()
{
}

void ArtificialModel::initializeModel() {
	// 1.regenerate data, changed contours and scalar field
	regenerateData();

	// 2.calculate signed distance
	_listFeature[0]->CalculateSDF(0);

	// 3.replace the scalar field using the SDF
	for (int i = 0; i < _nEnsembleLen; i++)
	{
		for (int j = 0; j < _nGrids; j++)
		{
			_pData->SetData(i, j, _listFeature[0]->GetSDF(i)[j]);
		}
	}

	// 4.update the data
	_pData->DoStatistic();

	// 5.update the feature
	FeatureSet* pFeature = new FeatureSet(_pData, 0, _nWidth, _nHeight, _nEnsembleLen, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	for each (FeatureSet* pFeature in _listFeature)
		delete pFeature;
	_listFeature.clear();
	_listFeature.append(pFeature);

}


void ArtificialModel::regenerateData() {
	// 1.clear contours
	QList<QList<ContourLine>>& listContour = _listFeature[0]->GetContours();
	listContour.clear();

	// 2.reset the scalar field to 1
	for (int i = 0; i < _nEnsembleLen; i++)
	{
		for (int j = 0; j < _nGrids; j++)
		{
			_pData->SetData(i, j, 1);
		}
	}

	// 3.regenerate a contour line, and set -1 of the grid points of the other side, and calculate sdf
	for (int i = 0; i < _nEnsembleLen; i++)
	{
		// generate contour
		ContourLine line;
		for (int j = 0; j < _nWidth; j++)
		{
			//double y = _nHeight / 2 + sin(j*i);

			//double y = _nHeight / 2 + sin(j+i)*2;

			//double y = _nHeight / 2 + sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX);
			//if (i == 0) y += 10;

			//qDebug() << "new data";
			//qDebug() << j / 100.0*(i > (_nEnsembleLen / 2) ?  1: -1);

			double y = _nHeight / 2
				+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX);
				//+ (j / 10.0*(i > (_nEnsembleLen / 2) ? 1 : -1));
			if (i > 40)
				y += j / 10.0;
			else if (i < 10)
				y -= j / 10.0;

			/*
			double y = _nHeight / 2
				+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX)
				+ (2.0*(i > (_nEnsembleLen / 2) ? 1 : -1));
				*/

			line._listPt.append(QPointF(j, y));
			for (int k = 0; k < y; k++)
			{
				_pData->SetData(i, k*_nWidth + j, -1);
			}
		}

		QList<ContourLine> contour;
		contour.append(line);
		listContour.push_back(contour);
	}
}