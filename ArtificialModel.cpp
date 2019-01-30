#include "ArtificialModel.h"
#include "DataField.h"
#include "FeatureSet.h"
#include "ContourGenerator.h"

#include <QDebug>


ArtificialModel::ArtificialModel()
{
	_type = MT_AF;
}


ArtificialModel::~ArtificialModel()
{
}


void ArtificialModel::generateDataFromContour() {
	// 1.regenerate data, changed contours and scalar field
	QList<QList<ContourLine>> listContour = regenerateData();

	// 2.calculate signed distance function
	double* arrSDF = new double[_nGrids*_nEnsembleLen];
	FeatureSet::CalculateSDF(0, _nEnsembleLen, _nWidth, _nHeight, _pData, arrSDF, listContour);

	// 3.replace the scalar field using the SDF
	for (int i = 0; i < _nEnsembleLen; i++)
	{
		for (int j = 0; j < _nGrids; j++)
		{
			_pData->SetData(i, j, arrSDF[i*_nGrids + j]);
		}
	}
	delete[] arrSDF;
}

const double c_dbK = 1.0 / sqrt(PI2d);
const double c_dbEpsilon = 1;
inline double KernelFun(double para) {
	return c_dbK * exp(-para * para / 2.0 / c_dbEpsilon / c_dbEpsilon);
}

void ArtificialModel::generateDataFromField() {
	double dbScale = 50.0;
	double arrVH[5];			// array of values of hight value
	double arrPH[5];			// array of horizontal positions of high value
	double arrVL[5];			// array of values of low value
	double arrPL[5];			// array of horizontal positions of low value
	double dbB = 5;				// band of patterns
	double dbB2 = 2;			// band of outliers
	double dbX, dbY, dbR, dbV;
	bool bPatternOutlier = false;

	// generate data for each ensemble emeber
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		// 1.generate the 5 pair of patterns
		for (int i = 0; i < 5; i++)
		{
			arrVH[i] = (double)rand() / RAND_MAX * 10+4;
			arrPH[i] = 2 + i * 20;
			arrVL[i] = (double)rand() / RAND_MAX * (-10)-4;
			arrPL[i] = 12 + i * 20;
			arrVH[i] *= dbScale;
			arrVL[i] *= dbScale;

		}

		for (int i = 0; i < _nHeight; i++)
		{
			for (int j = 0; j < _nWidth; j++)
			{
				double dbResult = 0;
				for (int k = 0; k < 5; k++)
				{
					// effect of high value
					dbX = j - arrPH[k];
					dbY = i - _nHeight + 1;
					if ((l < 10 && k == 1)
						|| (l >= 10 && l < 20 && k == 3))
					{
						dbY += 4;
					}
					dbR = sqrt(dbX*dbX + dbY * dbY);
					dbV = arrVH[k] * KernelFun(dbR / dbB) / dbB;
					dbResult += dbV;

					// effect of low value
					dbX = j - arrPL[k];
					dbY = i;
					if ((l >=40 && k == 1)
						|| (l >= 30 && l < 40 && k == 3))
					{
						dbY -= 15;
					}
					dbR = sqrt(dbX*dbX + dbY * dbY);
					dbV = arrVL[k] * KernelFun(dbR / dbB) / dbB;
					dbResult += dbV;
				}					

				// add a outlier
				if (bPatternOutlier) {
					if (l % 10 == 0)
					{
						dbX = j;
						dbY = i;
						dbR = sqrt(dbX*dbX + dbY * dbY);
						dbV = dbScale*10 * KernelFun(dbR / dbB2) / dbB2;
						dbResult += dbV;
					}
				}
				_pData->SetData(l, i*_nWidth+j, dbResult);
			}

		}
	}
}

void ArtificialModel::generateDataFromField_3() {
	double dbScale = 50.0;
	double arrVH[2];			// array of values of hight value
	double arrPH[2];			// array of horizontal positions of high value
	double arrVL[2];			// array of values of low value
	double arrPL[2];			// array of horizontal positions of low value
	double dbB = 5;				// band of patterns
	double dbB2 = 2;			// band of outliers
	double dbX, dbY, dbR, dbV;
	bool bPatternOutlier = false;

	// generate data for each ensemble emeber
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		// 1.generate the 5 pair of patterns
		for (int i = 0; i < 2; i++)
		{
			arrVH[i] = (double)rand() / RAND_MAX * 10 + 4;
			arrPH[i] = 12 + i * 40;
			arrVL[i] = (double)rand() / RAND_MAX * (-10) - 4;
			arrPL[i] = 22 + i * 40;
			arrVH[i] *= dbScale;
			arrVL[i] *= dbScale;

		}

		for (int i = 0; i < _nHeight; i++)
		{
			for (int j = 0; j < _nWidth; j++)
			{
				double dbResult = 0;
				for (int k = 0; k < 2; k++)
				{
					// effect of high value
					dbX = j - arrPH[k];
					dbY = i - _nHeight + 1;
					if (l < 10 && k == 0){
						dbY += 4;
					}
					if (l >= 10 && l < 20 && k == 1) {
						dbY += 4;
					}

					dbR = sqrt(dbX*dbX + dbY * dbY);
					dbV = arrVH[k] * KernelFun(dbR / dbB) / dbB;

				//	if (l >= 10 && l < 20 && k == 1) {
				//		dbV = 0;
				//	}
					dbResult += dbV;

					// effect of low value
					dbX = j - arrPL[k];
					dbY = i;
					if (l >= 40 && k == 0)
					{
						dbY -= 4;
					}
					if (l >= 30 && l < 40 && k == 1)
					{
						dbY -= 4;
					}
					dbR = sqrt(dbX*dbX + dbY * dbY);
					dbV = arrVL[k] * KernelFun(dbR / dbB) / dbB;
					dbResult += dbV;
				}

				// add a outlier
				if (bPatternOutlier) {
					if (l % 10 == 0)
					{
						dbX = j;
						dbY = i;
						dbR = sqrt(dbX*dbX + dbY * dbY);
						dbV = dbScale * 10 * KernelFun(dbR / dbB2) / dbB2;
						dbResult += dbV;
					}
				}
				_pData->SetData(l, i*_nWidth + j, dbResult);
			}

		}
	}
}

void ArtificialModel::generateDataFromField_4() {
	double dbScale = 50.0;
	double arrVH[2];			// array of values of hight value
	double arrPH[2];			// array of horizontal positions of high value
	double arrVL[2];			// array of values of low value
	double arrPL[2];			// array of horizontal positions of low value
	double dbB = 5;				// band of patterns
	double dbB2 = 2;			// band of outliers
	double dbX, dbY, dbR;
	bool bPatternOutlier = false;

	// generate data for each ensemble emeber
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		// 1.generate the 5 pair of patterns
		for (int i = 0; i < 2; i++)
		{
			arrVH[i] = (double)rand() / RAND_MAX * 10 + 4;
			arrPH[i] = 12 + i * 40;
			arrVL[i] = (double)rand() / RAND_MAX * (-10) - 4;
			arrPL[i] = 22 + i * 40;
			arrVH[i] *= dbScale;
			arrVL[i] *= dbScale;

		}

		for (int i = 0; i < _nHeight; i++)
		{
			for (int j = 0; j < _nWidth; j++)
			{
				double dbResult = 0;
				for (int k = 0; k < 2; k++)
				{
					// effect of high value
					dbX = j - arrPH[k];
					dbY = i - _nHeight + 1;

					dbR = sqrt(dbX*dbX + dbY * dbY);
					double dbVH = arrVH[k] * KernelFun(dbR / dbB) / dbB;

					//	if (l >= 10 && l < 20 && k == 1) {
					//		dbV = 0;
					//	}

					// effect of low value
					dbX = j - arrPL[k];
					dbY = i;

					dbR = sqrt(dbX*dbX + dbY * dbY);
					double dbVL = arrVL[k] * KernelFun(dbR / dbB) / dbB;

					if (l < 10 && k == 0) {
						dbVH *= -1;
					}
					if (l < 10 && k == 0) {
						dbVL *= -1;
					}
					if (l >= 10 && l < 20 && k == 1) {
						dbVL *= -1;
					}
					if (l >= 10 && l < 20 && k == 1) {
						dbVH *= -1;
					}

					dbResult += dbVH;
					dbResult += dbVL;
				}


				_pData->SetData(l, i*_nWidth + j, dbResult);
			}

		}
	}
}
void ArtificialModel::initializeModel() {
	// 2.Set isovalues
	QList<double> listIsoValue;
	listIsoValue.append(0);
	//listIsoValue.append(5);
	SetIsoValues(listIsoValue);


	//generateDataFromContour();
	//generateDataFromField();
	//generateDataFromField_2();
	//generateDataFromField_3();
	generateDataFromField_4();


	// 3.statistic
	_pData->DoStatistic();
	// EOF
	_pData->DoEOF();

	// 4.generate feature;
	for each (double isoValue in _listIsoValues)
	{
		_listFeature.append(new FeatureSet(_pData, isoValue, _nWidth, _nHeight, _nEnsembleLen));
	}


}

double generateY_1(int i, int j,int _nHeight,int _nEnsembleLen) {
	double y = _nHeight / 2
		+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX)
		+ (2.0*(i > (_nEnsembleLen / 2) ? 1 : -1));

	return y;
}

double generateY_2(int i, int j,int _nHeight,int _nEnsembleLen) {
	double y = _nHeight / 2
		+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX)
		+ (j / 10.0*(i > (_nEnsembleLen / 2) ? 1 : -1));

	return y;
}

double generateY_3(int i, int j, int _nHeight, int _nEnsembleLen) {
	double y = _nHeight / 2
		+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX);
	//+ (j / 10.0*(i > (_nEnsembleLen / 2) ? 1 : -1));
	if (i > 40)
		y += j / 10.0;
	else if (i < 10)
		y -= j / 10.0;

	return y;
}

double generateY_4(int i, int j,int _nWidth, int _nHeight, int _nEnsembleLen) {
	double y = _nHeight / 2
		+ sin((j + rand() / (double)RAND_MAX) / 2) * 2 * (rand() / (double)RAND_MAX);
	//+ (j / 10.0*(i > (_nEnsembleLen / 2) ? 1 : -1));

	if (i < 25) {
		if (j > _nWidth / 2)
			j = _nWidth / 2;
		if (i < 12)
			y += j / 8.0;
		else
			y -= j / 8.0;
	}
	else {
		j = _nWidth - 1 - j;
		if (j > _nWidth / 2)
			j = _nWidth / 2;
		if (i < 38)
			y += j / 12.0;
		else
			y -= j / 12.0;

	}

	return y;
}

QList<QList<ContourLine>> ArtificialModel::regenerateData(){
	QList<QList<ContourLine>> listContour;

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


			//double y = generateY_3(i, j, _nHeight, _nEnsembleLen);
			double y = generateY_4(i, j,_nWidth, _nHeight, _nEnsembleLen);




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
	return listContour;
}


void ArtificialModel::generateDataFromField_2() {
	// generate data for each ensemble emeber
	double dbB = 5;
	double dbMax = -100;
	double dbBg = -0.1;						// background value
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		double dbControlValue = (double)rand() / RAND_MAX * 100 + 45;

		double dbControlX = _nWidth / 2 + (double)rand() / RAND_MAX * 5 - 2.5;
		double dbControlY = _nHeight / 2 + (double)rand() / RAND_MAX * 5 - 2.5;

		for (int i = 0; i < _nHeight; i++)
		{
			for (int j = 0; j < _nWidth; j++)
			{
				double dbResult = dbBg;
				// effect of high value
				double dbX = j - dbControlX;
				double dbY = i - dbControlY;
				double dbR = sqrt(dbX*dbX + dbY * dbY);
				double dbV = dbControlValue * KernelFun(dbR / dbB) / dbB;
				dbResult += dbV;
				_pData->SetData(l, i*_nWidth + j, dbResult);
				//qDebug() << dbResult << endl;
				if (dbResult > dbMax) {
					qDebug() << dbResult;
					dbMax = dbResult;
				}
			}

		}
	}
}

















