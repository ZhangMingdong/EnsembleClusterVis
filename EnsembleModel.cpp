#include "EnsembleModel.h"
#include "def.h"
#include "DataField.h"
#include "MyPCA.h"

#include <qfile.h>
#include <qdebug.h>
#include <qmessagebox.h>

#include <fstream>


using namespace std;

EnsembleModel::EnsembleModel()
{
}


EnsembleModel::~EnsembleModel()
{
	for each (vector<DataField*> vec in _matrixData)
	{
		for each (DataField* pData in vec)
		{
			delete pData;
		}
	}
}

void EnsembleModel::readDataFromText() {
	int nTimeStep = g_nTimeStep;
	QFile file(_strFile);
	if (!file.open(QIODevice::ReadOnly)) {
		QMessageBox::information(0, "error", file.errorString());
	}
	int nCount = 0;
	QTextStream in(&file);
	int tt = 0;
	while (!in.atEnd()) {
		for (size_t iDay = 0; iDay < _nDay; iDay++)
		{
			vector<DataField*> vecDay;
			for (size_t iStep = 0; iStep < _nStep; iStep++)
			{
				DataField* pData = new DataField(_nWidth, _nHeight, _nEnsembleLen);
				for (size_t iEns = 0; iEns < _nEnsembleLen; iEns++)
				{
					QString line = in.readLine();
//					qDebug() << line;
					// every grid
					for (int i = 0; i < _nLen; i++)
					{
						QString line = in.readLine();
						pData->SetData(iEns, i, line.toFloat());
					}
				}
				vecDay.push_back(pData);
			}
			_matrixData.push_back(vecDay);
		}
	}
	qDebug() << "Finished";
	file.close();
	return;
	ofstream output("2D.txt");

	// PCA
	int mI = 961;
	int mO = 2;
	int n = _nEnsembleLen*_nStep *_nDay;
	double* arrInput = new double[mI*n];
	double* arrOutput = new double[mO*n];
	for (size_t iDay = 0; iDay < _nDay; iDay++)
	{
		for (size_t iStep = 0; iStep < _nStep; iStep++)
		{
			for (size_t iEns = 0; iEns < _nEnsembleLen; iEns++)
			{
				for (size_t i = 0; i < mI; i++)
				{
					arrInput[((iDay*_nStep + iStep)*_nEnsembleLen+iEns)*mI + i] = _matrixData[iDay][iStep]->GetData(iEns, i);
				}
			}
		}
	}
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);

	for (size_t i = 0; i < n; i++)
	{
		_vecPoints.push_back(DPoint3(arrOutput[i*mO], arrOutput[i*mO + 1], 0));
		output << arrOutput[i*mO]<<"\t"<<arrOutput[i*mO + 1] << endl;
	}


	delete[] arrInput;
	delete[] arrOutput;

}



GLubyte* EnsembleModel::generateTextureRange(int nIndex) {
	GLubyte* dataTexture = new GLubyte[4 * _nFocusLen];

//	const double *pMean = _matrixData[0][0]->GetLayer(nIndex);
	nIndex = 28;
	_matrixData[nIndex][0]->DoStatistic();
	const double *pMean = _matrixData[nIndex][0]->GetMean();

	double fMin = 100000;
	double fMax = -100000;

	for (int i = 0; i < _nLen; i++)
	{
		double f = _pData->GetData(0, i);
		if (pMean[i] > fMax) fMax = pMean[i];
		if (pMean[i] < fMin) fMin = pMean[i];
	}
	double fRange = fMax - fMin;

	for (int i = _nFocusY, iLen = _nFocusY + _nFocusH; i < iLen; i++) {
		for (int j = _nFocusX, jLen = _nFocusX + _nFocusW; j < jLen; j++) {

			int nIndex = i*_nWidth + j;
			int nIndexFocus = (i - _nFocusY)*_nFocusW + j - _nFocusX;
			// using transparency and the blue tunnel
			dataTexture[4 * nIndexFocus + 0] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 1] = (GLubyte)0;
			dataTexture[4 * nIndexFocus + 2] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
			dataTexture[4 * nIndexFocus + 3] = (GLubyte)((pMean[nIndex] - fMin) * 254 / fRange);
		}
	}

	return dataTexture;
}