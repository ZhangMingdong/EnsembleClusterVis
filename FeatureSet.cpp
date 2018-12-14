#include "FeatureSet.h"
#include "DataField.h"
#include "ContourGenerator.h"
#include "UnCertaintyArea.h"
#include "ContourStripGenerator.h"

#include "MyPCA.h"
#include <KMeansClustering.h>
#include <AHCClustering.h>

#include <qDebug>


double PointToSegDist(double x, double y, double x1, double y1, double x2, double y2)
{
	double cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
	if (cross <= 0)
		return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));

	double d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
	if (cross >= d2)
		return sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));

	double r = cross / d2;
	double px = x1 + (x2 - x1) * r;
	double py = y1 + (y2 - y1) * r;
	return sqrt((x - px) * (x - px) + (y - py) * (y - py));
}

FeatureSet::FeatureSet(DataField* pData, double dbIsoValue, int nWidth, int nHeight, int nEnsembleLen, int nFocusX, int nFocusY, int nFocusW, int nFocusH):
	_pData(pData)
	,_dbIsoValue(dbIsoValue)
{
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nGrids = nWidth * nHeight;
	_nFocusX = nFocusX;
	_nFocusY = nFocusY;
	_nFocusW = nFocusW;
	_nFocusH = nFocusH;
	_nEnsembleLen = nEnsembleLen;

	_gridHalfMax = new double[_nGrids];
	_gridHalfMin = new double[_nGrids];
	_gridValidMax = new double[_nGrids];
	_gridValidMin = new double[_nGrids];
	_pSDF = new double[_nGrids*_nEnsembleLen];
	_pICD = new double[_nGrids];
	_pSortedSDF = new double[_nGrids*_nEnsembleLen];
	_pSet = new bool[_nGrids*_nEnsembleLen];
	_pGridDiverse = new bool[_nGrids];
	_pSetBandDepth = new int[_nEnsembleLen];
	_pMemberType = new int[_nEnsembleLen];

	// 1.calculate set
	calculateSet();

	// 2.calculate band depth
	calculateBandDepth();

	// 3.calculate member type
	calculateMemberType();

	// 4.calculate valid max and min and mean
	doStatistics();

	// 5.generate all the contours and the bands
	generateContours();

	// 6.calculate signed distance function, according to a given isovalue
	CalculateSDF(_dbIsoValue);

	// 7.sort according to sdf and generate contours according to sdf
	buildSortedSDF();

	// 8.calculate ICD
	calculateICD();

//==Clustering==
	calculateDiverse();

	calculatePCA();	

	//doClustering();

	//doPCAClustering();

	//calculatePCABox();

	//calculatePCARecovery();


}

FeatureSet::~FeatureSet()
{
	delete[] _gridHalfMax;
	delete[] _gridHalfMin;
	delete[] _gridValidMax;
	delete[] _gridValidMin;
	delete[] _pSDF;
	delete[] _pSortedSDF;
	delete[] _pICD;
	delete[] _pSet;
	delete[] _pGridDiverse;
	delete[] _pSetBandDepth;
	delete[] _pMemberType;

	for each (UnCertaintyArea* pArea in _listAreaValid)
		delete pArea;
	for each (UnCertaintyArea* pArea in _listAreaHalf)
		delete pArea;
	for each (UnCertaintyArea* pArea in _listUnionAreaE)
		delete pArea;
}

void FeatureSet::generateContours() {
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// spaghetti
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetLayer(i), contour, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContour.push_back(contour);
		}
		// sorted spaghetti
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetSortedLayer(i), contour, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSorted.push_back(contour);
		}
	}
	ContourGenerator::GetInstance()->Generate(_pData->GetUMin(), _listContourMinE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(_pData->GetUMax(), _listContourMaxE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(_pData->GetMean(), _listContourMeanE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetValidMin(), _listContourMinValid, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetValidMax(), _listContourMaxValid, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetHalfMin(), _listContourMinHalf, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetHalfMax(), _listContourMaxHalf, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
	ContourGenerator::GetInstance()->Generate(GetMedian(), _listContourMedianE, _dbIsoValue, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);


	generateContourImp(_listContourMinValid, _listContourMaxValid, _listAreaValid);
	generateContourImp(_listContourMinHalf, _listContourMaxHalf, _listAreaHalf);
}

void FeatureSet::buildSortedSDF() {

	sortBuf(_pSDF, _pSortedSDF);

	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// contours from SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(getSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSDF.push_back(contour);
		}
		// contours from sorted SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(GetSortedSDF(i), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSortedSDF.push_back(contour);
		}
	}
}

void FeatureSet::sortBuf(const double* pS, double* pD) {
	for (int i = 0; i < _nGrids; i++) {
		for (int j = 0; j < _nEnsembleLen; j++) {
			double dbValue = pS[j*_nGrids + i];
			int k = 0;
			while (k < j && pD[k*_nGrids + i] < dbValue) k++;
			int l = j;
			while (l > k) {
				pD[l*_nGrids + i] = pD[(l - 1)*_nGrids + i];
				l--;
			}
			pD[k*_nGrids + i] = dbValue;
		}
	}
}

void FeatureSet::calculateMemberType() {
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pMemberType[i] = 1;
	}
	int nOutliers = 0;
	// find outlier
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		if (_pSetBandDepth[i] < _nOutlierThreshold) {
			_pMemberType[i] = 0;
			nOutliers++;
		}
	}
	int nValidLen = (_nEnsembleLen - nOutliers) / 2;
	// find valid and 50%
	for (int i = 0; i < nValidLen; i++)
	{
		int nMaxDepth = 0;
		int nIndex = -1;
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (_pMemberType[j] == 1 && _pSetBandDepth[j] > nMaxDepth)
			{
				nIndex = j;
				nMaxDepth = _pSetBandDepth[j];
			}
		}
		_pMemberType[nIndex] = 2;
	}
}

void FeatureSet::doStatistics() {
	for (int i = 0; i < _nGrids; i++) {
		std::vector<double> vecDataValid;
		std::vector<double> vecDataHalf;
		// 1.calculate mean
		for (int j = 0; j < _nEnsembleLen; j++)
		{
			if (_pMemberType[j])
			{
				vecDataValid.push_back(_pData->GetData(j, i));
				if (_pMemberType[j] == 2)
				{
					vecDataHalf.push_back(_pData->GetData(j, i));
				}
			}
		}


		std::sort(vecDataValid.begin(), vecDataValid.end());
		std::sort(vecDataHalf.begin(), vecDataHalf.end());
		_gridHalfMin[i] = vecDataHalf[0];
		_gridHalfMax[i] = vecDataHalf[vecDataHalf.size() - 1];
		_gridValidMin[i] = vecDataValid[0];
		_gridValidMax[i] = vecDataValid[vecDataValid.size() - 1];

	}

	// find median
	int nMaxDepth = 0;
	_nMedianIndex = -1;
	for (int i = 0; i < _nEnsembleLen; i++) {
		if (_pSetBandDepth[i] > nMaxDepth)
		{
			nMaxDepth = _pSetBandDepth[i];
			_nMedianIndex = i;
		}
	}
}

void FeatureSet::calculateSet() {
	// 1.calculate set
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		for (int i = 0; i < _nGrids; i++) {
			_pSet[l*_nGrids + i] = (_pData->GetData(l,i) > _dbIsoValue);
		}
	}
}

// caclulate sBand depth
void FeatureSet::calculateBandDepth() {
	int mSet = (_nEnsembleLen - 1) * (_nEnsembleLen - 2) / 2;	// count of j=2 sets
	int* arrMatrix = new int[_nEnsembleLen*mSet];				// the matrix
	// 1.calculate the matrix
	for (size_t i = 0; i < _nEnsembleLen; i++){					// for each member
		int nSetIndex = 0;										// index of current set
		for (size_t j = 0; j < _nEnsembleLen; j++){
			if (i == j) continue;
			for (size_t k = 0; k < j; k++){
				if (i == k) continue;
				// for each set
				int nIntersection = 0;				// count intersection breaks
				int nUnion = 0;						// count uniton breaks
				for (int l = 0; l < _nGrids; l++)
				{
					if ((_pSet[i*_nGrids + l] && !_pSet[j*_nGrids + l] && !_pSet[k*_nGrids + l])) nIntersection++;
					if ((!_pSet[i*_nGrids + l] && _pSet[j*_nGrids + l] && _pSet[k*_nGrids + l])) nUnion++;
				}
				arrMatrix[i*mSet + nSetIndex] = (nIntersection > nUnion) ? nIntersection: nUnion;
				nSetIndex++;
			}
		}
	}
	// 2.calculate epsilon
	int nMatrixLength = _nEnsembleLen * mSet;
	double dbThreshold = nMatrixLength / 6.0;
	int epsilon = 1;
	while (true) {	// break when find a proper epsilon
		int nCount = 0;
		for (size_t i = 0; i < nMatrixLength; i++)
		{
			if (arrMatrix[i] < epsilon) nCount++;
		}
		if (nCount > dbThreshold) break;
		else epsilon++;
	}
//	qDebug() <<"epsilon:"<< epsilon;
//	qDebug() << "band depth:";
	// 3.claculate band depth
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pSetBandDepth[i] = 0;
		for (size_t j = 0; j < mSet; j++) {
			if (arrMatrix[i*mSet + j] < epsilon)_pSetBandDepth[i]++;
		}
//		qDebug() << _pSetBandDepth[i];
	}
//	qDebug() << "~band depth:";
	/*
	qDebug() << "test:";
	for (size_t j = 0; j < mSet; j++) {
		qDebug() << arrMatrix[j];
	}
	qDebug() << "test finished";
	*/

	delete[] arrMatrix;
}

void FeatureSet::calculateBandDepth_v1() {
	int nThreshold = 18;
	qDebug() << "test:";
	// caclulate sBand depth
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pSetBandDepth[i] = 0;
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			if (i == j) continue;
			for (size_t k = 0; k < j; k++)
			{
				if (i == k) continue;
				int nFailed = 0;
				for (int l = 0; l < _nGrids; l++)
				{
					if ((_pSet[i*_nGrids + l] && !_pSet[j*_nGrids + l] && !_pSet[k*_nGrids + l])
						|| (!_pSet[i*_nGrids + l] && _pSet[j*_nGrids + l] && _pSet[k*_nGrids + l])) {
						nFailed++;
					}
					//if (nFailed == nThreshold) break;
				}
				if (i == 0) qDebug() << nFailed;
				if (nFailed < nThreshold) _pSetBandDepth[i]++;
			}
		}
	}
	qDebug() << "test finished";
}

const double* FeatureSet::GetMedian() { return _pData->GetData(_nMedianIndex); }

void FeatureSet::generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas) {
	ContourStripGenerator generator;
	generator.Generate(areas, contourMin, contourMax, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
}

void FeatureSet::CalculateSDF(double dbIsoValue) {
	for (size_t l = 0; l < _nEnsembleLen; l++) {
		// calculate sdf
		const double* arrData = _pData->GetData(l);
		double* arrSDF = getSDF(l);
		QList<ContourLine> contour = _listContour[l];
		for (size_t i = 0; i < _nHeight; i++)
		{
			for (size_t j = 0; j < _nWidth; j++)
			{
				double dbMinDis = _nGrids;
				for (size_t k = 0, lenK = contour.length(); k < lenK; k++)
				{
					for (size_t l = 0, lenL = contour[k]._listPt.length() - 1; l < lenL; l++)
					{
						double dbDis = PointToSegDist(j, i
							, contour[k]._listPt[l].x(), contour[k]._listPt[l].y()
							, contour[k]._listPt[l + 1].x(), contour[k]._listPt[l + 1].y());
						if (dbDis < dbMinDis)
						{
							dbMinDis = dbDis;
						}

					}
				}
				if (arrData[i*_nWidth + j] < dbIsoValue) dbMinDis = -dbMinDis;	// sign
				arrSDF[i*_nWidth + j] = dbMinDis;
			}
		}
	}
}

void FeatureSet::calculatePCA() {

	// 1.set parameter
	int mI = _nDiverseCount;
	int mO = 2;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = createDiverseArray();
	double* arrOutput = new double[mO*n];

	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);
	// 4.generate points from the output
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << arrOutput[2*i];
		//qDebug() << arrOutput[i * 2] << "," << arrOutput[i * 2 + 1];
		//qDebug() << arrOutput[i * 3] << "," << arrOutput[i * 3 + 1] << "," << arrOutput[i * 3 + 2];
		//qDebug() << arrOutput[i * 4] << "," << arrOutput[i * 4 + 1] << "," << arrOutput[i * 4 + 2] << "," << arrOutput[i * 4 + 3];
	}
	qDebug() << "===========================";
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << arrOutput[2 * i+1];
	}
	// 5.release the buffer

	delete[] arrOutput;
	delete[] arrInput;
}

void FeatureSet::calculateDiverse() {
	_nDiverseCount = 0;
	for (size_t i = 0; i < _nGrids; i++)
	{
		bool bDiverse = false;

		for (size_t l = 1; l < _nEnsembleLen; l++)
		{
			if (_pSet[l*_nGrids + i] != _pSet[(l - 1)*_nGrids + i]) {
				bDiverse = true;
				_nDiverseCount++;
				break;
			}
		}
		_pGridDiverse[i] = bDiverse;

	}
	qDebug() << "Diverse count: " << _nDiverseCount;
}

double* FeatureSet::createDiverseArray() {
	double* arrBuf = new double[_nEnsembleLen*_nDiverseCount];
	for (size_t l = 0; l < _nEnsembleLen; l++)
	{
		for (size_t i = 0, j = 0; i < _nGrids; i++)
		{
			if (_pGridDiverse[i]) {
				arrBuf[l * _nDiverseCount + j] = _pSDF[l*_nGrids + i];
				j++;
			}
		}
	}
	return arrBuf;
}

void FeatureSet::doClustering() {
	// 2.cluster
	CLUSTER::Clustering* pClusterer = new CLUSTER::KMeansClustering();
	int nN = _nEnsembleLen;			// number of data items
	int nM = _nDiverseCount;		// dimension
	int nK = 2;						// clusters
	int* arrLabel = new int[nN];
	double* arrBuf = createDiverseArray();


	pClusterer->DoCluster(nN, nM, nK, arrBuf, arrLabel);
	qDebug() << "Labels";
	for (size_t i = 0; i < nN; i++)
	{
		qDebug() << arrLabel[i];
	}

	// 4.record label

	// 5.release the resouse
	delete arrLabel;
	delete pClusterer;
	delete arrBuf;
}

void FeatureSet::doPCAClustering() {

	// 1.set parameter
	int mI = _nDiverseCount;
	int mO = 5;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = createDiverseArray();
	double* arrOutput = new double[mO*n];

	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);

	// clustering
	{
		// 2.cluster
		CLUSTER::Clustering* pClusterer = new CLUSTER::KMeansClustering();
		int nN = _nEnsembleLen;			// number of data items
		int nM = mO;		// dimension
		int nK = 2;						// clusters
		int* arrLabel = new int[nN];
		double* arrBuf = arrOutput;


		pClusterer->DoCluster(nN, nM, nK, arrBuf, arrLabel);
		qDebug() << "Labels";
		for (size_t i = 0; i < nN; i++)
		{
			qDebug() << arrLabel[i];
		}

		// 4.record label

		// 5.release the resouse
		delete arrLabel;
		delete pClusterer;
	}



	delete[] arrOutput;
	delete[] arrInput;
}

void FeatureSet::calculatePCARecovery() {

	// 1.set parameter
	int mI = _nDiverseCount;
	int mO = 1;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = createDiverseArray();
	double* arrOutput = new double[mO*n];

	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);

	/*
	double dbMax0 = -100;
	double dbMax1 = -100;
	double dbMin0 = 100;
	double dbMin1 = 100;
	// 4.generate points from the output
	for (size_t i = 0; i < n; i++)
	{
		if (arrOutput[i * 2] > dbMax0)
			dbMax0 = arrOutput[i * 2];
		if (arrOutput[i * 2] < dbMin0)
			dbMin0 = arrOutput[i * 2];
		if (arrOutput[i * 2 + 1] > dbMax1)
			dbMax1 = arrOutput[i * 2];
		if (arrOutput[i * 2 + 1] < dbMin1)
			dbMin1 = arrOutput[i * 2];
	}
	qDebug() << dbMin0 << "," << dbMin1 << "," << dbMax0 << "," << dbMax1;
	// recover
	double box[8] = {dbMin0, dbMin1
		, dbMin0, dbMax1
		, dbMax0, dbMin1
		, dbMax0, dbMax1 };
	double box[8] = { arrOutput[0],arrOutput[1]
		, arrOutput[2],arrOutput[3]
		, arrOutput[60],arrOutput[61]
		, arrOutput[62],arrOutput[63] };

	double *arrRecoveredBuf = new double[4 * _nDiverseCount];
	//pca.Recover(box, arrRecoveredBuf, 4,mO);
	for (size_t i = 0; i < 4; i++)
	{
		pca.Recover(box+i*2, arrRecoveredBuf+i*_nDiverseCount, 2);
	}
*/

/*
for (size_t i = 0; i < 10; i++)
{
	qDebug() << arrRecoveredBuf[i];
}
*/

	int nTestLen = 10;
	double *arrRecoveredBuf = new double[nTestLen* _nDiverseCount];
	//pca.Recover(box, arrRecoveredBuf, 4,mO);

	pca.Recover(arrOutput, arrRecoveredBuf, nTestLen,mO);

	for (size_t l = 0; l < nTestLen; l++)
	{
		//pca.Recover(arrOutput + l * mO, arrRecoveredBuf + l * _nDiverseCount, 1, mO);
		for (size_t i = 0, j = 0; i < _nGrids; i++)
		{
			if (_pGridDiverse[i]) {
				_pSDF[l*_nGrids + i] = arrRecoveredBuf[l * _nDiverseCount + j];
				j++;
			}
		}
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(getSDF(l), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSDF[l] = contour;
		}
	}

	// 5.release the buffer
	delete[] arrRecoveredBuf;

	delete[] arrOutput;
	delete[] arrInput;

}

void FeatureSet::calculatePCABox() {

	// 1.set parameter
	int mI = _nDiverseCount;
	int mO = 2;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = createDiverseArray();
	double* arrOutput = new double[mO*n];

	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, arrOutput, n, mI, mO, true);
	if(false)
		for (size_t i = 0; i < n; i++)
		{
			//qDebug() << arrOutput[i];
			qDebug() << arrOutput[i * 2] << "," << arrOutput[i * 2 + 1];
			//qDebug() << arrOutput[i * 3] << "," << arrOutput[i * 3 + 1] << "," << arrOutput[i * 3 + 2];
			//qDebug() << arrOutput[i * 4] << "," << arrOutput[i * 4 + 1] << "," << arrOutput[i * 4 + 2] << "," << arrOutput[i * 4 + 3];
		}

	double dbMax0 = -10000;
	double dbMax1 = -10000;
	double dbMin0 = 10000;
	double dbMin1 = 10000;
	// 4.generate points from the output
	for (size_t i = 0; i < n; i++)
	{
		if (arrOutput[i * 2] > dbMax0)
			dbMax0 = arrOutput[i * 2];
		if (arrOutput[i * 2] < dbMin0)
			dbMin0 = arrOutput[i * 2];
		if (arrOutput[i * 2 + 1] > dbMax1)
			dbMax1 = arrOutput[i * 2 + 1];
		if (arrOutput[i * 2 + 1] < dbMin1)
			dbMin1 = arrOutput[i * 2 + 1];
	}

	qDebug() << dbMin0 << "," << dbMin1 << "," << dbMax0 << "," << dbMax1;
	// recover
	double box[8] = {dbMin0, dbMin1
		, dbMin0, dbMax1
		, dbMax0, dbMin1
		, dbMax0, dbMax1 };
	/*
	double box[8] = { arrOutput[0],arrOutput[1]
		, arrOutput[2],arrOutput[3]
		, arrOutput[60],arrOutput[61]
		, arrOutput[62],arrOutput[63] };
		*/

	double *arrRecoveredBuf = new double[4 * _nDiverseCount];
	pca.Recover(box, arrRecoveredBuf, 4,mO);


	for (size_t l = 0; l < 4; l++)
	{
		for (size_t i = 0, j = 0; i < _nGrids; i++)
		{
			if (_pGridDiverse[i]) {
				_pSDF[l*_nGrids + i] = arrRecoveredBuf[l * _nDiverseCount + j];
				j++;
			}			
		}
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(getSDF(l), contour, 0, _nWidth, _nHeight, _nFocusX, _nFocusY, _nFocusW, _nFocusH);
			_listContourSDF[l]=contour;
		}
	}

	// 5.release the buffer
	delete[] arrRecoveredBuf;

	delete[] arrOutput;
	delete[] arrInput;

}

void FeatureSet::calculateICD() {
	double dbMin = 100;
	double dbMax = -100;
	for (size_t i = 0; i < _nGrids; i++)
	{
		double dbSum = 0;
		for (size_t l = 0; l < _nEnsembleLen; l++)
		{
			double dbValue = _pData->GetData(l, i);
			dbSum += dbValue;
		}
		double dbMean = dbSum / _nEnsembleLen;
		double dbVar = 0;
		for (size_t l = 0; l < _nEnsembleLen; l++)
		{
			double dbValue = _pData->GetData(l, i);
			double dbBias = dbValue - dbMean;
			dbVar += dbBias * dbBias;
		}

		double dbBias = _dbIsoValue - dbMean;
		_pICD[i] = pow(Ed,-dbBias*dbBias/2*dbVar);

		if (_pICD[i] > dbMax) dbMax = _pICD[i];
		if (_pICD[i] < dbMin) dbMin = _pICD[i];

	}
	qDebug() << "Min:" << dbMin;
	qDebug() << "Max:" << dbMax;


	return;
	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICD[i] = .5;
	}
}