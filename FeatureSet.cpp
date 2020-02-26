#include "FeatureSet.h"
#include "DataField.h"
#include "ContourGenerator.h"
#include "UnCertaintyArea.h"
#include "ContourStripGenerator.h"

#include "MyPCA.h"
#include <KMeansClustering.h>
#include <AHCClustering.h>

#include <qDebug>

#include <Eigen/Core>
#include <mathtoolbox/classical-mds.hpp>
#include <iostream>

#include "Switch.h"

using namespace Eigen;



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

int FeatureSet::GetLabel(int l) { 
	if(g_bClustering)
		return _arrLabels[l];
	else return 0;
}

FeatureSet::FeatureSet(DataField* pData, double dbIsoValue, int nWidth, int nHeight, int nEnsembleLen):
	_pData(pData)
	,_dbIsoValue(dbIsoValue)
{
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nGrids = nWidth * nHeight;
	_nEnsembleLen = nEnsembleLen;


	// 0.allocate resource
	_gridHalfMax = new double[_nGrids];
	_gridHalfMin = new double[_nGrids];
	_gridValidMax = new double[_nGrids];
	_gridValidMin = new double[_nGrids];

	_pSDF = new double[_nGrids*_nEnsembleLen];
	_pSortedSDF = new double[_nGrids*_nEnsembleLen];
	_pResampledSDF = new double[_nGrids*_nResampleLen];
	_pResampledSDF_C = new double[_nGrids*_nClusters*_nResampleLen_C];

	_pICDVX = new double[_nGrids*_nEnsembleLen];
	_pICDVY = new double[_nGrids*_nEnsembleLen];
	_pICDVW = new double[_nGrids*_nEnsembleLen];
	_pICD_LineKernel = new double[_nGrids];
	_pICDX = new double[_nGrids];
	_pICDY = new double[_nGrids];
	_pICDZ = new double[_nGrids];
	_pICD = new double[((_nWidth - 1)*_nDetailScale + 1)*((_nHeight - 1)*_nDetailScale + 1)];
	_pUpperSet = new bool[_nGrids*_nEnsembleLen];
	_pGridDiverse = new bool[_nGrids];
	_pSetBandDepth = new int[_nEnsembleLen];
	_pMemberType = new int[_nEnsembleLen];

	_arrLabels = new int[_nEnsembleLen];
	for (size_t i = 0; i < _nEnsembleLen; i++) _arrLabels[i] = 0;
	_arrMergeTarget = new int[_nEnsembleLen];
	_arrMergeSource = new int[_nEnsembleLen];
	_arrPC = new double[_nEnsembleLen*_nPCLen];





	if (g_bStatistic) {	
		// 1.calculate set
		calculateLevelSet();

		// 2.calculate band depth
		calculateBandDepth();

		// 3.calculate member type
		calculateMemberType();

		// 4.calculate valid max and min and mean
		doStatistics();

	}
	// 5.generate all the contours and the bands
	if(g_bGenerateContours) generateContours();

	// *5.generate smoothed contours
	if (g_bSmoothContours) smoothContours();

	// 6.generate bands
	if (g_bGenerateArea) {
		generateContourImp(_listContourMinValid, _listContourMaxValid, _listAreaValid);
		generateContourImp(_listContourMinHalf, _listContourMaxHalf, _listAreaHalf);
	}


	// 7.SDFs, sorted SDFs, contours from SDFs and sorted SDFs
	if(g_bCalculateSDF) buildSDF();

	// 8.Resample contours
	if(g_bResampleContours) resampleContours();

	// *8.Generate contours using domain sampled data
	if (g_bDomainResampleContours)
	{
		// domain space
		for (size_t l = 0; l < _nResampleLen; l++)
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetResampledData(l), contour, _dbIsoValue, _nWidth, _nHeight);
			_listContourDomainResampled.push_back(contour);
		}
	}

	if (g_bTestInfoLoseMeasure) {
		for (_nTestSamples=1;_nTestSamples<=_nEnsembleLen;_nTestSamples++)
		//for (size_t i = 0; i < 10; i++)
		{
			_nResampleLen = _nTestSamples;
			if(g_nTestTarget ==0) resampleBuf(_pSortedSDF, _pResampledSDF, _nResampleLen);
			g_bMutualInfoLose ? measureInfoLoseMI(): measureInfoLose();
		}
	}

	// 9. calculate iso-contour density
	if (g_bCalculateICD) {
		calculateICD();
		buildICD_LineKernel();
		buildICD_Vector();
		buildICDV();
	}

	// 10.Clustering
	if (g_bClustering)
	{
		// 10.1.calculate diversity
		calculateDiverse();

		// 10.2.calculate PCA

		//calculatePCA();
		//calculatePCA_Whole();
		//calculatePCA_MDS();
		//calculatePCA_MDS_Set();
		calculatePCA_MDS_Whole();
		//calculatePCA_MDS_Whole_Density();
		//calculatePCA_MutualInformation();

		// 10.3.clustering
		doPCAClustering();
		//doClustering();

		// 10.4.resample clusters
		if (g_bResampleForClusters)
		{
			resampleContours_C();
		}
	}

	//calculatePCABox();

	//calculatePCARecovery();


}

/*
	write these function according to paper "Isosurface Similarity Maps"
	2019/01/27
*/
double FeatureSet::calculateSimilarity(int l1, int l2,double dbMin, double dbMax) {
	const int nBin = 128;
	int arrBin[nBin][nBin];
	int arrBinX[nBin];
	int arrBinY[nBin];
	double dbRange = dbMax - dbMin + 1;

	for (size_t i = 0; i < nBin; i++) {
		arrBinX[i] = 0;
		arrBinY[i] = 0;
		for (size_t j = 0; j < nBin; j++) arrBin[i][j] = 0;
	}

	// calculate bin
	for (size_t i = 0; i < _nGrids; i++)
	{
		int nIndex1 = (_pSDF[l1*_nGrids + i] - dbMin) / dbRange * 128;
		int nIndex2 = (_pSDF[l2*_nGrids + i] - dbMin) / dbRange * 128;
		arrBin[nIndex1][nIndex2]++;
		arrBinX[nIndex1]++;
		arrBinY[nIndex2]++;
	}

	// calculate mutual information
	double dbHx = 0;
	double dbHy = 0;
	double dbHxy = 0;
	for (size_t i = 0; i < nBin; i++)
	{
		if (arrBinX[i]) {
			double dbPx = arrBinX[i] / (double)_nGrids;
			dbHx += dbPx * log(dbPx);
		}
		if(arrBinY[i]){
			double dbPy = arrBinY[i] / (double)_nGrids;
			dbHy += dbPy * log(dbPy);
		}
		for (size_t j = 0; j < nBin; j++)
		{
			if(arrBin[i][j]){
				double dbPxy = arrBin[i][j] / (double)_nGrids;
				dbHxy += dbPxy * log(dbPxy);
			}
		}
	}
	return -(dbHx + dbHy - dbHxy);
}

void FeatureSet::calculateSimilarityMatrix() {
	// calculate the range of sdf;
	double dbMax = -10000;
	double dbMin = 10000;
	for (size_t i = 0, length = _nEnsembleLen * _nGrids; i < length; i++)
	{
		if (_pSDF[i] > dbMax) dbMax = _pSDF[i];
		if (_pSDF[i] < dbMin) dbMin = _pSDF[i];
	}
	qDebug() << "calculateSimilarityMatrix";
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		double dbSimilarityDistribution = 0;
		for (size_t j = 0; j < _nEnsembleLen; j++)
		{
			double dbSimilarity = calculateSimilarity(i, j,dbMin,dbMax);
			dbSimilarityDistribution += dbSimilarity;
			//qDebug() << dis;
		}
		qDebug() << dbSimilarityDistribution;
	}
	qDebug() << "~calculateSimilarityMatrix";
}

FeatureSet::~FeatureSet()
{
	delete[] _gridHalfMax;
	delete[] _gridHalfMin;
	delete[] _gridValidMax;
	delete[] _gridValidMin;
	delete[] _pSDF;
	delete[] _pICDVX;
	delete[] _pICDVY;
	delete[] _pICDVW;
	delete[] _pICD_LineKernel;
	delete[] _pICDX;
	delete[] _pICDY;
	delete[] _pICDZ;
	delete[] _pSortedSDF;
	delete[] _pICD;
	delete[] _pUpperSet;
	delete[] _pGridDiverse;
	delete[] _pSetBandDepth;
	delete[] _pMemberType;
	delete[] _pResampledSDF;
	delete[] _pResampledSDF_C;
	delete[] _arrLabels;
	delete[] _arrMergeTarget;
	delete[] _arrMergeSource;

	delete[] _arrPC;

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
			ContourGenerator::GetInstance()->Generate(_pData->GetLayer(i), contour, _dbIsoValue, _nWidth, _nHeight);
			_listContour.push_back(contour);
		}
		// sorted spaghetti
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(_pData->GetSortedLayer(i), contour, _dbIsoValue, _nWidth, _nHeight);
			_listContourSorted.push_back(contour);
		}
	}
	ContourGenerator::GetInstance()->Generate(_pData->GetUMin(), _listContourMinE, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(_pData->GetUMax(), _listContourMaxE, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(_pData->GetMean(), _listContourMeanE, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(GetValidMin(), _listContourMinValid, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(GetValidMax(), _listContourMaxValid, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(GetHalfMin(), _listContourMinHalf, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(GetHalfMax(), _listContourMaxHalf, _dbIsoValue, _nWidth, _nHeight);
	ContourGenerator::GetInstance()->Generate(GetMedian(), _listContourMedianE, _dbIsoValue, _nWidth, _nHeight);

}
/*
	smooth contours
	using contours
*/
void FeatureSet::smoothContours() {
	double dbScale = .2;
	for each (QList<ContourLine> contours in _listContour)
	{
		QList<ContourLine> contours_s;
		for each (ContourLine contour in contours)
		{
			int nPoints = contour._listPt.length();
			if (nPoints<4)
			{
				contours_s.append(contour);
			}
			else {
				ContourLine contour_s;
				contour_s._listPt.append(contour._listPt[0]);
				contour_s._listPt.append(contour._listPt[1]);
				for (size_t i = 1; i < nPoints-2; i++)
				{
					QPointF pt1 = contour._listPt[i-1];
					QPointF pt2 = contour._listPt[i];
					QPointF pt3 = contour._listPt[i+1];
					QPointF pt4 = contour._listPt[i+2];

					QPointF pt14 = (pt1 + pt4) / 2;
					QPointF pt23 = (pt2 + pt3) / 2;
					QPointF pt = pt23 + (pt23 - pt14)*dbScale;
					contour_s._listPt.append(pt);
					contour_s._listPt.append(pt3);
				}



				contour_s._listPt.append(contour._listPt[nPoints - 1]);
				contours_s.append(contour_s);

			}
		}

		_listContourSmooth.append(contours_s);
	}
}

/*
	calculate signed distance function, according to a given isovalue
	sort according to sdf and generate contours according to sdf
*/
void FeatureSet::buildSDF() {
	// 1.build SDF
	CalculateSDF(_dbIsoValue, _nEnsembleLen, _nWidth, _nHeight, _pData, _pSDF, _listContour);

	// 2.Sort SDF
	sortBuf(_pSDF, _pSortedSDF);

	// 3.generate contours from SDF and sorted SDF
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		// contours from SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(getSDF(i), contour, 0, _nWidth, _nHeight);
			_listContourSDF.push_back(contour);
		}
		// contours from sorted SDF
		{
			QList<ContourLine> contour;
			ContourGenerator::GetInstance()->Generate(GetSortedSDF(i), contour, 0, _nWidth, _nHeight);
			_listContourSortedSDF.push_back(contour);
		}
	}
}

/*
	resample the sdf and generate resampled contours
	using Sorted SDF
*/
void FeatureSet::resampleContours() {
	// 1.Resample buffer
	resampleBuf(_pSortedSDF, _pResampledSDF, _nResampleLen);	

	// 2.generate contours for each grid point
	for (size_t l = 0; l < _nResampleLen; l++)
	{
		QList<ContourLine> contour;
		ContourGenerator::GetInstance()->Generate(_pResampledSDF + l * _nGrids, contour, 0, _nWidth, _nHeight);
		_listContourResampled.push_back(contour);
	}
}

void FeatureSet::resampleContours_C() {
	// 1.Resample buffer
	resampleBuf_C(_pSDF, _arrLabels, _pResampledSDF_C, _nClusters);

	// 2.generate contours for each grid point
	for (size_t l = 0; l < _nClusters*11; l++)
	{
		QList<ContourLine> contour;
		ContourGenerator::GetInstance()->Generate(_pResampledSDF_C + l * _nGrids, contour, 0, _nWidth, _nHeight);
		_listContourResampled_C.push_back(contour);
	}
	qDebug() << "resampleContours_C:" << _listContourResampled_C.length();
}

/*
	calculate type of each member: 0-outlier; 1-100%; 2-50%
	using SetBandDepth
*/
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
	int nValidHalfLen = (_nEnsembleLen - nOutliers) / 2;
	// find valid and 50%
	for (int i = 0; i < nValidHalfLen; i++)
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

/*
	calculate meadian, valid, and 50%
	using MemberType and BandDepth
*/
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

void FeatureSet::calculateLevelSet() {
	// 1.calculate set
	for (int l = 0; l < _nEnsembleLen; l++)
	{
		for (int i = 0; i < _nGrids; i++) {
			_pUpperSet[l*_nGrids + i] = (_pData->GetData(l,i) > _dbIsoValue);
		}
	}
}

/*
	caclulate sBand depth
	using UpperSet;
*/
void FeatureSet::calculateBandDepth() {
	int mSet = (_nEnsembleLen - 1) * (_nEnsembleLen - 2) / 2;	// count of j=2 sets
	int nMatrixLength = _nEnsembleLen * mSet;
	int* arrMatrix = new int[nMatrixLength];					// the matrix
	// 1.calculate the matrix, n*C^2_n
	for (size_t i = 0; i < _nEnsembleLen; i++){					// for each member
		int nSetIndex = 0;										// index of current set
		for (size_t j = 0; j < _nEnsembleLen; j++){
			if (i == j) continue;
			for (size_t k = 0; k < j; k++){
				if (i == k || j==k) continue;
				// for each set
				int nBreaks = 0;				// count of breaks(i not between j and k)
				for (int l = 0; l < _nGrids; l++)
				{
					if ((_pUpperSet[i*_nGrids + l] && !_pUpperSet[j*_nGrids + l] && !_pUpperSet[k*_nGrids + l])) nBreaks++;
					if ((!_pUpperSet[i*_nGrids + l] && _pUpperSet[j*_nGrids + l] && _pUpperSet[k*_nGrids + l])) nBreaks++;
				}
				arrMatrix[i*mSet + nSetIndex] = nBreaks;
				nSetIndex++;
			}
		}
	}
	// 2.calculate epsilon
	double dbThreshold = nMatrixLength / 6.0;
	/*
		6 is from "Contour Boxplots: A Method for Characterizing Uncertainty in Feature Sets from Simulation Ensembles"
		number of breaks that more than 1/6 grid points have can be though as threshold.
	*/
	int epsilon = 1;				// breaks less than epsilon can be considered as be sandwiched
	while (true) {	// break when find a proper epsilon
		int nCount = 0;
		for (size_t i = 0; i < nMatrixLength; i++)
		{
			if (arrMatrix[i] < epsilon) nCount++;
		}
		if (nCount > dbThreshold) break;
		else epsilon++;
	}
	// 3.claculate band depth
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		_pSetBandDepth[i] = 0;
		for (size_t j = 0; j < mSet; j++) {
			if (arrMatrix[i*mSet + j] < epsilon)_pSetBandDepth[i]++;
		}
	}
	delete[] arrMatrix;
}

const double* FeatureSet::GetMedian() { return _pData->GetData(_nMedianIndex); }

void FeatureSet::generateContourImp(const QList<ContourLine>& contourMin, const QList<ContourLine>& contourMax, QList<UnCertaintyArea*>& areas) {
	ContourStripGenerator generator;
	generator.Generate(areas, contourMin, contourMax, _nWidth, _nHeight);
}

void FeatureSet::CalculateSDF(double dbIsoValue,int nEnsembleLen,int nWidth,int nHeight
	,DataField* pData
	,double* pSDF
	,QList<QList<ContourLine>> listContour) 
{
	//qDebug() << "CalculateSDF";
	for (size_t l = 0; l < nEnsembleLen; l++) {
		// calculate sdf
		const double* arrData = pData->GetData(l);
		double* arrSDF = pSDF + l * nWidth*nHeight;
		QList<ContourLine> contour=listContour[l];
		for (size_t i = 0; i < nHeight; i++)
		{
			for (size_t j = 0; j < nWidth; j++)
			{
				double dbMinDis = nHeight*nWidth;
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
				if (arrData[i*nWidth + j] < dbIsoValue) dbMinDis = -dbMinDis;	// sign
				arrSDF[i*nWidth + j] = dbMinDis;
				//if (l == 0) qDebug() << dbMinDis;
			}
		}
	}
	//qDebug() << "~CalculateSDF";
}

/*
	calculate pca using diverse grids
*/
void FeatureSet::calculatePCA() {
	// 1.set parameter
	int mI = _nDiverseCount;
	int mO = _nPCLen;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = createDiverseArray();

	// 3.pca
	MyPCA pca;
	pca.DoPCA(arrInput, _arrPC, n, mI, mO, true);

	delete[] arrInput;
}

/*
	calculate pca using whole grids
*/
void FeatureSet::calculatePCA_Whole() {
	// 1.set parameter
	int mI = _nGrids;
	int mO = _nPCLen;
	int n = _nEnsembleLen;
	// 2.allocate input and output buffer
	double* arrInput = _pSDF;

	// 3.pca
	MyPCA pca;

	unsigned long t1 = GetTickCount();

	pca.DoPCA(arrInput, _arrPC, n, mI, mO, true);
	unsigned long t2 = GetTickCount();
	qDebug() << "Time of calculatePCA_Whole: " << t2-t1;

	delete[] arrInput;
}

/*
	calculate PCA using MDS for the diverse grids
*/
void FeatureSet::calculatePCA_MDS() {

	// 2.allocate input buffer
	double* arrInput = createDiverseArray();

	int n = _nEnsembleLen;
	std::vector<VectorXd> points(n);
	for (size_t i = 0; i < n; i++)
	{
		points[i] = VectorXd(_nDiverseCount);
		for (size_t j = 0; j < _nDiverseCount; j++)
		{
			points[i](j) = arrInput[i*_nDiverseCount + j];
		}
	}

	MatrixXd D(n, n);
	for (unsigned i = 0; i < n; ++i)
	{
		for (unsigned j = i; j < n; ++j)
		{
			const double d = (points[i] - points[j]).norm();
			D(i, j) = d;
			D(j, i) = d;
		}
	}

	// Compute metric MDS (embedding into a 2-dimensional space)
	const MatrixXd X = mathtoolbox::ComputeClassicalMds(D, _nPCLen);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < _nPCLen; j++)
		{
			_arrPC[i*_nPCLen + j] = X(j, i);
		}
	}

	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(0, i);
	}
	qDebug() << "======================";
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(1, i);
	}
	delete[] arrInput;
}

/*
	calculate PCA using MDS
	using SDFs.
*/
void FeatureSet::calculatePCA_MDS_Whole() {

	int n = _nEnsembleLen;

	// 1.build vectors
	std::vector<VectorXd> points(n);
	for (size_t l = 0; l < n; l++)
	{
		points[l] = VectorXd(_nHeight*_nWidth);

		for (size_t i = 0; i < _nHeight; i++)
		{
			for (size_t j = 0; j < _nWidth; j++)
			{
				points[l](i*_nWidth + j) = _pSDF[l*_nWidth*_nHeight + i * _nWidth + j];
			}

		}
	}

	// 2.calculate distance matrix of vectors
	MatrixXd D(n, n);
	for (unsigned i = 0; i < n; ++i)
	{
		for (unsigned j = i; j < n; ++j)
		{
			const double d = (points[i] - points[j]).norm();
			D(i, j) = d;
			D(j, i) = d;
		}
	}

	// 3.Compute metric MDS (embedding into a 2-dimensional space)
	const MatrixXd X = mathtoolbox::ComputeClassicalMds(D, _nPCLen);
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < _nPCLen; j++)
		{
			_arrPC[i*_nPCLen + j] = X(j, i);
		}
	}
}

/*
	calculate PCA usinng MDS and density
*/
void FeatureSet::calculatePCA_MDS_Whole_Density() {
	int nClusterWidth = _nWidth;// / 4;
	int nStartX = 0;
	//nClusterWidth = 10;
	//nStartX = 20;
	qDebug() << "Width:" << _nWidth;
	int n = _nEnsembleLen;
	std::vector<VectorXd> points(n);
	for (size_t l = 0; l < n; l++)
	{
		points[l] = VectorXd(_nHeight*nClusterWidth*2);
		//X
		for (size_t i = 0; i < _nHeight; i++)
		{
			for (size_t j = 0; j < nClusterWidth; j++)
			{
				points[l](i*nClusterWidth + j) = _pICDVX[l*_nWidth*_nHeight + i * _nWidth + nStartX + j];
			}
		}
		//Y
		for (size_t i = 0; i < _nHeight; i++)
		{
			for (size_t j = 0; j < nClusterWidth; j++)
			{
				points[l](_nHeight*nClusterWidth+i*nClusterWidth + j) = _pICDVY[l*_nWidth*_nHeight + i * _nWidth + nStartX + j];
			}
		}
	}

	MatrixXd D(n, n);
	for (unsigned i = 0; i < n; ++i)
	{
		for (unsigned j = i; j < n; ++j)
		{
			const double d = (points[i] - points[j]).norm();
			D(i, j) = d;
			D(j, i) = d;
		}
	}

	// Compute metric MDS (embedding into a 2-dimensional space)
	const MatrixXd X = mathtoolbox::ComputeClassicalMds(D, _nPCLen);
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < _nPCLen; j++)
		{
			_arrPC[i*_nPCLen + j] = X(j, i);
		}
	}
	/*
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(0, i);
	}
	qDebug() << "======================";
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(1, i);
	}
	*/
}

/*
	calculate PCA using mutual information
*/
void FeatureSet::calculatePCA_MutualInformation() {
	// calculate the range of sdf;
	double dbMax = -10000;
	double dbMin = 10000;
	for (size_t i = 0,length=_nEnsembleLen*_nGrids; i < length; i++)
	{
		if (_pSDF[i] > dbMax) dbMax = _pSDF[i];
		if (_pSDF[i] < dbMin) dbMin = _pSDF[i];
	}

	double dbSimilarityMax = calculateSimilarity(0, 0, dbMin, dbMax);

	MatrixXd D(_nEnsembleLen, _nEnsembleLen);
	for (unsigned i = 0; i < _nEnsembleLen; ++i)
	{
		for (unsigned j = i; j < _nEnsembleLen ; ++j)
		{
			const double d = dbSimilarityMax-calculateSimilarity(i, j, dbMin, dbMax);
			D(i, j) = d;
			D(j, i) = d;
		}
	}

	// Compute metric MDS (embedding into a 2-dimensional space)
	const MatrixXd X = mathtoolbox::ComputeClassicalMds(D, _nPCLen);
	for (size_t i = 0; i < _nEnsembleLen; i++)
	{
		for (size_t j = 0; j < _nPCLen; j++)
		{
			_arrPC[i*_nPCLen + j] = X(j, i);
		}
	}
}

/*
	calculate pca using MDS by set differences
*/
void FeatureSet::calculatePCA_MDS_Set() {
	int n = _nEnsembleLen;

	MatrixXd D(n, n);
	for (unsigned i = 0; i < n; ++i)
	{
		for (unsigned j = i; j < n; ++j)
		{
			int dis = 0;
			for (size_t k = 0; k < _nGrids; k++)
			{
				if (_pUpperSet[i*_nGrids + k] != _pUpperSet[j*_nGrids + k])
					dis++;
			}
			D(i, j) = dis;
			D(j, i) = dis;
		}
	}

	// Compute metric MDS (embedding into a 2-dimensional space)
	const MatrixXd X = mathtoolbox::ComputeClassicalMds(D, _nPCLen);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < _nPCLen; j++)
		{
			_arrPC[i*_nPCLen + j] = X(j, i);
		}
	}

	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(0, i);
	}
	qDebug() << "======================";
	for (size_t i = 0; i < n; i++)
	{
		qDebug() << X(1, i);
	}
}

/*
	calculate diverse grids and diverse count
	using UpperSet
*/
void FeatureSet::calculateDiverse() {
	_nDiverseCount = 0;
	for (size_t i = 0; i < _nGrids; i++)
	{
		bool bDiverse = false;

		for (size_t l = 1; l < _nEnsembleLen; l++)
		{
			if (_pUpperSet[l*_nGrids + i] != _pUpperSet[(l - 1)*_nGrids + i]) {
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

/*
	cluster using diverse grids
*/
void FeatureSet::doClustering() {
	// 1.create cluster instance
	CLUSTER::Clustering* pClusterer = new CLUSTER::AHCClustering();
	int nN = _nEnsembleLen;			// number of data items
	int nK = 1;

	// 2.clustering
	bool bUseDiversity = false;		// whether using diversity
	if (bUseDiversity) {
		// 2.cluster
		int nM = _nDiverseCount;		// dimension
		double* arrBuf = createDiverseArray();

		pClusterer->DoCluster(nN, nM, nK, arrBuf, _arrLabels, _arrMergeSource, _arrMergeTarget);

		delete[] arrBuf;
	}
	else {
		int nM = _nGrids;
		pClusterer->DoCluster(nN, nM, nK, _pSDF, _arrLabels, _arrMergeSource, _arrMergeTarget);
	}

	// 3.reset labels
	resetLabels();

	// 4.release cluster instance
	delete pClusterer;
}

/*
	cluster using PCA
*/
void FeatureSet::doPCAClustering() {

	// 2.cluster
	CLUSTER::Clustering* pClusterer = new CLUSTER::AHCClustering();
	int nN = _nEnsembleLen;			// number of data items
	int nM = _nPCLen;		// dimension
	int nK = 1;						// clusters
	double* arrBuf = _arrPC;

	pClusterer->DoCluster(nN, nM, nK, arrBuf, _arrLabels,_arrMergeSource,_arrMergeTarget);
	
	resetLabels();

	// 5.release the resouse
	delete pClusterer;
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
			ContourGenerator::GetInstance()->Generate(getSDF(l), contour, 0, _nWidth, _nHeight);
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
			ContourGenerator::GetInstance()->Generate(getSDF(l), contour, 0, _nWidth, _nHeight);
			_listContourSDF[l]=contour;
		}
	}

	// 5.release the buffer
	delete[] arrRecoveredBuf;

	delete[] arrOutput;
	delete[] arrInput;

}

// 感觉这个方法实现的还是不太对啊
void FeatureSet::calculateICD() {
	qDebug() << "calculateICD" << _dbIsoValue;
	double dbMax = 0;

	// 0.allocate resource
	double* pMeans = new double[_nGrids];
	double* pVari2 = new double[_nGrids];
	// 1.calculate means and variances
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
		pMeans[i] = dbMean;
		pVari2[i] = dbVar;
	}
	// 2.calculate probability
	int nWidth = (_nWidth - 1) * _nDetailScale + 1;
	int nHeight = (_nHeight - 1) * _nDetailScale + 1;
	for (size_t i = 0; i < nHeight; i++)
	{
		for (size_t j = 0; j < nWidth; j++)
		{
			int i0 = i / _nDetailScale;
			int j0 = j / _nDetailScale;
			int i1 = i0 + 1;
			int j1 = j0 + 1;
			int id = i % _nDetailScale;
			int jd = j % _nDetailScale;

			double dbi1 = id / (double)((_nDetailScale - 1) > 0 ? (_nDetailScale - 1) : 1);
			double dbi0 = 1.0 - dbi1;
			double dbj1 = jd / (double)((_nDetailScale - 1) > 0 ? (_nDetailScale - 1) : 1);
			double dbj0 = 1.0 - dbj1;

			double dbM00 = pMeans[i0*_nWidth + j0];
			double dbM01 = i1 < _nHeight ? pMeans[i1*_nWidth + j0] : 0;
			double dbM10 = j1 < _nWidth ? pMeans[i0*_nWidth + j1] : 0;
			double dbM11 = (i1 < _nHeight && j1 < _nWidth) ? pMeans[i1*_nWidth + j1] : 0;

			double dbV00 = pVari2[i0*_nWidth + j0];
			double dbV01 = i1 < _nHeight ? pVari2[i1*_nWidth + j0] : 0;
			double dbV10 = j1 < _nWidth ? pVari2[i0*_nWidth + j1] : 0;
			double dbV11 = (i1 < _nHeight && j1 < _nWidth) ? pVari2[i1*_nWidth + j1] : 0;

			double dbMean = dbM00 * dbi0*dbj0
				+ dbM01 * dbi1*dbj0
				+ dbM10 * dbi0*dbj1
				+ dbM11 * dbi1*dbj1;

			double dbVar2 = dbV00 * dbi0*dbj0
				+ dbV01 * dbi1*dbj0
				+ dbV10 * dbi0*dbj1
				+ dbV11 * dbi1*dbj1;

			double dbBias = _dbIsoValue - dbMean;
			double dbVar = sqrt(dbVar2);		// adjust the effect of variance
			//qDebug() << dbMean << "\t" << dbVar2 << "\t" << i << "\t" << j;
			//_pICD[i*nWidth + j] = pow(Ed, -dbBias * dbBias / 2 * dbVari);
			
			//_pICD[i*nWidth + j] = pow(Ed, -dbBias * dbBias / (2 * dbVar2)) /(sqrt(PI2d*dbVar2));

			_pICD[i*nWidth + j] = pow(Ed, -dbBias * dbBias / (2 * dbVar2)) / ((sqrt(PI2d))*dbVar);
			//if (_pICD[i*nWidth + j] > dbMax) {
			//	dbMax = _pICD[i*nWidth + j];
			//	qDebug() << "var:" << dbVar;
			//	qDebug() << "dbBias:" << dbBias;
			//}



		}
	}
	qDebug() << "Max of ICD: " << dbMax;
	// free resource
	delete pMeans;
	delete pVari2;
}

void FeatureSet::buildICD_LineKernel() {
	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICD_LineKernel[i] = 0;
	}
	int nRadius = 10;
	double dbH = 1.0;
	double dbMax = 0;
	for (size_t l = 0; l < _nEnsembleLen; l++) {
		for each (QList<ContourLine> contours in _listContour)
		{
			for each (ContourLine contour in contours)
			{
				for each (QPointF pt in contour._listPt)
				{
					double x = pt.x();
					double y = pt.y();
					for (int iX=x- nRadius;iX<x+ nRadius;iX++)
					{
						if (iX < 0 || iX >= _nWidth) continue;
						for (int iY = y - nRadius; iY < y + nRadius; iY++)
						{
							if (iY < 0 || iY >= _nHeight) continue;
							double dbX = iX - x;
							double dbY = iY - y;
							double dbDis2 = dbX*dbX + dbY * dbY;
							double dbV = pow(Ed, -dbDis2 / 2);
							_pICD_LineKernel[iY*_nWidth + iX] += dbV;
							if (_pICD_LineKernel[iY*_nWidth + iX] > dbMax)
								dbMax = _pICD_LineKernel[iY*_nWidth + iX];
						}
					}
				}
			}
		}
	}

	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICD_LineKernel[i] /= (dbMax + 1);
	}
}

void FeatureSet::buildICD_Vector() {
	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICDX[i] = 0;
		_pICDY[i] = 0;
		_pICDZ[i] = 0;
	}
	int nRadius = 40;
	double dbH = 1.0;
	double dbMaxX = 0;
	double dbMaxY = 0;
	double dbMaxZ = 0;
	for (size_t l = 0; l < _nEnsembleLen; l++) {
		for each (QList<ContourLine> contours in _listContour)
		{
			for each (ContourLine contour in contours)
			{
				int nLen = contour._listPt.length();
				for (size_t i = 1; i < nLen; i++)
				{
					double x1 = contour._listPt[i - 1].x();
					double y1 = contour._listPt[i - 1].y();
					double x2 = contour._listPt[i].x();
					double y2 = contour._listPt[i].y();

					for (int iX = x1 - nRadius; iX < x1 + nRadius; iX++)
					{
						if (iX < 0 || iX >= _nWidth) continue;
						for (int iY = y1 - nRadius; iY < y1 + nRadius; iY++)
						{
							if (iY < 0 || iY >= _nHeight) continue;
							double dbX = iX - x1;
							double dbY = iY - y1;
							double dbDis2 = dbX * dbX + dbY * dbY;
							double dbV = pow(Ed, -dbDis2 / 2);
							int nIndex = iY * _nWidth + iX;
							double dbVX = x2 - x1;
							double dbVY = y2 - y1;
							_pICDZ[nIndex] += abs(_pICDX[nIndex] * dbVY + _pICDY[nIndex] * dbVX)*dbV;

							_pICDX[nIndex] += dbVX * dbV;
							_pICDY[nIndex] += dbVY * dbV;

						//	if (_pICDX[nIndex] > dbMaxX) dbMaxX = _pICDX[nIndex];
						//	if (_pICDY[nIndex] > dbMaxY) dbMaxY = _pICDY[nIndex];
							if (_pICDZ[nIndex] > dbMaxZ) dbMaxZ = _pICDZ[nIndex];
						}
					}
				}
			}
		}
	}


	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICDX[i] = sqrt(_pICDX[i] * _pICDX[i] + _pICDY[i] * _pICDY[i]);
		if (_pICDX[i] > .1) _pICDY[i] = _pICDZ[i] / _pICDX[i];
		else _pICDY[i] = 0;

		if (_pICDX[i] > dbMaxX) dbMaxX = _pICDX[i];
		if (_pICDY[i] > dbMaxY) dbMaxY = _pICDY[i];
	}
	for (size_t i = 0; i < _nGrids; i++)
	{
		_pICDX[i] /= (dbMaxX + 1);
		_pICDY[i] /= (dbMaxY + 1);
		_pICDZ[i] /= (dbMaxZ + 1);
	}
}

double Kernel(double x) {
	return 1 / SQRTP2d * pow(Ed, -x * x / 2);
}
void FeatureSet::buildICDV() {

	for (size_t i = 0, len = _nEnsembleLen*_nGrids; i < len; i++)
	{
		_pICDVX[i] = 0;
		_pICDVY[i] = 0;
	}
	int nRadius = 50;
	double dbH = 1;
	double dbMaxX = 0;
	double dbMaxY = 0;
	double dbMaxZ = 0;
	for (size_t l = 0; l < _nEnsembleLen; l++) {
		QList<ContourLine> contours  = _listContour[l];

		for each (ContourLine contour in contours)
		{
			int nLen = contour._listPt.length();
			for (size_t i = 1; i < nLen; i++)
			{
				double x1 = contour._listPt[i - 1].x();
				double y1 = contour._listPt[i - 1].y();
				double x2 = contour._listPt[i].x();
				double y2 = contour._listPt[i].y();
				double x0 = (x1 + x2) / 2.0;
				double y0 = (y1 + y2) / 2.0;

				for (int iX = x1 - nRadius; iX < x1 + nRadius; iX++)
				{
					if (iX < 0 || iX >= _nWidth) continue;
					for (int iY = y1 - nRadius; iY < y1 + nRadius; iY++)
					{
						if (iY < 0 || iY >= _nHeight) continue;
						double dbX = iX - x0;
						double dbY = iY - y0;
						double dbDis = sqrt(dbX * dbX + dbY * dbY);
						double dbV = Kernel(dbDis / dbH) / dbH;

						double dbVX = x2 - x1;
						double dbVY = y2 - y1;
						int nIndex = l * _nGrids + iY * _nWidth + iX;

						_pICDVX[nIndex] += dbVX * dbV;
						_pICDVY[nIndex] += dbVY * dbV;
					}
				}
			}
		}
		for (size_t i = 0; i < _nGrids; i++)
		{
			_pICDVX[l * _nGrids + i] *= 20;
			_pICDVY[l * _nGrids + i] *= 20;
		}
	}
	for (size_t i = 0, len = _nEnsembleLen * _nGrids; i < len; i++)
	{
		_pICDVW[i] = sqrt(_pICDVX[i] * _pICDVX[i] + _pICDVY[i] * _pICDVY[i]);
	}

}


void FeatureSet::SetClustersLen(int nClustersLen) {
	_nClusters = nClustersLen;
	resetLabels();
}

void FeatureSet::resetLabels() {
	// initialize labels
	for (size_t i = 0; i < _nEnsembleLen; i++) {
		_arrLabels[i] = i;
	}

	// merge
	for (int i = 0; i < _nEnsembleLen-_nClusters; i++)
	{
		int nSource = _arrMergeSource[i];
		int nTarget = _arrMergeTarget[i];

		// update source and target
		for (int j = 0; j < _nEnsembleLen; j++)
		{
			if (_arrLabels[j] == nSource) {
				_arrLabels[j] = nTarget;
			}
		}
	}

	// align
	CLUSTER::Clustering::AlignLabels(_arrLabels, _nEnsembleLen);

}

void FeatureSet::generateRandomIndices() {
	// generate a random indices list
	int arrIndices[50];
	for (int i = 0; i < 50; i++) arrIndices[i] = i;
	for (int l = 50; l > 0; l--)
	{
		// random choose an index
		int nIndex = rand() / (double)RAND_MAX*l;
		_arrRandomIndices[l - 1] = arrIndices[nIndex];
		arrIndices[nIndex] = arrIndices[l - 1];
	}
}
void FeatureSet::generateSequentialIndices() {
	for (int i = 0; i < 50; i++) _arrRandomIndices[i] = i;
}

void FeatureSet::generateCentralIndices() {
	int nBias = (_nEnsembleLen - _nTestSamples) / 2;
	for (size_t i = 0; i < _nTestSamples; i++)
	{
		_arrRandomIndices[i] = i + nBias;
	}
}

void FeatureSet::measureInfoLose() {
	srand(time(NULL));
	// generate indices list
	switch (g_nTestTarget)
	{
	case 1:
		generateRandomIndices();
		break;
	case 2:
		generateCentralIndices();
		break;
	case 3:
		generateSequentialIndices();
		break;
	default:
		break;
	}

	int nMK = 5000000;
	double dbLose = 0;
	for (size_t i = 0; i < nMK; i++)
	{
		double x = 1.0 * rand() / RAND_MAX * _nWidth;
		double y = 1.0 * rand() / RAND_MAX * _nHeight ;
		dbLose += measureDis(x, y);
	}

	qDebug()<<"Samples: " << _nTestSamples << "measureInfoLose: " << dbLose;
}

double FeatureSet::measureDis(double x, double y) {

	double dbPO = getUpperProperty(x, y, _pSDF, _nEnsembleLen);
	double dbPS;

	switch (g_nTestTarget)
	{
	case 0:
		dbPS=getUpperProperty(x, y, _pResampledSDF, _nTestSamples);
		break;
	case 1:
		dbPS = getUpperProperty(x, y, _pSDF, _nTestSamples, true);
		break;
	case 2:
		dbPS = getUpperProperty(x, y, _pSortedSDF, _nTestSamples, true);
		break;
	case 3:
		dbPS = getUpperProperty(x, y, _pSortedSDF, _nTestSamples, true);
		break;
	default:
		break;
	}

	//qDebug() << dbPO << dbPS << abs(dbPO - dbPS);
	return abs(dbPO - dbPS);
}

double FeatureSet::getUpperProperty(double x, double y, const double* pSDF, int nLen, bool bUseIndices) {
	int nX = x;
	int nY = y;
	double dbX = x - nX;
	double dbY = y - nY;

	int nCount = 0;

	for (size_t l = 0; l < nLen; l++)
	{
		if (getFieldValue(nX, nY, dbX, dbY, pSDF + (bUseIndices? _arrRandomIndices[l] :l) * _nGrids) > 0) nCount++;
	}
	return nCount / (double)nLen;
}

double FeatureSet::getFieldValue(int nX, int nY, double dbX, double dbY, const double* pField) {
	double dbV00 = pField[nY*_nWidth + nX];
	double dbV01 = pField[(nY+1)*_nWidth + nX];
	double dbV10 = pField[nY*_nWidth + nX+1];
	double dbV11 = pField[(nY+1)*_nWidth + nX+1];

	//qDebug() << dbV00 << dbV01 << dbV10 << dbV11;

	return dbX * dbY*dbV00
		+ dbX * (1 - dbY)*dbV01
		+ (1 - dbX) * dbY*dbV10
		+ (1 - dbX) * (1 - dbY)*dbV11;
}


void FeatureSet::measureInfoLoseMI() {
	srand(time(NULL));
	// generate indices list
	generateRandomIndices();
	//generateCentralIndices();
	//generateUniformIndices();



	const int nBin = 128;
	int arrBin[nBin][nBin];
	int arrBinX[nBin];
	int arrBinY[nBin];

	for (size_t i = 0; i < nBin; i++) {
		arrBinX[i] = 0;
		arrBinY[i] = 0;
		for (size_t j = 0; j < nBin; j++) arrBin[i][j] = 0;
	}

	// calculate bin
	int nMK = 5000000;
	for (size_t i = 0; i < nMK; i++)
	{
		double x = 1.0 * rand() / RAND_MAX * _nWidth;
		double y = 1.0 * rand() / RAND_MAX * _nHeight;
		double dbPO = getUpperProperty(x, y, _pSDF, _nEnsembleLen);
		double dbPS = getUpperProperty(x, y, _pSDF, _nTestSamples, g_nTestTarget ==0);
		int nX = dbPO * 100;
		int nY = dbPS * 100;
		if (nX == 100) nX--;
		if (nY == 100) nY--;
		arrBin[nX][nY]++;
		arrBinX[nX]++;
		arrBinY[nY]++;
	}

	// calculate mutual information
	double dbHx = 0;
	double dbHy = 0;
	double dbHxy = 0;
	for (size_t i = 0; i < nBin; i++)
	{
		if (arrBinX[i]) {
			double dbPx = arrBinX[i] / (double)nMK;
			dbHx += dbPx * log(dbPx);
		}
		if (arrBinY[i]) {
			double dbPy = arrBinY[i] / (double)nMK;
			dbHy += dbPy * log(dbPy);
		}
		for (size_t j = 0; j < nBin; j++)
		{
			if (arrBin[i][j]) {
				double dbPxy = arrBin[i][j] / (double)nMK;
				dbHxy += dbPxy * log(dbPxy);
			}
		}
	}
	double dbLose = -(dbHx + dbHy - dbHxy);

	qDebug() << "Samples: " << _nTestSamples << "measureInfoLose: " << dbLose;
}













