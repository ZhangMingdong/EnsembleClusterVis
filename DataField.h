#pragma once
#include <QList>
#include "def.h"
/*
	represent a data field
	author: Mingdong
	2017/06/26
	functionality: provide grid related informations
	2017/10/24
*/
class DataField
{
public:
	DataField(int w,int h,int l);
	~DataField();
private:	
	double* _pBuf;					// store the data	
	double* _pDipBuf;				// store the dip data	
	double* _gridVari;				// variance of each grid among different ensemble results	
	double* _gridVV;				// variance of variance of each grid among different ensemble results	
	double* _gridVarSmooth[10];		// smoothed variance
	double* _gridMean;				// mean of each grid among different ensemble results
	double* _gridUMax;				// maximum of union of each grid among different ensemble results
	double* _gridUMin;				// minimum of union of each grid among different ensemble results


	double* _gridHalfMax;				// maximum of half valid ensemble results
	double* _gridHalfMin;				// minimum of half valid ensemble results
	double* _gridValidMax;				// maximum of valid ensemble results
	double* _gridValidMin;				// minimum of valid ensemble results
	int _nMedianIndex = -1;				// index of median of the contour

	double* _pSortedBuf;			// data sorted at each grid point
	double* _pSDF;					// data of signed distance function
	double* _pSortedSDF;			// data of signed distance function

	bool * _pSet;					// set state of the grid point given iso-value
	int *_pSetBandDepth;			// sBandDepth
	int *_pRegionType;				// region type of each member.0:outlier,1-100%,2-50%.
	int _nOutlierThreshold = 1;		// threshold for outliars


	int _nW;						// width
	int _nH;						// height
	int _nGrids;					// width*height
	int _nL;						// number of ensemble members
	int _nSmooth = 0;				// level of smooth
	int _nEOF = 0;					// level of EOF
	double* _gridEOF[g_nEOFLen];	// eof fields
public:
	// get the l'th layer
	const double* GetLayer(int l);
	const double* GetSortedLayer(int l);
	// get variance, nSmooth: level of smooth
	const double* GetVari(int nSmooth=0);
	const double* GetEOF(int nSeq=0);
	const double* GetDipValue();
	const double* GetMean();
	const double* GetMedian();
	const double* GetUMax();
	const double* GetUMin();
	const double* GetValidMax() { return _gridValidMax; };
	const double* GetValidMin() { return _gridValidMin; };
	const double* GetHalfMax() { return _gridHalfMax; };
	const double* GetHalfMin() { return _gridHalfMin; };
	double* GetEditableLayer(int l);
	// set the data value at a given position
	void SetData(int l, int bias, double dbValue);
	void SetData(int l, int y,int x, double dbValue);
	// set the dip value
	void SetDipValue(int bias, double dbValue);
	// get the data value at a given position
	double GetData(int l, int bias);
	double GetData(int l, int r, int c);
	const double* GetData() { return _pBuf; }
	const double* GetData(int l) { return _pBuf+l* _nGrids; }
	double* GetSDF() { return _pSDF; }
	double* GetSDF(int l) { return _pSDF + l * _nGrids; }
	double* GetSortedSDF(int l) { return _pSortedSDF + l * _nGrids; }
	bool* GetSet(int l) { return _pSet + l * _nGrids; }
	void DoStatistic();
	void CalculateSet(double dbIsoValue);
	int GetDepth(int l) { return _pSetBandDepth[l]; }
	int GetRegionType(int l) { return _pRegionType[l]; }
	/*
		generate clustered data using the labels
		params:
			listClusterLens: length of each cluster
			arrLabels: cluster labels of each cluster member
			arrData: return generate data for each cluster
		written before, used 2017/11/07
	*/
	void GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData);

	// perform eof analysis
	void DoEOF();

private:
	// smooth the variance to level nSmooth
	void smoothVar(int nSmooth);

private:
	void sortBuf(const double* pS, double* pD);	// sort source to target
	void buildSortedBuf();			// build the sorted buf from buf
public:
	void BuildSortedSDF();

};

