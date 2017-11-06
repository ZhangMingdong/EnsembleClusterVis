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
	int _nW;						// width
	int _nH;						// height
	int _nL;						// number of ensemble members
	int _nSmooth = 0;				// level of smooth
	int _nEOF = 0;					// level of EOF
	double* _gridEOF[g_nEOFLen];	// eof fields
public:
	// get the l'th layer
	const double* GetLayer(int l);
	// get variance, nSmooth: level of smooth
	const double* GetVari(int nSmooth=0);
	const double* GetEOF(int nSeq=0);
	const double* GetDipValue();
	const double* GetMean();
	const double* GetUMax();
	const double* GetUMin();
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
	void DoStatistic();
	// generate clustered data using the labels
	void GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData);

	// perform eof analysis
	void DoEOF();
private:
	// smooth the variance to level nSmooth
	void smoothVar(int nSmooth);
};

