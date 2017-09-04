#pragma once
#include <QList>
/*
	represent a data field
	Mingdong
	before 2017/06/26
*/
class DataField
{
public:
	DataField(int w,int h,int l);
	~DataField();
private:	
	double* _pBuf;			// store the data	
	double* _gridVari;		// variance of each grid among different ensemble results	
	double* _gridMean;		// mean of each grid among different ensemble results
	double* _gridUMax;		// maximum of union of each grid among different ensemble results
	double* _gridUMin;		// minimum of union of each grid among different ensemble results
	int _nW;				// width
	int _nH;				// height
	int _nL;				// number of ensemble members
public:
	// get the l'th layer
	const double* GetLayer(int l);
	const double* GetVari();
	const double* GetMean();
	const double* GetUMax();
	const double* GetUMin();
	double* GetEditableLayer(int l);
	// set the data value at a given position
	void SetData(int l, int bias, double dbValue);
	void SetData(int l, int y,int x, double dbValue);
	// get the data value at a given position
	double GetData(int l, int bias);
	double GetData(int l, int r, int c);
	const double* GetData() { return _pBuf; }
	void DoStatistic();
	// generate clustered data using the labels
	void GenerateClusteredData(const QList<int> listClusterLens, const int* arrLabels, QList<DataField*>& arrData);
};

