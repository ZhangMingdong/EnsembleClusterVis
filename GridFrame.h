#pragma once
#include <QObject>
class GridFrame:public QObject
{
public:
	GridFrame();
	~GridFrame();

protected:
	int _nWidth;
	int _nHeight;
	int _nGrids;								// _nWidth*_nHeight

	int _nWest;
	int _nEast;
	int _nSouth;
	int _nNorth;

	int _nEnsembleLen;						// number of ensemble members

protected:
	void sortBuf(const double* pS, double* pD);	// sort source to target


	void resampleBuf(const double* pSortedBuf, double* pResampledBuf, int nResampleLen);	// resample buffer
	void resampleBuf_C(const double* pUnSortedBuf, const int* pLabels, double* pResampledBuf, int nClusters);	// resample buffer for each clusters


public:
	int GetW() { return _nWidth; }
	int GetH() { return _nHeight; }
	int GetWest() { return _nWest; }
	int GetEast() { return _nEast; }
	int GetSouth() { return _nSouth; }
	int GetNorth() { return _nNorth; }
	int GetEnsembleLen() { return _nEnsembleLen; }
};

