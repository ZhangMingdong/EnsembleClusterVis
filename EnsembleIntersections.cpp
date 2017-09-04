#include "EnsembleIntersections.h"



Intersection::Intersection(){

}

Intersection::Intersection(double fMin, double fMax, int nOverlap) :_fMax(fMax), _fMin(fMin), _nOverlap(nOverlap){

}

void EnsembleIntersections::Intersect(double fMin, double fMax){
	int nLen = _listIntersection.length();
	if (_listIntersection.empty()){
		// insert the first intersection
		_listIntersection.append(Intersection(fMin, fMax, 1));
	}
	else if (fMin >= _listIntersection[nLen - 1]._fMax){
		// append at the end, also insert an empty intersection
		if (fMin > _listIntersection[nLen - 1]._fMax)
		{
			_listIntersection.append(Intersection(_listIntersection[nLen - 1]._fMax, fMin, 0));
		}
		_listIntersection.append(Intersection(fMin, fMax, 1));
	}
	else if (fMax <= _listIntersection[0]._fMin)
	{
		// insert at the beginning, also the empty intersection between
		if (fMax < _listIntersection[0]._fMin)
		{
			_listIntersection.insert(_listIntersection.begin(), Intersection(fMax, _listIntersection[0]._fMin, 0));
		}
		_listIntersection.insert(_listIntersection.begin(), Intersection(fMin, fMax, 1));
	}
	else{
		// calculate new intersection
		int nMinL = -1;		// left index of min
		int nMaxL = -1;		// left index of max
		for (int i = 0; i < nLen; i++)
		{
			if (_listIntersection[i]._fMin < fMin)
			{
				nMinL = i;
			}
			if (_listIntersection[i]._fMin < fMax)
			{
				nMaxL = i;
			}
		}
		// if fMax is bigger than the max value of last intersection, set it to nLen+1
		if (fMax>_listIntersection[nLen-1]._fMax)
		{
			nMaxL = nLen;
		}
		if (nMinL == nMaxL){
			// both in one intersection, divide this intersection into 3 piece
			// piece 1
			_listIntersection.insert(_listIntersection.begin() + nMinL
				, Intersection(_listIntersection[nMinL]._fMin, fMin, _listIntersection[nMinL]._nOverlap));
			// piece 2
			_listIntersection.insert(_listIntersection.begin() + nMinL + 1
				, Intersection(fMin, fMax, _listIntersection[nMinL]._nOverlap + 1));
			// piece 3
			_listIntersection[nMinL + 2]._fMin = fMax;
		}
		else{
			// increment the overlap
			for (int i = nMinL + 1; i < nMaxL; i++)
			{
				_listIntersection[i]._nOverlap++;
			}
			// handle the right end
			if (nMaxL == nLen){
				_listIntersection.insert(_listIntersection.end(), Intersection(_listIntersection[nLen - 1]._fMax, fMax, 1));
			}
			else{
				_listIntersection.insert(_listIntersection.begin() + nMaxL
					, Intersection(_listIntersection[nMaxL]._fMin, fMax, _listIntersection[nMaxL]._nOverlap));
				_listIntersection[nMaxL + 1]._fMin = fMax;
			}
			// handle the left end
			if (nMinL == -1){
				_listIntersection.insert(_listIntersection.begin(), Intersection(fMin, _listIntersection[0]._fMin, 1));
			}
			else{
				_listIntersection.insert(_listIntersection.begin() + nMinL
					, Intersection(_listIntersection[nMinL]._fMin, fMin, _listIntersection[nMinL]._nOverlap));
				_listIntersection[nMinL + 1]._fMin = fMin;
				_listIntersection[nMinL + 1]._nOverlap++;
			}

		}


	}
}

void EnsembleIntersections::GetMostOverlapped(double& fMin, double& fMax){
	int nOverlap = 0;
	for each (Intersection intersection in _listIntersection)
	{
		if (intersection._nOverlap > nOverlap)
		{
			nOverlap = intersection._nOverlap;
			fMin = intersection._fMin;
			fMax = intersection._fMax;
		}
		else if (intersection._nOverlap == nOverlap){
			if (intersection._fMin<fMin)
			{
				fMin = intersection._fMin;
			}
			if (intersection._fMax>fMax)
			{
				fMax = intersection._fMax;
			}
		}
	}
}

void EnsembleIntersections::GetMostOverlapped(double& fMin, double& fMax, int nThreshold){
	fMin = 100;
	fMax = -100;
	for each (Intersection intersection in _listIntersection)
	{
		if (intersection._nOverlap > nThreshold)
		{
			if (intersection._fMin<fMin)
			{
				fMin = intersection._fMin;
			}
			if (intersection._fMax>fMax)
			{
				fMax = intersection._fMax;
			}
		}
	}
}
