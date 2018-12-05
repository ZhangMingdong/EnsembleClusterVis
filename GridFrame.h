#pragma once
class GridFrame
{
public:
	GridFrame();
	~GridFrame();

protected:
	int _nWidth;
	int _nHeight;
	int _nGrids;								// _nWidth*_nHeight
	int _nFocusX;
	int _nFocusY;
	int _nFocusW;
	int _nFocusH;
	int _nFocusGrids;							//_nFocusW*_nFocusH
	int _nWest;
	int _nEast;
	int _nSouth;
	int _nNorth;
	int _nFocusWest;
	int _nFocusEast;
	int _nFocusSouth;
	int _nFocusNorth;
	int _nEnsembleLen;						// number of ensemble members
public:
	int GetW() { return _nWidth; }
	int GetH() { return _nHeight; }
	int GetFocusW() { return _nFocusW; }
	int GetFocusH() { return _nFocusH; }
	int GetFocusX() { return _nFocusX; }
	int GetFocusY() { return _nFocusY; }
	int GetWest() { return _nWest; }
	int GetEast() { return _nEast; }
	int GetSouth() { return _nSouth; }
	int GetNorth() { return _nNorth; }
	int GetFocusWest() { return _nFocusWest; }
	int GetFocusEast() { return _nFocusEast; }
	int GetFocusSouth() { return _nFocusSouth; }
	int GetFocusNorth() { return _nFocusNorth; }
	int GetEnsembleLen() { return _nEnsembleLen; }
};

