#include "ColorMap.h"


MYGLColor::MYGLColor(){
	_rgb[0] = _rgb[1] = _rgb[2] = 0;
}

MYGLColor::MYGLColor(GLubyte r, GLubyte g, GLubyte b){
	_rgb[0] = r;
	_rgb[1] = g;
	_rgb[2] = b;
}


ColorMap* ColorMap::GetInstance(Enum_Color_Pallet cp){
	if (!_bInitialized)
	{
		_bInitialized = true;
		for (size_t i = 0; i < CP_Length; i++)
		{
			_pInstance[i] = NULL;
		}
	}
	if (!_pInstance[cp])
	{
		_pInstance[cp] = new ColorMap();
		switch (cp)
		{
		case ColorMap::CP_Perception:
			_pInstance[cp]->_nLen = 9;
			_pInstance[cp]->_nStep = 2;					// for mean
										// 	_nStep = 7;					// for variance
			_pInstance[cp]->_pValues = new double[_pInstance[cp]->_nLen];
			_pInstance[cp]->_pColors = new MYGLColor[_pInstance[cp]->_nLen];

			// 	_pValues[0] = 0;
			// 	_pValues[1] = 2;
			// 	_pValues[2] = 5;
			// 	_pValues[3] = 10;
			// 	_pColors[0] = MYGLColor(255,255,255);
			// 	_pColors[1] = MYGLColor(0, 0, 255);
			// 	_pColors[2] = MYGLColor(0, 255, 0);
			// 	_pColors[3] = MYGLColor(255, 0, 0);

			for (int i = 0; i < _pInstance[cp]->_nLen; i++)
			{
				_pInstance[cp]->_pValues[i] = i * _pInstance[cp]->_nStep;
			}
			_pInstance[cp]->_pColors[0] = MYGLColor(165, 206, 227);
			_pInstance[cp]->_pColors[1] = MYGLColor(81, 154, 165);
			_pInstance[cp]->_pColors[2] = MYGLColor(86, 177, 70);
			_pInstance[cp]->_pColors[3] = MYGLColor(249, 141, 141);
			_pInstance[cp]->_pColors[4] = MYGLColor(239, 105, 68);
			_pInstance[cp]->_pColors[5] = MYGLColor(255, 135, 16);
			_pInstance[cp]->_pColors[6] = MYGLColor(178, 149, 200);
			_pInstance[cp]->_pColors[7] = MYGLColor(198, 180, 153);
			_pInstance[cp]->_pColors[8] = MYGLColor(176, 88, 40);
			break;
		case ColorMap::CP_RB:
			_pInstance[cp]->_nLen = 3;
			_pInstance[cp]->_nStep = 10;					// for mean
														// 	_nStep = 7;					// for variance
			_pInstance[cp]->_pValues = new double[_pInstance[cp]->_nLen];
			_pInstance[cp]->_pColors = new MYGLColor[_pInstance[cp]->_nLen];


			_pInstance[cp]->_pValues[0] = -10;
			_pInstance[cp]->_pValues[1] = 0;
			_pInstance[cp]->_pValues[2] = 10;
			_pInstance[cp]->_pColors[0] = MYGLColor(0, 0, 255);
			_pInstance[cp]->_pColors[1] = MYGLColor(255, 255, 255);
			_pInstance[cp]->_pColors[2] = MYGLColor(255, 0, 0);
		case ColorMap::CP_RainBow:
			_pInstance[cp]->_nLen = 8;
			_pInstance[cp]->_nStep = 2;					// for mean
															// 	_nStep = 7;					// for variance
			_pInstance[cp]->_pValues = new double[_pInstance[cp]->_nLen];
			_pInstance[cp]->_pColors = new MYGLColor[_pInstance[cp]->_nLen];


			_pInstance[cp]->_pValues[0] = -7;
			_pInstance[cp]->_pValues[1] = -5;
			_pInstance[cp]->_pValues[2] = -3;
			_pInstance[cp]->_pValues[3] = -1;
			_pInstance[cp]->_pValues[4] = 1;
			_pInstance[cp]->_pValues[5] = 3;
			_pInstance[cp]->_pValues[6] = 5;
			_pInstance[cp]->_pValues[7] = 7;

			_pInstance[cp]->_pColors[0] = MYGLColor(87, 0, 255);
			_pInstance[cp]->_pColors[1] = MYGLColor(43, 0, 255);
			_pInstance[cp]->_pColors[2] = MYGLColor(0, 0, 255);
			_pInstance[cp]->_pColors[3] = MYGLColor(0, 255, 255);
			_pInstance[cp]->_pColors[4] = MYGLColor(0, 255, 0);
			_pInstance[cp]->_pColors[5] = MYGLColor(255, 255, 0);
			_pInstance[cp]->_pColors[6] = MYGLColor(255, 165, 0);
			_pInstance[cp]->_pColors[7] = MYGLColor(255, 0, 0);
			break;
		case ColorMap::CP_Length:
			break;
		default:
			break;
		}
		
	}
	return _pInstance[cp];
}

bool ColorMap::_bInitialized = false;

ColorMap* ColorMap::_pInstance[CP_Length];

ColorMap::ColorMap()
{

}


ColorMap::~ColorMap()
{
	delete _pColors;
	delete _pValues;
}

MYGLColor ColorMap::GetColor(double fValue){
	int i = 1;
	while (fValue > _pValues[i] && i < _nLen) i++;


	if (i == _nLen) return _pColors[_nLen - 1];
	else{
		MYGLColor color0 = _pColors[i - 1];
		MYGLColor color1 = _pColors[i];
		return interpolateColor(color0, color1, fValue - _pValues[i - 1], _pValues[i] - fValue);
	}
}
MYGLColor ColorMap::interpolateColor(MYGLColor color0, MYGLColor color1, double fBias0, double fBias1){
	MYGLColor result;
	for (int i = 0; i < 3;i++)
	{
		result._rgb[i] = (color0._rgb[i] * fBias1 + color1._rgb[i] * fBias0) / (fBias0 + fBias1);
	}
	return result;
}