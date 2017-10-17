// ============================================================
// used to generate color map
// Written by Mingdong
// 2017/05/07
// ============================================================
#pragma once
#include <QGLWidget>
// struct of the color
struct MYGLColor
{
	GLubyte _rgb[3];	// 0~255
public:
	MYGLColor();
	MYGLColor(GLubyte r, GLubyte g, GLubyte b);
};
class ColorMap
{
public:
	enum Enum_Color_Pallet
	{
		CP_Perception	// perception color, support value 0~16, using 9 colors
		, CP_RB			// red/blue color, support value -10~10, using 2 colors only
		, CP_RainBow		// rainbow color, support value -9~9, using 7 colors 
		, CP_12			// 12 colors
		, CP_Length		// length of color pallet
	};
public:
	static ColorMap* GetInstance(Enum_Color_Pallet cp = CP_Perception);
private:
	static ColorMap* _pInstance[CP_Length];
	// weather initialized
	static bool _bInitialized;
private:
	ColorMap();
public:
	// get a color given an value
	MYGLColor GetColor(double fValue);
	~ColorMap();
	int GetLength() { return _nLen; }
	int GetStep() { return _nStep; }
private:
	MYGLColor* _pColors;		// list of color
	double* _pValues;		// list of value
	int _nLen;				// length of color mapping
	int _nStep;
private:
	MYGLColor interpolateColor(MYGLColor color0, MYGLColor color1, double fBias0, double fBias1);

private:
	// color categoryt20 from d3
	static int s_arrCategoryColor20[20][3];

	// self defined 20 colors
	static double s_arrSelfDefinedColor20[20][3];

	// threshold color
	static double s_arrThreshold[3];
public:
	// get color category20: 0~255
	static int GetCategory20I(int nColorIndex, int nComponentIndex);
	// get color category20: 0.0~1.0
	static double GetCategory20D(int nColorIndex, int nComponentIndex);
	// get color of red, green, and blue
	static double GetRGB(int nColorIndex, int nComponentIndex);
	// get the color used for threshhold
	static int GetThresholdColorI(int nComponentIndex);
	// get the color used for threshold
	static double GetThresholdColorD(int nComponentIndex);

	// get 10 colors in category20: 0~255
	static int GetCategory10I(int nColorIndex, int nComponentIndex, int nBias = 0);
	// get 10 colors in category20: 0.0~1.0
	static double GetCategory10D(int nColorIndex, int nComponentIndex, int nBias = 0);

};

