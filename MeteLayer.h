/*
	Meteorology display layer
	Mingdong
	2017/05/05
*/
#pragma once

#include <QGLWidget>

class MeteModel;

// callback interface
class ILayerCB{
public:
	virtual void DrawText(char* pText, double fX, double fY) = 0;
};

class MeteLayer
{
public:
	static ILayerCB* _pCB;
public:
	struct DisplayStates 
	{
		bool _bShowBackground = false;
		bool _bShowContourLineMin = false;
		bool _bShowContourLineMax = false;
		bool _bShowContourLineMean = false;
		bool _bShowUnionE = false;
		bool _bShowGradientE = false;
		bool _bShowContourLine = false;
		bool _bShowBeliefEllipse = false;
	};
public:
	MeteLayer();
	virtual ~MeteLayer();
public:
	// interface for render
	void DrawLayer(DisplayStates states);
	// interface for initialization
	void InitLayer(double fLeft, double fRight, double fTop, double fBottom, double fScaleW, double fScaleH);
	// set model for the layer
	void SetModel(MeteModel* pModel){ _pModel = pModel; };
	// show or hide
	void Show(bool bShow){ _bShow = bShow; }
	// reload texture
	virtual void ReloadTexture() {};
	// brushing
	virtual void Brush(int nLeft, int nRight, int nTop, int nBottom) {};

	// select an index in variance widget, total range is 100
	virtual void OnSelectVar(int nIndex) {}
private:
	// implementation for render
	virtual void draw(DisplayStates states) = 0;
	// implementation for initialization
	virtual void init() = 0;


protected:
	MeteModel* _pModel;
	GLuint _gllist;                           // display index
	GLuint _gllistC;                           // display index for clustering
	GLuint _gllistG;                           // display index for gradient
	double _fLeft;
	double _fRight;
	double _fBottom;
	double _fTop;
	double _fScaleW;
	double _fScaleH;
	bool _bShow;
};

