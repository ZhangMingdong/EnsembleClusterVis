/*
	ensemble result display layer
	Mingdong
	2017/05/05
*/
#pragma once

#include "MeteLayer.h"
#include <QGLWidget>
#include <gl/GLU.h>

#include "BasicStruct.h"

class UnCertaintyArea;

class EnsembleLayer : public MeteLayer
{
public:
	EnsembleLayer();
	virtual ~EnsembleLayer();
protected:
	virtual void draw(DisplayStates states);
	virtual void init();
private:
	// draw a contour line
	void drawContourLine(const QList<ContourLine>& contours);

	// tessellation the area segmentation, generate three display list start from _gllist
	void tessSegmentation(GLuint gllist, QList<UnCertaintyArea*> areas);

	// draw the color bar
	void drawColorBar();

	// generate the texture for the color bar
	void generateColorBarTexture();

	// build the tess for uncertainty regions
	void buildTess();
private:
	// texture of the color bar, 160*2
	GLubyte* _colorbarTexture;
	// texture id:0-data,1-colorbar
	GLuint _uiTexID[2];
	// 网格化对象指针
	GLUtesselator *_tobj;  
	
public:
	// reload texture
	virtual void ReloadTexture();
	// brushing
	virtual void Brush(int nLeft, int nRight, int nTop, int nBottom);

private:
	void drawContourLineOutlier();
	void drawContourLineMin();
	void drawContourLineMax();
	void drawContourLineMean();
	void drawContourLineMedian();
	void drawContourLine();
	void drawContourLineSorted();
	void drawContourLineSortedSDF();
	void drawContourLineResampled();
	void drawContourLineDomainResampled();
	void drawContourLineSDF();
	void drawBand();
public:
	bool _arrIsoValueState[10] = { true,true,true,true,true,true,true,true,true,true };		// states of each isovalue
};

