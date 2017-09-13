/*
	ensemble result display layer
	Mingdong
	2017/05/05
*/
#pragma once
#include "MeteLayer.h"
#include "ContourGenerator.h"
#include <QGLWidget>
#include <gl/GLU.h>
#include "GLFont.h"
class UnCertaintyArea;
class EnsembleLayer :
	public MeteLayer
{
public:
	EnsembleLayer();
	virtual ~EnsembleLayer();
protected:
	virtual void draw(DisplayStates states);
	virtual void init();
	// select an index in variance widget, total range is 100
	virtual void OnSelectVar(int nIndex);
private:
	// draw a contour line
	void drawContourLine(const QList<ContourLine>& contours);

	// tessellation the area segmentation, generate three display list start from _gllist
	void tessSegmentation(GLuint gllist, QList<UnCertaintyArea*> areas);

	// draw the chart of variance
	void drawVarChart();

	// draw the pca points
	void drawPCAPoints();

	// draw the cluster bars
	void drawClusterBars();
private:
	// truth texture, generated from truth data
	GLubyte* _dataTexture;
	// texture of the color bar, 160*2
	GLubyte* _colorbarTexture;
	// texture id,1-data,2-colorbar
	GLuint texID[2];

	GLUtesselator *_tobj;                    /**< 网格化对象指针 */

	// threshold of the variance
	double _dbVarThreshold = 1.5;
	
public:
	// reload texture
	virtual void ReloadTexture();
	// brushing
	virtual void Brush(int nLeft, int nRight, int nTop, int nBottom);

};

