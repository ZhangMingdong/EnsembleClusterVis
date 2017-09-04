#include "CostLineLayer.h"


#include <gl/GLU.h>
#include <iostream>
#include "../Shapefile/Shapefile.h"

using namespace std;

CostLineLayer::CostLineLayer()
{
}


CostLineLayer::~CostLineLayer()
{
}

void CostLineLayer::draw(DisplayStates states){

	drawCostline();
}
void CostLineLayer::init(){

	readCostline();
}



void CostLineLayer::drawCostline(){
	glColor3f(0, 0, 0);
	glCallList(_gllist);

	// draw the frame
	glBegin(GL_LINE_LOOP);
	glVertex2f(_fLeft, _fTop);
	glVertex2f(_fRight, _fTop);
	glVertex2f(_fRight, _fBottom);
	glVertex2f(_fLeft, _fBottom);

	glEnd();
}

void CostLineLayer::readCostline(){

	const char* path = "E:\\Meteorology\\coastline\\ne_10m_coastline\\ne_10m_coastline.shp";

	FILE* pShapefile = open_shapefile(path);

	if (pShapefile == 0) {
		return;
	}
	double fWidth = _fRight - _fLeft;
	double fHeight = _fTop - _fBottom;

	SFShapes* pShapes = read_shapes(pShapefile);

	glNewList(_gllist, GL_COMPILE);
	for (uint32_t x = 0; x < pShapes->num_records; ++x) {
		const SFShapeRecord* record = pShapes->records[x];
		SFPolyLine* polyline = get_polyline_shape(pShapefile, record);
		printf("Polyline %d: num_parts = %d, num_points = %d\n", x, polyline->num_parts, polyline->num_points);
		/* Do things with the polygon. */
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < polyline->num_points; i++)
		{
			double x = polyline->points[i].x*_fScaleW;
			double y = polyline->points[i].y*_fScaleH;
			glVertex2f(x, y);
		}
		glEnd();

		free_polyline_shape(polyline);
		polyline = 0;

		fflush(stdout);
	}
	glEndList();

	free_shapes(pShapes);
	close_shapefile(pShapefile);

	/*
	glNewList(_gllist + g_gradient_l * 3 + g_temperatureLen * 3 + g_temperatureLen * 3 + 2 + 1, GL_COMPILE);
	glBegin(GL_LINE_LOOP);
	glVertex2f(0.0, 0.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(1.0, 1.0);
	glEnd();
	glEndList();
	*/

	cout << "read cost line" << endl;

}