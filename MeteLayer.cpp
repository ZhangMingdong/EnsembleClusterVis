#include "MeteLayer.h"
#include "LayerLayout.h"
#include "MeteModel.h"
#include "VarAnalysis.h"
#include "ColorMap.h"

#include <iostream>
#include <math.h>

ILayerCB* MeteLayer::_pCB = NULL;

MeteLayer::MeteLayer() :_bShow(true), _pModel(NULL)
{
}

MeteLayer::~MeteLayer()
{
}

void MeteLayer::DrawLayer(DisplayStates states){
	if (_bShow) {
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glPushMatrix();
		this->draw(states);
		glPopMatrix();
		glPopAttrib();
	}
}

void MeteLayer::InitLayer(const LayerLayout* pLayout, double fScaleW, double fScaleH){
	_pLayout = pLayout;

	_fScaleW = fScaleW;
	_fScaleH = fScaleH;

	_gllist = glGenLists(1);	// generate the display lists
	this->init();
}

// set color according to group index
void MeteLayer::SetGroupColor(int nIndex) {
	if (nIndex<20)
	{
		glColor3f(ColorMap::GetCategory20D(nIndex, 0), ColorMap::GetCategory20D(nIndex, 1), ColorMap::GetCategory20D(nIndex, 2));
	}
	return;
	// old code
	switch (nIndex)
	{
	case 0:
		glColor3f(1, 0, 0);
		break;
	case 1:
		glColor3f(0, 1, 0);
		break;
	case 2:
		glColor3f(0, 0, 1);
		break;
	case 3:
		glColor3f(0, 1, 1);
		break;
	case 4:
		glColor3f(1, 0, 1);
		break;
	default:
		break;
	}
}

void MeteLayer::OnSelectVar(int nIndex) {
	if (_pModel)
	{
		std::vector<double> vecVar = _pModel->GetVariance();
		double dbVarThreshold = vecVar[nIndex*vecVar.size() / VarAnalysis::_nVarBars];
		_pModel->SetVarThreshold(dbVarThreshold);
		ReloadTexture();
	}
}

void MeteLayer::drawCircle(double x, double y, double dbRadius, bool bFill) {
	int nSegs = 20;
	if (bFill) glBegin(GL_POLYGON);
	else glBegin(GL_LINE_LOOP);

	for (size_t i = 0; i < nSegs; i++)
	{
		double dbAngle = 2.0*M_PI*i / nSegs;
		double dbX = x + dbRadius*cos(dbAngle);
		double dbY = y + dbRadius*sin(dbAngle);
		std::cout << dbX << ": "<<dbY << std::endl;
		glVertex3f(dbX, dbY, 0);
	}

	glEnd();
}