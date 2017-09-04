#include "MeteLayer.h"


ILayerCB* MeteLayer::_pCB = NULL;

MeteLayer::MeteLayer() :_bShow(true), _pModel(NULL)
{
}

MeteLayer::~MeteLayer()
{
}

void MeteLayer::DrawLayer(DisplayStates states){
	if (_bShow)
		this->draw(states);
}

void MeteLayer::InitLayer(double fLeft, double fRight, double fTop, double fBottom, double fScaleW, double fScaleH){
	_fLeft = fLeft;
	_fRight = fRight;
	_fBottom = fBottom;
	_fTop = fTop;

	_fScaleW = fScaleW;
	_fScaleH = fScaleH;

	_gllist = glGenLists(1);	// generate the display lists
	this->init();
}
