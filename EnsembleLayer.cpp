#include "EnsembleLayer.h"
#include "MeteModel.h"

#include "ColorMap.h"
#include "def.h"
#include "LayerLayout.h"

#include "SpatialCluster.h"


#include <QDebug>


// for tess

/// 初始化用户自定义的程序变量

/** 环绕规则 */
GLdouble TessProperty[5] =
{
	GLU_TESS_WINDING_ODD,        /**< 环绕数为奇数 */
	GLU_TESS_WINDING_NONZERO,    /**< 环绕数为非0 */
	GLU_TESS_WINDING_POSITIVE,   /**< 环绕数为正数 */
	GLU_TESS_WINDING_NEGATIVE,   /**< 环绕数为负数 */
	GLU_TESS_WINDING_ABS_GEQ_TWO /**< 环绕数绝对值大于等于2 */
};
GLuint nProperty = 0;              /**< 环绕规则索引 */

/** gluTessCallback注册的回调函数 */
void CALLBACK beginCallback(GLenum which)
{
	glBegin(which);
}

void CALLBACK errorCallback(GLenum errorCode)
{
	const GLubyte *string;
	///输出错误信息
	string = gluErrorString(errorCode);
	// 	fprintf(stderr, "Tessellation Error: %s\n", string);
	// 	exit(0);
}

void CALLBACK endCallback(void)
{
	glEnd();
}

void CALLBACK vertexCallback(GLvoid *vertex)
{

	GLdouble* pt;
	GLubyte red, green, blue;
	pt = (GLdouble*)vertex;
	/** 随即产生颜色值 */
	red = (GLubyte)rand() & 0xff;
	green = (GLubyte)rand() & 0xff;
	blue = (GLubyte)rand() & 0xff;
	glColor3ub(red, green, blue);
	glVertex3dv(pt);

}

/** 用于处理检测轮廓线交点，并决定是否合并顶点，
新创建的顶点最多可以是4个已有顶点线性组合，这些定点坐标存储在data中
其中weight为权重，weight[i]的总合为1.0
*/
void CALLBACK combineCallback(GLdouble coords[3],
	GLdouble *vertex_data[4],
	GLfloat weight[4],
	GLdouble **dataOut)
{
	GLdouble *vertex;
	int i;
	/** 分配存储新顶点的内存 */
	vertex = (GLdouble *)malloc(6 * sizeof(GLdouble));

	/** 存储坐标值 */
	vertex[0] = coords[0];
	vertex[1] = coords[1];
	vertex[2] = coords[2];

	/** 通过插值计算RGB颜色值 */
	for (i = 3; i < 6; i++)
		vertex[i] = weight[0] * vertex_data[0][i] +
		weight[1] * vertex_data[1][i] +
		weight[2] * vertex_data[2][i] +
		weight[3] * vertex_data[3][i];
	*dataOut = vertex;
}

// ~for tess
EnsembleLayer::EnsembleLayer()
{
}

EnsembleLayer::~EnsembleLayer()
{
	if (_colorbarTexture)
	{
		delete[]_colorbarTexture;
	}
	
}

void EnsembleLayer::draw(DisplayStates states){
	// border
	glBegin(GL_LINE_LOOP);
	glVertex2f(_pLayout->_dbLeft, _pLayout->_dbBottom);
	glVertex2f(_pLayout->_dbRight, _pLayout->_dbBottom);
	glVertex2f(_pLayout->_dbRight, _pLayout->_dbTop);
	glVertex2f(_pLayout->_dbLeft, _pLayout->_dbTop);
	glEnd();



	int nClusterIndex = -1;	// -1 means show all clusters

	// the x radius and y radius of the map in drawing space
	double biasX = (_pModel->GetWest() + _pModel->GetEast()) / 2.0 *_fScaleW;
	double biasY = (_pModel->GetSouth() + _pModel->GetNorth()) / 2.0*_fScaleH;

	// the focused x radius and y radius of the map in drawing space
	double biasFocusX = (_pModel->GetFocusWest() + _pModel->GetFocusEast()) / 2.0 *_fScaleW;
	double biasFocusY = (_pModel->GetFocusSouth() + _pModel->GetFocusNorth()) / 2.0*_fScaleH;

	double scaleFocusX = (_pModel->GetFocusEast() - _pModel->GetFocusWest()) / 360.0;
	double scaleFocusY = (_pModel->GetFocusNorth() - _pModel->GetFocusSouth()) / 180.0;

	// draw background
	if (states._bShowBackground){
		glPushMatrix();
		glEnable(GL_TEXTURE_2D);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glBindTexture(GL_TEXTURE_2D, _uiTexID[0]);

		glTranslatef(biasFocusX, biasFocusY, 0);			// 位置偏移
		glScalef(scaleFocusX, scaleFocusY, 0);				// 改变尺寸

		glBegin(GL_QUADS);

		glTexCoord2f(0.0f, 0.0f); glVertex2f(_pLayout->_dbLeft, _pLayout->_dbBottom);
		glTexCoord2f(1.0f, 0.0f); glVertex2f(_pLayout->_dbRight, _pLayout->_dbBottom);
		glTexCoord2f(1.0f, 1.0f); glVertex2f(_pLayout->_dbRight, _pLayout->_dbTop);
		glTexCoord2f(0.0f, 1.0f); glVertex2f(_pLayout->_dbLeft, _pLayout->_dbTop);

		glEnd();

		glDisable(GL_TEXTURE_2D);
		glPopMatrix();

		// color bar
		drawColorBar();
	}

	glTranslatef(biasX, biasY, 0);				// 移动到指定位置
	
	glTranslatef(-(_pModel->GetW()-1)*_fScaleW/2, -(_pModel->GetH() - 1)*_fScaleH / 2, 0);		// 移动到中间

	glScalef(_fScaleW, _fScaleH, 0);


	if (states._bShowContourLineMin)
	{
		glColor4f(1, 1, 0.0, 1.0);
		drawContourLine(_pModel->GetContourMin());
	}
	if (states._bShowContourLineMax){
		glColor4f(0.0, 1.0, 1.0, 1.0);
		drawContourLine(_pModel->GetContourMax());
	}
	if (states._bShowContourLineMean)
	{
		glColor4f(1.0, 0.0, 1.0, .5);
		drawContourLine(_pModel->GetContourMean());
	}
	if (states._bShowContourLine)
	{
		bool bBrush = false;
		if (bBrush) {
			QList<QList<ContourLine>> contoursBrushed = _pModel->GetContourBrushed();
			QList<QList<ContourLine>> contoursNotBrushed = _pModel->GetContourNotBrushed();
			glColor4f(.8, 0.2, 0.0, .8);
			for (int i = 0; i < contoursBrushed.size(); i++)
			{
				drawContourLine(contoursBrushed[i]);
			}
			glColor4f(.2, 0.8, 0.0, .2);
			for (int i = 0; i < contoursNotBrushed.size(); i++)
			{
				drawContourLine(contoursNotBrushed[i]);
			}
		}
		else {
			QList<QList<ContourLine>> contours = _pModel->GetContour();
			const ClusterResult* pCR=_pModel->GetClusterResultOfFocusedRegion();
			for (int i = 0; i < contours.size(); i++)
			{
//				SetGroupColor(pCR->_arrLabels[i]);
				drawContourLine(contours[i]);
			}
		}		
	}
	if (states._bShowGradientE)
	{
		for (size_t i = 0; i < g_nEnsembles; i++)
		{
			glColor4f(1.0,1.0,0.0,0.06);
			glCallList(_gllistG + i * 3 + 1);
		}
	}
	// union
	if (states._bShowUnionE)
	{
		double fTransparency = .2;
		// render each area
		if (g_bShowUncertaintyOnly) {

			glColor4f(0, 1, 0, .5);
			glCallList(_gllist + 1);
		}
		else
			for (size_t j = 0; j < 3; j++)
			{
				glColor4f(ColorMap::GetRGB(j, 0), ColorMap::GetRGB(j, 1), ColorMap::GetRGB(j, 2), fTransparency);
				glCallList(_gllist + j);
			}
	}
}

void EnsembleLayer::ReloadTexture() {
//	_dataTexture = _pModel->generateTexture();
//	_dataTexture = _pModel->generateTextureGridCluster();
//	_dataTexture = _pModel->generateTextureRange(0);
//	_dataTexture = _pModel->generateTextureMean();
//	_dataTexture = _pModel->generateTextureDiscreteSummary();
	_dataTexture = _pModel->generateTextureNew();
//	_dataTexture = _pModel->generateTextureSDF();

	glBindTexture(GL_TEXTURE_2D, _uiTexID[0]);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _pModel->GetFocusW(), _pModel->GetFocusH(), 0, GL_RGBA, GL_UNSIGNED_BYTE, _dataTexture);

	generateColorBarTexture();
}

void EnsembleLayer::init(){
	// 0.generate texture
	glGenTextures(2, &_uiTexID[0]);

	// 1.load texture
	ReloadTexture();

	// 2.generate color bar
	generateColorBarTexture();

	// 3.initialize tess
	buildTess();
}

void EnsembleLayer::drawContourLine(const QList<ContourLine>& contours){

	for each (ContourLine contour in contours)
	{
		glBegin(GL_LINE_STRIP);
		for each (QPointF pt in contour._listPt)
		{
			double x = pt.x();// *_fScaleW;
			double y = pt.y();// *_fScaleH;
			glVertex2f(x, y);
		}
		glEnd();
	}
}

void EnsembleLayer::tessSegmentation(GLuint gllist, QList<UnCertaintyArea*> areas){
	// create buffer for the uncertainty areas	
	int nAreaLen = areas.size();										// how many area
	GLdouble**** contourBuffer = new GLdouble***[areas.size()];
	// for each area
	for (int i = 0; i < nAreaLen; i++)
	{
		int nContourLen = areas[i]->_listEmbeddedArea.size() + 1;		// how many contour
		contourBuffer[i] = new GLdouble**[nContourLen];
		// first contour
		int nLen = areas[i]->_contour._listPt.size();
		contourBuffer[i][0] = new GLdouble*[nLen];
		for (int j = 0; j < nLen; j++)
		{
			contourBuffer[i][0][j] = new GLdouble[3];
			contourBuffer[i][0][j][0] = areas[i]->_contour._listPt[j].x();
			contourBuffer[i][0][j][1] = areas[i]->_contour._listPt[j].y();
			contourBuffer[i][0][j][2] = 0;
		}
		// other contours
		for (int k = 0; k < nContourLen - 1; k++)
		{
			int nLen = areas[i]->_listEmbeddedArea[k]->_contour._listPt.size();
			contourBuffer[i][k + 1] = new GLdouble*[nLen];
			for (int j = 0; j < nLen; j++)
			{
				contourBuffer[i][k + 1][j] = new GLdouble[3];
				contourBuffer[i][k + 1][j][0] = areas[i]->_listEmbeddedArea[k]->_contour._listPt[j].x();
				contourBuffer[i][k + 1][j][1] = areas[i]->_listEmbeddedArea[k]->_contour._listPt[j].y();
				contourBuffer[i][k + 1][j][2] = 0;
			}
		}
	}
	/*
		state==-1: negative
		state==0:  uncertain
		state==1:  positive
	*/
	for (int state = -1; state < 2; state++)
	{
		glNewList(gllist + 1 + state, GL_COMPILE);
		gluTessProperty(_tobj, GLU_TESS_WINDING_RULE, TessProperty[nProperty]);
		gluTessBeginPolygon(_tobj, NULL);
		for (int i = 0; i < nAreaLen; i++)
		{
			// 		if (i != 1) continue;
			if (areas[i]->_nState == state)
			{
				int nContourLen = areas[i]->_listEmbeddedArea.size() + 1;
				for (int k = 0; k < nContourLen; k++)
				{
					int nLen;
					if (k == 0) nLen = areas[i]->_contour._listPt.size();
					else nLen = areas[i]->_listEmbeddedArea[k - 1]->_contour._listPt.size();
					gluTessBeginContour(_tobj);
					for (int j = 0; j < nLen; j++)
					{
						gluTessVertex(_tobj, contourBuffer[i][k][j], contourBuffer[i][k][j]);
					}
					gluTessEndContour(_tobj);
				}
			}
		}
		gluTessEndPolygon(_tobj);
		glEndList();
	}

	for (int i = 0; i < nAreaLen; i++)
	{
		int nContourLen = areas[i]->_listEmbeddedArea.size() + 1;		// how many contour
		// first contour
		int nLen = areas[i]->_contour._listPt.size();
		for (int j = 0; j < nLen; j++)
		{
			delete contourBuffer[i][0][j];
		}
		for (int k = 0; k < nContourLen - 1; k++)
		{
			int nLen = areas[i]->_listEmbeddedArea[k]->_contour._listPt.size();
			for (int j = 0; j < nLen; j++)
			{
				delete[]contourBuffer[i][k + 1][j];
			}
			delete[]contourBuffer[i][k + 1];
		}
		delete[]contourBuffer[i];
	}

	delete[]contourBuffer;
}

void EnsembleLayer::Brush(int nLeft, int nRight, int nTop, int nBottom) {
	_pModel->Brush(nLeft, nRight, nTop, nBottom);
}

void EnsembleLayer::drawColorBar() {
	// border
	glColor3f(0, 0, 0);		// black
	glBegin(GL_LINE_LOOP);
	glVertex2f(_pLayout->_dbColorBarLeft, _pLayout->_dbBottom);
	glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbBottom);
	glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbTop);
	glVertex2f(_pLayout->_dbColorBarLeft, _pLayout->_dbTop);
	glEnd();
	MeteModel::enumBackgroundFunction bgFun = _pModel->GetBgFunction();
	switch (bgFun)
	{
	case MeteModel::bg_mean:
	case MeteModel::bg_vari:
	case MeteModel::bg_vari_smooth:
	case MeteModel::bg_dipValue:
	case MeteModel::bg_EOF:
	{
		// scale
		ColorMap* colormap = bgFun== MeteModel::bg_EOF?ColorMap::GetInstance(ColorMap::CP_EOF): ColorMap::GetInstance();
		int nLen = colormap->GetLength();
		int nStep = colormap->GetStep();
		int nMin = colormap->GetMin();
		for (int i = 0; i < nLen; i++)
		{
			glBegin(GL_LINES);
			glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbBottom + (_pLayout->_dbTop - _pLayout->_dbBottom)*i / (nLen - 1));
			glVertex2f(_pLayout->_dbColorBarRight-.01, _pLayout->_dbBottom + (_pLayout->_dbTop - _pLayout->_dbBottom)*i / (nLen - 1));
			glEnd();
			// draw text

			char buf[10];
			sprintf_s(buf, "%d", nMin+i*nStep);
			_pCB->DrawText(buf, _pLayout->_dbColorBarRight + .02, _pLayout->_dbBottom + (_pLayout->_dbTop - _pLayout->_dbBottom)*i / (nLen - 1) - .01);
		}
		// colors
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, _uiTexID[1]);
		glBegin(GL_QUADS);
		glTexCoord2f(0.1f, 0.0f); glVertex2f(_pLayout->_dbColorBarLeft, _pLayout->_dbBottom);
		glTexCoord2f(0.9f, 0.0f); glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbBottom);
		glTexCoord2f(0.9f, 1.0f); glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbTop);
		glTexCoord2f(0.1f, 1.0f); glVertex2f(_pLayout->_dbColorBarLeft, _pLayout->_dbTop);
		glEnd();
		glDisable(GL_TEXTURE_2D);
	}
		break;
	case MeteModel::bg_cluster:
	{
		const std::vector<UncertaintyRegion> vecR = _pModel->GetUncertaintyRegions();
		int nRegions = _pModel->GetUncertaintyAreas();
		int nLen = 0;
		for (size_t i = 0; i < nRegions; i++)
		{
			nLen += vecR[i]._nArea;
		}
		// draw

		char buf[20];
		double dbTop = _pLayout->_dbTop;		
		int nLenAccumulate = 0;
		for (size_t i = 0; i < nRegions; i++)
		{
			nLenAccumulate+= vecR[i]._nArea;
			double dbBottom = _pLayout->_dbTop-_pLayout->_dbHeight*nLenAccumulate / nLen;
			double dbMid = (dbTop + dbBottom) / 2.0;
			SetRegionColor(i);
			glBegin(GL_QUADS);
			glVertex2f(_pLayout->_dbColorBarLeft, dbBottom);
			glVertex2f(_pLayout->_dbColorBarRight, dbBottom);
			glVertex2f(_pLayout->_dbColorBarRight, dbTop);
			glVertex2f(_pLayout->_dbColorBarLeft, dbTop);
			glEnd();

			glColor3f(0, 0, 0);
			sprintf_s(buf, "%.2f%%", vecR[i]._nArea*100.0/nLen);
			_pCB->DrawText(buf, _pLayout->_dbColorBarRight + .02, dbMid);

			dbTop = dbBottom;
		}
	}
		break;
	case MeteModel::bg_varThreshold:
	case MeteModel::bg_dipValueThreshold:
	{
		int nLen = _pModel->GetGridLength();
		int nTLen = _pModel->GetThresholdedGridLength();
		double dbMid = _pLayout->_dbTop-_pLayout->_dbHeight*nTLen / nLen;
		glColor3f(ColorMap::GetThresholdColorD(0), ColorMap::GetThresholdColorD(1), ColorMap::GetThresholdColorD(2));
		glBegin(GL_QUADS);
		glVertex2f(_pLayout->_dbColorBarLeft, dbMid);
		glVertex2f(_pLayout->_dbColorBarRight, dbMid);
		glVertex2f(_pLayout->_dbColorBarRight, _pLayout->_dbTop);
		glVertex2f(_pLayout->_dbColorBarLeft, _pLayout->_dbTop);
		glEnd();

		glColor3f(0, 0, 0);
		int nPercentage = nTLen * 100 / nLen;
		char buf[10];
		sprintf_s(buf, "%d%%", nPercentage);
		_pCB->DrawText(buf, _pLayout->_dbColorBarRight + .02, (_pLayout->_dbTop + dbMid) / 2.0);


		sprintf_s(buf, "%d%%", 100-nPercentage);
		_pCB->DrawText(buf, _pLayout->_dbColorBarRight + .02, (_pLayout->_dbBottom + dbMid) / 2.0);

	}
		break;
	default:
		break;
	}
}

// generate the texture for the color bar
void EnsembleLayer::generateColorBarTexture() {
	// 1.generate color bar data
	MeteModel::enumBackgroundFunction bgFun = _pModel->GetBgFunction();
	ColorMap* colormap = bgFun == MeteModel::bg_EOF ? ColorMap::GetInstance(ColorMap::CP_EOF) : ColorMap::GetInstance();
	int nLen = colormap->GetLength();
	int nStep = colormap->GetStep();
	int nMin = colormap->GetMin();
	int nColorBarW = 1;
	int nColorBarH = (nLen - 1) * 10 + 1;
	_colorbarTexture = new GLubyte[nColorBarH*nColorBarW * 4];
	for (int i = 0; i < nColorBarH; i++)
	{
//		MYGLColor color = colormap->GetColor(i / 10.0*nStep);
		MYGLColor color = colormap->GetColor(nMin+i*nStep/10.0);
		for (size_t j = 0; j < nColorBarW; j++)
		{
			int nIndex = nColorBarW*i + j;
			_colorbarTexture[4 * nIndex + 0] = color._rgb[0];
			_colorbarTexture[4 * nIndex + 1] = color._rgb[1];
			_colorbarTexture[4 * nIndex + 2] = color._rgb[2];
			_colorbarTexture[4 * nIndex + 3] = (GLubyte)255;
		}
	}
	// 2.build texture
	glBindTexture(GL_TEXTURE_2D, _uiTexID[1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nColorBarW, nColorBarH, 0, GL_RGBA, GL_UNSIGNED_BYTE, _colorbarTexture);
	// ~build texture
}

// build the tess for uncertainty regions
void EnsembleLayer::buildTess() {
	_gllist = glGenLists(3);					// generate the display lists
	_gllistG = glGenLists(g_gradient_l * 3);	// generate the display lists


	_tobj = gluNewTess();
	gluTessCallback(_tobj, GLU_TESS_VERTEX,
		(void(__stdcall*)())glVertex3dv);
	gluTessCallback(_tobj, GLU_TESS_BEGIN,
		(void(__stdcall*)())beginCallback);
	gluTessCallback(_tobj, GLU_TESS_END,
		(void(__stdcall*)())endCallback);
	gluTessCallback(_tobj, GLU_TESS_ERROR,
		(void(__stdcall*)())errorCallback);


	// 4.tess the areas
	tessSegmentation(_gllist, _pModel->GetUncertaintyArea());
	/*
	QList<QList<UnCertaintyArea*> > listAreas = _pModel->GetUncertaintyAreaG();
	for (size_t i = 0, len = listAreas.length(); i < len; i++)
	{
	tessSegmentation(_gllistG + i * 3, listAreas[i]);
	}
	*/
}

	