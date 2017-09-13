#include "EnsembleLayer.h"
#include "MeteModel.h"

#include "ColorMap.h"
#include "def.h"
#include <QDebug>

// set color according to group index
void SetGroupColor(int nIndex) {
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
EnsembleLayer::EnsembleLayer() :_dataTexture(NULL)
{
}

EnsembleLayer::~EnsembleLayer()
{
	if (_dataTexture)
	{
		delete[]_dataTexture;
	}
	if (_colorbarTexture)
	{
		delete[]_colorbarTexture;
	}
	
}

void EnsembleLayer::draw(DisplayStates states){
	glPushAttrib(GL_LINE_BIT);
	// border
	glBegin(GL_LINE_LOOP);
	glVertex2f(_fLeft, _fBottom);
	glVertex2f(_fRight, _fBottom);
	glVertex2f(_fRight, _fTop);
	glVertex2f(_fLeft, _fTop);
	glEnd();


	//*
	double g_arrColors[9][3] = {
		{ 1, 0, 0 },			// R
		{ 0, 1, 0 },			// G
		{ 0, 0, 1 },			// B
		{ 1, 1, 0 },			// R
		{ 0, 1, 1 },			// G
		{ 1, 0, 1 },			// B
		{ 1, .5, 0 },			// R
		{ 0, 1, .5 },			// G
		{ .5, 0, 1 },			// B
	};


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

		glBindTexture(GL_TEXTURE_2D, texID[0]);

		glTranslatef(biasFocusX, biasFocusY, 0);			// 位置偏移
		glScalef(scaleFocusX, scaleFocusY, 0);				// 改变尺寸

		glBegin(GL_QUADS);

		glTexCoord2f(0.0f, 0.0f); glVertex2f(_fLeft, _fBottom);
		glTexCoord2f(1.0f, 0.0f); glVertex2f(_fRight, _fBottom);
		glTexCoord2f(1.0f, 1.0f); glVertex2f(_fRight, _fTop);
		glTexCoord2f(0.0f, 1.0f); glVertex2f(_fLeft, _fTop);

		glEnd();

		glPopMatrix();

		// color bar
		glBindTexture(GL_TEXTURE_2D, texID[1]);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f); glVertex2f(_fRight+.03, _fBottom);
		glTexCoord2f(1.0f, 0.0f); glVertex2f(_fRight+.06, _fBottom);
		glTexCoord2f(1.0f, 1.0f); glVertex2f(_fRight+.06, _fTop);
		glTexCoord2f(0.0f, 1.0f); glVertex2f(_fRight+.03, _fTop);
		glEnd();
		glDisable(GL_TEXTURE_2D);

		// border
		glBegin(GL_LINE_LOOP);
		glVertex2f(_fRight + .03, _fBottom);
		glVertex2f(_fRight + .06, _fBottom);
		glVertex2f(_fRight + .06, _fTop);
		glVertex2f(_fRight + .03, _fTop);
		glEnd();
		ColorMap* colormap = ColorMap::GetInstance();
		int nLen = colormap->GetLength();
		int nStep = colormap->GetStep();
		// scale
		for (int i = 1; i < nLen; i++)
		{
			glBegin(GL_LINES);
			glVertex2f(_fRight + .06, _fBottom + (_fTop - _fBottom)*i / (nLen-1));
			glVertex2f(_fRight + .05, _fBottom + (_fTop - _fBottom)*i / (nLen-1));
			glEnd();
			// draw text

			char buf[10];
			sprintf_s(buf, "%d", i*nStep);
// 			font.PrintText(buf, _fRight + .06, _fBottom + (_fTop - _fBottom)*i / 8);
			_pCB->DrawText(buf, _fRight + .063, _fBottom + (_fTop - _fBottom)*i / (nLen-1));
		}
	}


	glPushMatrix();

	glTranslatef(biasX, biasY, 0);				// 移动到指定位置

	/*
	if (!_pModel->GetFilter())
		glScalef(.5, .5, 1);						// 变换大小
	glScalef(2, 2, 0);
	*/

	glTranslatef(-(_pModel->GetW()-1)*_fScaleW/2, -(_pModel->GetH() - 1)*_fScaleH / 2, 0);		// 移动到中间


	glScalef(_fScaleW, _fScaleH, 0);

	// draw contour line
	if (states._bShowContourLineMin)
	{
		if (g_bClustering)
		{
			if (nClusterIndex>-1)	// show single cluster
			{
				glColor4f(.5, .5, 0.0, 1.0);
				drawContourLine(_pModel->GetContourMin(nClusterIndex));

			}
			else {

				double fOpacity = 1.0;
				for (size_t i = 0, length = _pModel->GetClusters(); i < length; i++)
				{
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity);
					drawContourLine(_pModel->GetContourMin(i));
				}
			}
		}
		else {
			glColor4f(1, 1, 0.0, 1.0);
			drawContourLine(_pModel->GetContourMin());
		}
	}
	if (states._bShowContourLineMax){
		if (g_bClustering)
		{
			if (nClusterIndex>-1)	// show single cluster
			{
				int i = 0;
				glColor4f(0.0, 1.0, 1.0, 1.0);
				drawContourLine(_pModel->GetContourMax(nClusterIndex));
			}
			else {
				double fOpacity = 1.0;
				for (size_t i = 0, length = _pModel->GetClusters(); i < length; i++)
				{
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity);
					drawContourLine(_pModel->GetContourMax(i));
				}

			}
		}
		else {
			glColor4f(0.0, 1.0, 1.0, 1.0);
			drawContourLine(_pModel->GetContourMax());
		}
	}
	if (states._bShowContourLineMean)
	{
		if (g_bClustering)
		{
			double fOpacity = 1.0;
			for (size_t i = 0,length=_pModel->GetClusters(); i < length; i++)
			{
				if (false) { // show two kinds of means
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity);
					drawContourLine(_pModel->GetContourMeanPCA(i));
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity / 2);
					drawContourLine(_pModel->GetContourMean(i));
				}
				else {
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity);
					drawContourLine(_pModel->GetContourMean(i));
				}
			}			
		}
		else {
			glColor4f(1.0, 0.0, 1.0, .5);
			drawContourLine(_pModel->GetContourMean());
		}
	}
	// show contour line
	if (states._bShowContourLine)
	{
		if (g_bClustering)
		{
			if (false)		// show first one
			{
				QList<QList<ContourLine>> contours = _pModel->GetContour();
				glColor4f(1.0, 0.0, 0.0, 1.0);
				drawContourLine(contours[0]);
			}
			else {
				float fOpacity = .5;
				for (size_t i = 0, length = _pModel->GetClusters(); i < length; i++)
				{
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fOpacity);
					QList<QList<ContourLine>> contours = _pModel->GetContour(i);
					for (int j = 0; j < contours.size(); j++)
					{
						drawContourLine(contours[j]);
					}
				}
			}
			/*
			QList<QList<ContourLine>>* contours = _pModel->GetContour();
			glColor4f(1.0, 0.0, 0.0, .1);
			for (int i = 0; i < contours[1].size(); i++)
			{
				drawContourLine(contours[1][i]);
			}
			glColor4f(0.0, 1.0, 0.0, .1);
			for (int i = 0; i < contours[2].size(); i++)
			{
				drawContourLine(contours[2][i]);
			}
			*/
		}
		else {
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
	}
	if (states._bShowGradientE)
	{
		for (size_t i = 0; i < g_lenEnsembles; i++)
		{
			glColor4f(1.0,1.0,0.0,0.06);
			glCallList(_gllistG + i * 3 + 1);
		}
	}

	// union
	if (states._bShowUnionE)
	{

		double fTransparency = .2;
		if (g_bClustering) {
			if (nClusterIndex>-1)
			{
				// show single cluster
				for (size_t j = 0; j < 3; j++)
				{
					glColor4f(g_arrColors[j][0], g_arrColors[j][1], g_arrColors[j][2], fTransparency);
					glCallList(_gllistC + nClusterIndex * 3 + j);
				}
			}
			else {
				for (size_t i = 0, length = _pModel->GetClusters(); i < length; i++)
				{
					glColor4f(g_arrColors[i][0], g_arrColors[i][1], g_arrColors[i][2], fTransparency);
					glCallList(_gllistC + i* 3 + 1);
				}

			}
		}
		else {
			// render each area
			if (g_bShowUncertaintyOnly) {

				glColor4f(1, 1, 1, .5);
				glCallList(_gllist + 1);
			}
			else
			for (size_t j = 0; j < 3; j++)
			{
				glColor4f(g_arrColors[j][0], g_arrColors[j][1], g_arrColors[j][2], fTransparency);
				glCallList(_gllist + j);
			}
		}
	}
	if (states._bShowBeliefEllipse)
	{
//		glColor3f(.5, .5, 0);
		glColor3f(1, 0, 0);
		glPointSize(2.0f);
		glBegin(GL_POINTS);
		const std::vector<DPoint3> pts = _pModel->GetPoints();
		for (size_t i = 0; i < pts.size(); i++)
		{

			glVertex3d(pts[i].x* _fScaleW, pts[i].y* _fScaleH, 0);
		}
		glEnd();
	}
	glPopMatrix();

	drawVarChart();

	// draw pca point
	drawPCAPoints();

	// draw cluster bars
	drawClusterBars();

	glPopAttrib();
}

void EnsembleLayer::drawVarChart() {
	// draw chart
	double dbLayer = .05;							// distance between two layers
	int nLayers = 10;								// layers of the chart

	double dbChartBottom = _fTop + .01;				// bottom position of the chart, on top of the map
	double dbChartLeft = 0.2;						// left position of the chart, 0.2
	double dbChartRight = _fRight;
	double dbChartTop = dbChartBottom + dbLayer*nLayers;

	// draw chart framework
	glLineWidth(2.0f);
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	// bottom axis
	glVertex2d(dbChartLeft, dbChartTop);
	glVertex2d(dbChartRight, dbChartTop);
	// top axis
	glVertex2d(dbChartLeft, dbChartBottom);
	glVertex2d(_fRight, dbChartBottom);
	// left axis
	glVertex2d(dbChartLeft, dbChartBottom);
	glVertex2d(dbChartLeft, dbChartTop);
	// right axis
	glVertex2d(dbChartRight, dbChartBottom);
	glVertex2d(dbChartRight, dbChartTop);
	glEnd();

	glLineWidth(1.0f);
	glColor3f(0, 0, .5);
	// scale
	for (size_t i = 1; i < nLayers; i++)
	{
		glBegin(GL_LINES);
		glVertex2d(dbChartLeft, dbChartBottom + i * dbLayer);
		glVertex2d(dbChartRight, dbChartBottom + i * dbLayer);
		glEnd();
		char buf[10];
		sprintf_s(buf, "%d", i);
		_pCB->DrawText(buf, dbChartLeft + .01, dbChartBottom + i * dbLayer + .01);
	}

	// draw the chart line
	glLineWidth(3.0f);
	glColor3f(0.0, .5, .5);
	std::vector<double> vecVar = _pModel->GetVariance();
	int nLen = vecVar.size();
	double dbStep = (_fRight - dbChartLeft) / (nLen - 1);
	double dbThresholdX = -1;
	glBegin(GL_LINE_STRIP);
	for (size_t i = 0; i < nLen; i += 100)
	{
		double dbYBias = vecVar[i] * dbLayer;
		if (dbYBias>dbLayer*nLayers)
		{
			dbYBias = dbLayer*nLayers;
		}
		glVertex2d(dbChartLeft + dbStep*i, dbChartBottom + dbYBias);
		if (dbThresholdX<0 && vecVar[i]>_dbVarThreshold) dbThresholdX = dbChartLeft + dbStep*i;
	}
	glEnd();

	// draw the threshold line
	if (dbThresholdX > 0) {
		glColor3f(1.0, 0, 0);
		glBegin(GL_LINES);
		glVertex2d(dbThresholdX, dbChartBottom);
		glVertex2d(dbThresholdX, dbChartBottom + nLayers*dbLayer);
		glEnd();

		char buf[10];
		sprintf_s(buf, "%.2f", _dbVarThreshold);
		_pCB->DrawText(buf, dbThresholdX, dbChartBottom + nLayers*dbLayer + .01);
	}
}

void EnsembleLayer::drawPCAPoints() {
	double dbRadius = .5;
	double dbSpace = .01;
	double dbPCAChartLeft = _fLeft + dbRadius;
	double dbPCAChartBottom = _fTop + dbRadius;
	double dbScale = .01;
	glPointSize(3.0f);

//	glPushMatrix();
//	glTranslated(_fLeft / 2.0, _fTop + .5, 0);

	std::vector<DPoint3>* vecPoints = _pModel->GetPCAPoints();
	for (size_t clusterIndex = 0; clusterIndex < 2; clusterIndex++)
	{
		// draw border
		// set color according to the color of the spatial cluster
		switch (clusterIndex)
		{
		case 0:
			glColor3f(1, 0, 0);
			break;			
		case 1:
			glColor3f(0, 1, 0);
			break;
		default:
			break;
		}
		glBegin(GL_LINE_LOOP);
		glVertex3f(dbPCAChartLeft - (dbRadius - dbSpace), dbPCAChartBottom - (dbRadius - dbSpace), 0);
		glVertex3f(dbPCAChartLeft + (dbRadius - dbSpace), dbPCAChartBottom - (dbRadius - dbSpace), 0);
		glVertex3f(dbPCAChartLeft + (dbRadius - dbSpace), dbPCAChartBottom + (dbRadius - dbSpace), 0);
		glVertex3f(dbPCAChartLeft - (dbRadius - dbSpace), dbPCAChartBottom + (dbRadius - dbSpace), 0);
		glEnd();
//		glPushMatrix();
//		glScaled(.01, .01, 1);
//		qDebug() << vecPoints[clusterIndex].size();
		glBegin(GL_POINTS);
		for (size_t i = 0, length = vecPoints[clusterIndex].size(); i < length; i++)
		{
			int nZ = vecPoints[clusterIndex][i].z;
			switch (nZ)
			{
			case 0:
				glColor3f(1.0, 0, 0);		// red
				break;
			case 1:
				glColor3f(0, 1.0, 0);		// green
				break;
			case 2:
				glColor3f(0, 0, 1.0);		// blue
				break;
			case 3:
				glColor3f(0, 1, 1);			// Cyan
				break;
			case 4:
				glColor3f(1, 0, 1);			// purple
				break;
			default:
				glColor3f(.5, .5, .5);
				break;
			}
			glVertex3f(dbPCAChartLeft + vecPoints[clusterIndex][i].x*dbScale, dbPCAChartBottom + vecPoints[clusterIndex][i].y*dbScale, 0);
		}
		glEnd();
		dbPCAChartLeft += 2*dbRadius;
//		glPopMatrix();

//		glTranslated(1,0,0);
	}

//	glPopMatrix();
}


void EnsembleLayer::drawClusterBars() {
	const ClusterResult* pCR = _pModel->GetClusterResults();
	if (pCR->_nK == 0) return;

	double dbRadius = .5;
	double dbSpace = .01;
	double dbSpaceII = .05;
	double dbChartLeft = _fLeft + dbSpace;
	double dbChartRight = dbChartLeft + dbRadius * 2 - dbSpace;
	double dbChartMid = _fLeft + dbRadius;
	double dbChartBottom = _fTop + dbRadius*2+ dbSpaceII;		// on top of the pca chart

	int nClusters = _pModel->GetClusters();

	double arrBaseY[50];			// the y position of the element in the first cluster

	// first region
	double dbY = dbChartBottom;
	for (size_t i = 0; i < pCR[0]._nK; i++)
	{
		SetGroupColor(i);

		for (size_t j = 0; j < pCR[0]._vecItems[i].size(); j++) {
			glBegin(GL_LINES);
			glVertex3f(dbChartLeft, dbY,0);
			glVertex3f(dbChartRight, dbY, 0);
			glEnd();
			arrBaseY[pCR[0]._vecItems[i][j]] = dbY;
			dbY += dbSpace;
		}
		dbY += dbSpaceII;
	}
	// second region
	dbChartLeft += dbRadius * 2;
	dbChartMid += dbRadius * 2;
	dbChartRight += dbRadius * 2;
	/*
	// left part
	dbY = dbChartBottom;
	for (size_t i = 0; i < pCR[0]._nK; i++)
	{
		for (size_t j = 0; j < pCR[0]._vecItems[i].size(); j++) {

			SetGroupColor(pCR[1]._arrLabels[pCR[0]._vecItems[i][j]]);

			glBegin(GL_LINES);
			glVertex3f(dbChartLeft, dbY, 0);
			glVertex3f(dbChartMid, dbY, 0);
			glEnd();
			dbY += dbSpace;
		}
		dbY += dbSpaceII;
	}
	*/
	// right part
	dbY = dbChartBottom;
	for (size_t i = 0; i < pCR[1]._nK; i++)
	{
		for (size_t j = 0; j < pCR[1]._vecItems[i].size(); j++) {
			SetGroupColor(i);

			glBegin(GL_LINES);
			glVertex3f(dbChartLeft, arrBaseY[pCR[1]._vecItems[i][j]], 0);
			glVertex3f(dbChartMid, dbY, 0);
			glEnd();

			glBegin(GL_LINES);
			glVertex3f(dbChartMid, dbY, 0);
			glVertex3f(dbChartRight, dbY, 0);
			glEnd();
			dbY += dbSpace;
		}
		dbY += dbSpaceII;
	}
}

void EnsembleLayer::ReloadTexture() {
	// Enable texturing
	// 	generateBackground();
	// 	glEnable(GL_TEXTURE_2D);
//	_dataTexture = _pModel->generateTexture();
//	_dataTexture = _pModel->generateTextureGridCluster();
//	_dataTexture = _pModel->generateTextureRange(0);
//	_dataTexture = _pModel->generateTextureMean();
//	_dataTexture = _pModel->generateTextureDiscreteSummary();
	_dataTexture = _pModel->generateTextureNew();
	//	_dataTexture = _pModel->generateTextureSDF();
	glGenTextures(1, &texID[0]);
	glBindTexture(GL_TEXTURE_2D, texID[0]);
	// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _pModel->GetFocusW(), _pModel->GetFocusH(), 0, GL_RGBA, GL_UNSIGNED_BYTE, _dataTexture);

}

void EnsembleLayer::init(){
	ReloadTexture();
	// color bar
	ColorMap* colormap=ColorMap::GetInstance();
	int nLen = colormap->GetLength();
	int nStep = colormap->GetStep();
	int nColorBarW = 1;
	int nColorBarH = (nLen-1)*10+1;
	_colorbarTexture = new GLubyte[nColorBarH * 4];
	for (int i = 0; i < nColorBarH; i++)
	{
		MYGLColor color = colormap->GetColor(i/10.0*nStep);
		// using transparency and the blue tunnel
		_colorbarTexture[4 * i + 0] = color._rgb[0];
		_colorbarTexture[4 * i + 1] = color._rgb[1];
		_colorbarTexture[4 * i + 2] = color._rgb[2];
		_colorbarTexture[4 * i + 3] = (GLubyte)255;
	}
	glGenTextures(1, &texID[1]);
	glBindTexture(GL_TEXTURE_2D, texID[1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nColorBarW, nColorBarH, 0, GL_RGBA, GL_UNSIGNED_BYTE, _colorbarTexture);
	// ~build texture


	// initialize tess
	int nClusterLen = _pModel->GetClusters();
	_gllist = glGenLists(3);					// generate the display lists
	_gllistG = glGenLists(g_gradient_l * 3);	// generate the display lists
	_gllistC = glGenLists(nClusterLen * 3);		// generate the display lists


	_tobj = gluNewTess();
	gluTessCallback(_tobj, GLU_TESS_VERTEX,
		(void(__stdcall*)())glVertex3dv);
	gluTessCallback(_tobj, GLU_TESS_BEGIN,
		(void(__stdcall*)())beginCallback);
	gluTessCallback(_tobj, GLU_TESS_END,
		(void(__stdcall*)())endCallback);
	gluTessCallback(_tobj, GLU_TESS_ERROR,
		(void(__stdcall*)())errorCallback);



	tessSegmentation(_gllist, _pModel->GetUncertaintyArea());
	if (g_bClustering)
	{
		for (size_t i = 0; i < nClusterLen; i++)
		{
			tessSegmentation(_gllistC + i * 3, _pModel->GetUncertaintyArea(i));
		}
	}
	QList<QList<UnCertaintyArea*> > listAreas = _pModel->GetUncertaintyAreaG();
	for (size_t i = 0, len = listAreas.length(); i < len; i++)
	{
		tessSegmentation(_gllistG + i * 3, listAreas[i]);
	}
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
			contourBuffer[i][0][j][0] = areas[i]->_contour._listPt[j].x()*_fScaleW;
			contourBuffer[i][0][j][1] = areas[i]->_contour._listPt[j].y()*_fScaleH;
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
				contourBuffer[i][k + 1][j][0] = areas[i]->_listEmbeddedArea[k]->_contour._listPt[j].x()*_fScaleW;
				contourBuffer[i][k + 1][j][1] = areas[i]->_listEmbeddedArea[k]->_contour._listPt[j].y()*_fScaleH;
				contourBuffer[i][k + 1][j][2] = 0;
			}
		}
	}
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

void EnsembleLayer::OnSelectVar(int nIndex) {
	std::vector<double> vecVar = _pModel->GetVariance();
	_dbVarThreshold = vecVar[nIndex*vecVar.size() / 100];
	_pModel->SetVarThreshold(_dbVarThreshold);
	ReloadTexture();
}