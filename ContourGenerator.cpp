#include "ContourGenerator.h"
#include "def.h"
#include <QDebug>




LineSeg::LineSeg(int i11, int j11, int i12, int j12, int i21, int j21, int i22, int j22){
	_arrIndices[0] = i11;
	_arrIndices[1] = j11;
	_arrIndices[2] = i12;
	_arrIndices[3] = j12;
	_arrIndices[4] = i21;
	_arrIndices[5] = j21;
	_arrIndices[6] = i22;
	_arrIndices[7] = j22;
}


ContourGenerator::ContourGenerator()
{
}


ContourGenerator::~ContourGenerator()
{
}

void ContourGenerator::Generate(const double* data
	, QList<ContourLine>& listContour
	, double fValue
	, int nWidth
	, int nHeight
	, int nFocusX
	, int nFocusY
	, int nFocusW
	, int nFocusH
	){
	// 0.set value
	_data = data;
	_fValue = fValue+0.0001;
	_nWidth = nWidth;
	_nHeight = nHeight;
	_nFocusX = nFocusX;
	_nFocusY = nFocusY;
	_nFocusW = nFocusW;
	_nFocusH = nFocusH;
	// 1.generate line strips
// 	generateLineStrip_old();
	generateLineStrip();
	// 2 calculate segment in domain
	combineLineStrip(listContour);
}


void ContourGenerator::generateLineStrip(){
	// 1.generate segment
// 	for (int i = 1; i < _nHeight; i++){
// 		for (int j = 1; j < _nWidth; j++){
	for (int i = _nFocusY + 1, iLen = _nFocusY + _nFocusH; i < iLen; i++){
		for (int j = _nFocusX + 1, jLen = _nFocusX + _nFocusW; j < jLen; j++){
			int grid[4][2] = {
				{ i - 1, j - 1 }
				, { i - 1, j }
				, { i, j }
				, { i, j - 1 }
			};
			double v[4] = {
				_data[((i - 1)) * _nWidth + j - 1],
				_data[((i - 1)) * _nWidth + j],
				_data[(i) * _nWidth + j],
				_data[(i) * _nWidth + j - 1],
			};
			int t = 0;
			for (int k = 0; k<4; k++) if (v[k]>_fValue) t++;
			if (t == 1){
				int corner = 0;
				while (v[corner] < _fValue) corner++;
				int c1 = (corner + 1) % 4;
				int c2 = (corner + 3) % 4;
				_listSeg.append(LineSeg(grid[corner][0], grid[corner][1]
					, grid[c1][0], grid[c1][1]
					, grid[corner][0], grid[corner][1]
					, grid[c2][0], grid[c2][1]
					));
			}
			else if (t == 2)
			{
				bool state[4] = {
					v[0] > _fValue
					, v[1] > _fValue
					, v[2] > _fValue
					, v[3] > _fValue
				};
				if (state[1] != state[2] && state[2] != state[3]){
					// 1,3 : 2,4
					_listSeg.append(LineSeg(grid[0][0], grid[0][1]
						, grid[1][0], grid[1][1]
						, grid[1][0], grid[1][1]
						, grid[2][0], grid[2][1]
						));
					_listSeg.append(LineSeg(grid[2][0], grid[2][1]
						, grid[3][0], grid[3][1]
						, grid[3][0], grid[3][1]
						, grid[0][0], grid[0][1]
						));
				}
				else{
					int c0 = 0;
					while (state[c0] != state[c0 + 1]) c0++;
					int c1 = (c0 + 1) % 4;
					int c2 = (c0 + 2) % 4;
					int c3 = (c0 + 3) % 4;
					_listSeg.append(LineSeg(grid[c0][0], grid[c0][1]
						, grid[c3][0], grid[c3][1]
						, grid[c2][0], grid[c2][1]
						, grid[c1][0], grid[c1][1]
						));

				}
			}
			else if (t==3)
			{
				int corner = 0;
				while (v[corner] > _fValue) corner++;
				int c1 = (corner + 1) % 4;
				int c2 = (corner + 3) % 4;
				_listSeg.append(LineSeg(grid[corner][0], grid[corner][1]
					, grid[c1][0], grid[c1][1]
					, grid[corner][0], grid[corner][1]
					, grid[c2][0], grid[c2][1]
					));
			}

		}
	}
}

void ContourGenerator::generateLineStrip_old(){
	// 1.generate segment
	for (int i = 1; i < _nHeight; i++)
	{
		for (int j = 1; j < _nWidth; j++)
		{
			int grid[2][3][2] = {
				{
					{ i - 1, j - 1 }
					, { i - 1, j }
					, { i, j }

				}
				, {
					{ i - 1, j - 1 }
					, { i, j - 1 }
					, { i, j }

				}
			};
			double v[2][3] = {
				{
					_data[(_nHeight - i) * _nWidth + j - 1],
					_data[(_nHeight - i) * _nWidth + j],
					_data[(_nHeight - i - 1) * _nWidth + j],

				}, {
					_data[(_nHeight - i) * _nWidth + j - 1],
					_data[(_nHeight - i - 1) * _nWidth + j - 1],
					_data[(_nHeight - i - 1) * _nWidth + j],

				}
			};
			// treat two triangles
			for (int tri = 0; tri < 2; tri++)
			{
				int t = 0;
				for (int k = 0; k<3; k++) if (v[tri][k]>_fValue) t++;
				int c1, c2, corner;
				if (t == 1){
					int corner = 0;
					while (v[tri][corner] < _fValue) corner++;
					int c1 = (corner + 1) % 3;
					int c2 = (corner + 2) % 3;
					_listSeg.append(LineSeg(grid[tri][corner][0], grid[tri][corner][1]
						, grid[tri][c1][0], grid[tri][c1][1]
						, grid[tri][corner][0], grid[tri][corner][1]
						, grid[tri][c2][0], grid[tri][c2][1]
						));
				}
				else if (t == 2)
				{
					int corner = 0;
					while (v[tri][corner] > _fValue) corner++;
					int c1 = (corner + 1) % 3;
					int c2 = (corner + 2) % 3;
					_listSeg.append(LineSeg(grid[tri][corner][0], grid[tri][corner][1]
						, grid[tri][c1][0], grid[tri][c1][1]
						, grid[tri][corner][0], grid[tri][corner][1]
						, grid[tri][c2][0], grid[tri][c2][1]
						));
				}
			}

		}
	}
}

void ContourGenerator::combineLineStrip(QList<ContourLine>& listContour){
	// process until the strip list is empty
	while (!_listSeg.empty())
	{
		// 1.get the first segment as the start point
		LineSeg lineSeg = _listSeg[0];
		_listSeg.removeFirst();
		ContourLine contour;

		// 2.add the first two points
		double v[4] = {
			_data[(lineSeg._arrIndices[0]) * _nWidth + lineSeg._arrIndices[1]],
			_data[(lineSeg._arrIndices[2]) * _nWidth + lineSeg._arrIndices[3]],
			_data[(lineSeg._arrIndices[4]) * _nWidth + lineSeg._arrIndices[5]],
			_data[(lineSeg._arrIndices[6]) * _nWidth + lineSeg._arrIndices[7]]
		};

		double bias0 = abs(v[0] - _fValue);
		double bias1 = abs(v[1] - _fValue);
		double bias2 = abs(v[2] - _fValue);
		double bias3 = abs(v[3] - _fValue);

		double x0 = bias1 / (bias0 + bias1)*lineSeg._arrIndices[0] + bias0 / (bias0 + bias1)*lineSeg._arrIndices[2];
		double y0 = bias1 / (bias0 + bias1)*lineSeg._arrIndices[1] + bias0 / (bias0 + bias1)*lineSeg._arrIndices[3];
		double x1 = bias3 / (bias2 + bias3)*lineSeg._arrIndices[4] + bias2 / (bias2 + bias3)*lineSeg._arrIndices[6];
		double y1 = bias3 / (bias2 + bias3)*lineSeg._arrIndices[5] + bias2 / (bias2 + bias3)*lineSeg._arrIndices[7];
		contour._listPt.append(QPointF(y0, x0));
		contour._listPt.append(QPointF(y1, x1));

		// 3.record the position of the first two points
		int nLastIndices[8] = {
			lineSeg._arrIndices[0],
			lineSeg._arrIndices[1],
			lineSeg._arrIndices[2],
			lineSeg._arrIndices[3],
			lineSeg._arrIndices[4],
			lineSeg._arrIndices[5],
			lineSeg._arrIndices[6],
			lineSeg._arrIndices[7]
		};
		// 4.go on to finished the current contour
		bool bContinue = true; // this contour is continue
		while (bContinue&&!_listSeg.empty())
		{
			bContinue = false;
			for (int i = 0, len = _listSeg.length(); i < len; i++)
			{
				LineSeg lineSeg = _listSeg[i];
				bool bTail = true;				// add the new point to tail
				// check whether and how to connect
				if ((lineSeg._arrIndices[0] == nLastIndices[0]
					&& lineSeg._arrIndices[1] == nLastIndices[1]
					&& lineSeg._arrIndices[2] == nLastIndices[2]
					&& lineSeg._arrIndices[3] == nLastIndices[3])
					||
					(lineSeg._arrIndices[0] == nLastIndices[2]
					&& lineSeg._arrIndices[1] == nLastIndices[3]
					&& lineSeg._arrIndices[2] == nLastIndices[0]
					&& lineSeg._arrIndices[3] == nLastIndices[1])
					)
				{
					nLastIndices[0] = lineSeg._arrIndices[4];
					nLastIndices[1] = lineSeg._arrIndices[5];
					nLastIndices[2] = lineSeg._arrIndices[6];
					nLastIndices[3] = lineSeg._arrIndices[7];
					bTail = false;
					bContinue = true;
				}
				else if ((
					lineSeg._arrIndices[4] == nLastIndices[0]
					&& lineSeg._arrIndices[5] == nLastIndices[1]
					&& lineSeg._arrIndices[6] == nLastIndices[2]
					&& lineSeg._arrIndices[7] == nLastIndices[3])
					||
					(
					lineSeg._arrIndices[4] == nLastIndices[2]
					&& lineSeg._arrIndices[5] == nLastIndices[3]
					&& lineSeg._arrIndices[6] == nLastIndices[0]
					&& lineSeg._arrIndices[7] == nLastIndices[1])

					)
				{
					nLastIndices[0] = lineSeg._arrIndices[0];
					nLastIndices[1] = lineSeg._arrIndices[1];
					nLastIndices[2] = lineSeg._arrIndices[2];
					nLastIndices[3] = lineSeg._arrIndices[3];
					bTail = false;
					bContinue = true;
				}
				else if ((
					lineSeg._arrIndices[0] == nLastIndices[4]
					&& lineSeg._arrIndices[1] == nLastIndices[5]
					&& lineSeg._arrIndices[2] == nLastIndices[6]
					&& lineSeg._arrIndices[3] == nLastIndices[7])
					||
					(
					lineSeg._arrIndices[0] == nLastIndices[6]
					&& lineSeg._arrIndices[1] == nLastIndices[7]
					&& lineSeg._arrIndices[2] == nLastIndices[4]
					&& lineSeg._arrIndices[3] == nLastIndices[5])

					)
				{
					nLastIndices[4] = lineSeg._arrIndices[4];
					nLastIndices[5] = lineSeg._arrIndices[5];
					nLastIndices[6] = lineSeg._arrIndices[6];
					nLastIndices[7] = lineSeg._arrIndices[7];
					bContinue = true;
				}
				else if ((
					lineSeg._arrIndices[4] == nLastIndices[4]
					&& lineSeg._arrIndices[5] == nLastIndices[5]
					&& lineSeg._arrIndices[6] == nLastIndices[6]
					&& lineSeg._arrIndices[7] == nLastIndices[7])
					||
					(
					lineSeg._arrIndices[4] == nLastIndices[6]
					&& lineSeg._arrIndices[5] == nLastIndices[7]
					&& lineSeg._arrIndices[6] == nLastIndices[4]
					&& lineSeg._arrIndices[7] == nLastIndices[5])

					)
				{
					nLastIndices[4] = lineSeg._arrIndices[0];
					nLastIndices[5] = lineSeg._arrIndices[1];
					nLastIndices[6] = lineSeg._arrIndices[2];
					nLastIndices[7] = lineSeg._arrIndices[3];
					bContinue = true;
				}


				if (bContinue){
					// add the first two points
					int nIndices[4];
					if (bTail)
					{
						nIndices[0] = nLastIndices[4];
						nIndices[1] = nLastIndices[5];
						nIndices[2] = nLastIndices[6];
						nIndices[3] = nLastIndices[7];
					}
					else
					{
						nIndices[0] = nLastIndices[0];
						nIndices[1] = nLastIndices[1];
						nIndices[2] = nLastIndices[2];
						nIndices[3] = nLastIndices[3];
					}
					
					double v[2] = {
						_data[(nIndices[0]) * _nWidth + nIndices[1]],
						_data[(nIndices[2]) * _nWidth + nIndices[3]]
					};

					double bias0 = abs(v[0] - _fValue);
					double bias1 = abs(v[1] - _fValue);

					double x0 = bias1 / (bias0 + bias1)*nIndices[0] + bias0 / (bias0 + bias1)*nIndices[2];
					double y0 = bias1 / (bias0 + bias1)*nIndices[1] + bias0 / (bias0 + bias1)*nIndices[3];
					if (bTail)
					{
						contour._listPt.append(QPointF(y0, x0));
					}
					else{
						contour._listPt.insert(0,QPointF(y0, x0));
					}
					
					_listSeg.removeAt(i);
					break;
				}
			}
		}
		// check whether contour is closed
		{
			QPointF pt1 = contour._listPt[0];
			QPointF pt2 = contour._listPt[contour._listPt.size() - 1];
			if (abs(pt1.x() - pt2.x()) < g_fDifference && (abs(pt1.y() - pt2.y()) < g_fDifference)) contour._nClosed = 1;
			else contour._nClosed = 0;			
		}
		// ignore the contour shorter than 4
		if (contour._listPt.size()>4)
		{
			listContour.append(contour);
		}
		

	}
}