#pragma once
//20111214
struct DPoint3 {
	DPoint3() { x = y = z = 0.0f; }
	DPoint3(double fx, double fy, double fz) { x = fx; y = fy; z = fz; }
	DPoint3 &operator=(const DPoint3 &v) { x = v.x; y = v.y; z = v.z; return *this; }
	double x, y, z;
};//--------------------------------------------------------------------------------
class Point : public DPoint3
{
public:
	Point() {};
	Point(DPoint3& p) :DPoint3(p), _ulPos(0) {}
	int _ulPos;
};
