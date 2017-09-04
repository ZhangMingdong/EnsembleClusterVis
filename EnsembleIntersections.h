#pragma once
#include <QList>
// an intersection
struct Intersection
{
	double _fMin;		// the list of intersection min
	double _fMax;		// the list of intersection max
	int _nOverlap;		// the overlaps for each intersection
public:
	Intersection();
	Intersection(double fMin, double fMax, int nOverlap);
};
// store the intersections
struct EnsembleIntersections
{
	QList<Intersection> _listIntersection;
public:
	// intersect with a new segment
	void Intersect(double fMin, double fMax);
	// get the most overlapped intersection
	void GetMostOverlapped(double& fMin, double& fMax);
	// get the most overlapped intersection by threshold of intersection ensembles
	void GetMostOverlapped(double& fMin, double& fMax, int nThreshold);
};