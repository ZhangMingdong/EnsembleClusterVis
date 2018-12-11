// ============================================================
// PCA
// Written by Mingdong
// 2017/05/13
// ============================================================
#pragma once

#include "tnt_array2d.h"
#include <MathTypes.hpp>
#include <vector>
class MyPCA
{
public:
	MyPCA();
	~MyPCA();
public:


	/*
		do the pca
		n:			number of vectors
		mI:			length of the input vectors
		mO:			length of the output vectors
		arrInput:	input data: n*m
		arrOutput:	output data: n*(n-1)
		bNewData:	whether recalculate for new data
	*/
	void DoPCA(const double* arrInput, double* arrOutput, int n, int mI, int mO, bool bNewData = false);

	/*
		recover nRow of vectors
		nRow: how many data item
		nCol: length of pca
	*/
	void Recover(double* arrInput, double* arrOutput,int nRow,int nCol);
	/*
		recover 1 data item
		nCol: length of pca
	*/
	void Recover(double* arrInput, double* arrOutput, int nCol);

	// test pca and recover
	void TestPCARecover(double* arrInput, int n, int mI, int mO, bool bNewData = false);
private:
	int _nRow;
	int _nCol;
	double* _arrDataMean;					// array of data means
	TNT::Array2D<double>* _pR;				// function to transform back
	TNT::Array2D<double>* _pCovar_matrix;	// covariance matrix
	TNT::Array2D<double>* _pEigenvector;	// eigen vector
	TNT::Array2D<double>* _pEigenvalue;		// eigen value
	TNT::Array2D<double>* _pD;				// input data
private:
	// read the four matrices
	void readMatrices();
	// calculate the four matrices
	void calculateMatrices();
	// write the four matrices
	void writeMatrices();
public:
	/*
		points:			return the generated points
		arrGridLabel:	label of the grid points
		nWidth,nHeight:	the width and height of the field
		dbDis:			the threshold
	*/
	static void generateEllipse(std::vector<DPoint3>& points, int* arrGridLabel, int label, int nWidth, int nHeight,double dbDis);
};

