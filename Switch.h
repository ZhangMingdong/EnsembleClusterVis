#pragma once
/*
	switch definitions
*/

//====================feature================
const bool g_bStatistic = true;

const bool g_bGenerateContours = true;
const bool g_bGenerateArea = true; // my algorithm cannot handle the case there's closed areas only

const bool g_bSmoothContours = false;

const bool g_bCalculateSDF = true;

const bool g_bResampleContours = false;

const bool g_bDomainResampleContours = false;

const bool g_bCalculateICD = false;

const bool g_bClustering = true;
const bool g_bResampleForClusters = true;		// resample for each clusters



// for test
const bool g_bTestInfoLoseMeasure = false;
const bool g_bMutualInfoLose = false;			// use mutual information methods
const int g_nTestTarget = 0;					// which to test: 0-reconplot;1-raw_random;2-sorted_central;3:sorted_sequential


// ==================Model=========================
const bool g_bEOF = false;

const bool g_bArtificialModel = false;

// ==================Frame=========================

const bool g_bChart = false;	// if show the chart

// ==================Chart=========================
const bool g_bDrawHierarchy = false;
const bool g_bDrawCluster = true;

// ==================Debug==========================
const bool g_bDebugDataField = false;