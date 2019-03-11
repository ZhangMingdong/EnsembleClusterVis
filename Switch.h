#pragma once
/*
	switch definitions
*/

//====================feature================
const bool g_bCalculateSet = true;
const bool g_bCalculateBandDepth = true;
const bool g_bCalculateMemberType = true;
const bool g_bStatistic = true;
const bool g_bGenerateContours = true;
const bool g_bSmoothContours = false;
const bool g_bCalculateSDF = true;
const bool g_bBuildSortedSDF = true;
const bool g_bResampleContours = true;
const bool g_bCalculateICD = false;
const bool g_bCalculateICD_L = false;
const bool g_bCalculateICD_V = false;
const bool g_bCalculateICDV = false;
const bool g_bCalculateDiverse = false;
const bool g_bCalculatePCA = false;
const bool g_bClustering = false;
const bool g_bGenerateArea = false; // my algorithm cannot handle the case there's closed areas only

// ==================Model=========================
const bool g_bEOF = false;

const bool g_bArtificialModel = true;

// ==================Frame=========================

const bool g_bChart = false;	// if show the chart