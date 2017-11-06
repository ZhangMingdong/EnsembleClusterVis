#pragma once

const double g_fDifference = 0.0001;

/*
*/
#define GLOBAL_SCOPE

// width and height of the grid
// int const g_grid_width = 181;
// // int const g_grid_height = 101;
int const g_grid_width = 360;
int const g_grid_height = 181;
// int const g_grid_width = 181;
// int const g_grid_height = 360;
int const g_dataLen = g_grid_width*g_grid_height;


//double const g_temperature = -10;


int const g_focus_x = 0;
int const g_focus_y = 0;
int const g_focus_w_raw = g_grid_width;
int const g_focus_h_raw = g_grid_height;
// width and height of the view
int const g_width = 1810;
int const g_height = 1010;



// threshold of intersection numbers
const int g_nThreshold = 40;

// space while using the data
const int g_nSpace = 1;
const int g_focus_w = g_focus_w_raw;
const int g_focus_h = g_focus_h_raw;


// the data length
int const g_focus_l = g_focus_w*g_focus_h;

// the gradient length
const int g_gradient_l = 15;



const bool g_bSubArea = false;

// whether calculate matrix for new data or use stored one
const bool g_bNewData = false;

const int g_nClusters = 6;




const int g_globalW = 361;
const int g_globalH = 181;

enum enumMeteModel
{
	  PRE_CMA
	, PRE_CPTEC
	, PRE_ECCC
	, PRE_ECMWF
	, PRE_JMA
	, PRE_KMA
	, PRE_NCEP
	, T2_ECMWF
	, PRE_ECMWF_2017
};



// just use white to show the uncertainty area
const bool g_bShowUncertaintyOnly = false;

// calculate the uncertainty band based on the signed distance function, otherwise calculate directly
const bool g_bSDFBand = false;

// the index of the cluster to generate confidence Ellipse, -1 means all
const int g_nConfidenceEllipseIndex = -1;

// the Mahalanobis distance used in calculating the confidence ellipse
const double g_dbMDis = 2.0;



// 2017/09/14
const int g_nEnsembles = 50;					// number of ensemble members
const int g_nUncertaintyAreaMax = 6;			// max number of uncertainty area
const int g_nClusterMax = 10;					// max number of clusters

const int g_nIsoValuesLen = 5;					// length of array of isovalues
// 2017/10/17
//# define GLOBAL_PRE

#ifdef GLOBAL_PRE			
const double g_fThreshold = 2.0;
const bool g_bGlobalArea = true;
const int g_arrIsoValues[5] = { 1,2,3,4,5 };

// used model
//const enumMeteModel g_usedModel = PRE_CMA;
//const enumMeteModel g_usedModel = PRE_CPTEC;
//const enumMeteModel g_usedModel = PRE_ECCC;
//const enumMeteModel g_usedModel = PRE_ECMWF;
//const enumMeteModel g_usedModel = PRE_JMA;
//const enumMeteModel g_usedModel = PRE_KMA;
//const enumMeteModel g_usedModel = PRE_NCEP;
const enumMeteModel g_usedModel = PRE_ECMWF_2017;

#else
const double g_fThreshold = 273.16 - 15;
const bool g_bGlobalArea = false;
const enumMeteModel g_usedModel = T2_ECMWF;
const int g_arrIsoValues[5] = { 273.16-20,273.16-10,273.16,273.16+10,273.16+20 };

#endif

// 2017/10/23
const int g_nEOFLen = 10;