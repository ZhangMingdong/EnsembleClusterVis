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


// length of the ensembles
const int g_lenEnsembles = 50;

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



const int g_temperatureLen = 1;


// new switch
const bool g_bClustering = false;			// whether do clustering (PCA+AHC)
const bool g_bSpatialClustering = false;	// whether do spatial clustering (DBSCAN)

const double g_fThreshold = 1;// 273.16 - 15;
//const double g_fThreshold = 2.0;
const int g_nMinPts = 110;
const double g_dbEps = 6;

const bool g_bGlobalArea = true;

// whether filter the half degree data
const bool g_bFilter = false;
//	bFilter = false;
const bool g_bSubArea = false;

// whether calculate matrix for new data or use stored one
const bool g_bNewData = false;

const int g_nClusters = 6;

// time steps in the data file
const int g_nTimeStep = 0;



const int g_globalW = 360;
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
	, T2_Reanalysis
};

// used model
//const enumMeteModel g_usedModel = PRE_CMA;
//const enumMeteModel g_usedModel = PRE_CPTEC;
//const enumMeteModel g_usedModel = PRE_ECCC;
const enumMeteModel g_usedModel = PRE_ECMWF;
//const enumMeteModel g_usedModel = PRE_JMA;
//const enumMeteModel g_usedModel = PRE_KMA;
//const enumMeteModel g_usedModel = PRE_NCEP;
//const enumMeteModel g_usedModel = T2_ECMWF;
//const enumMeteModel g_usedModel = T2_Reanalysis;

// just use white to show the uncertainty area
const bool g_bShowUncertaintyOnly = false;

// calculate the uncertainty band based on the signed distance function, otherwise calculate directly
const bool g_bSDFBand = false;

// the index of the cluster to generate confidence Ellipse, -1 means all
const int g_nConfidenceEllipseIndex = -1;

// the Mahalanobis distance used in calculating the confidence ellipse
const double g_dbMDis = 2.0;

// whether use ensemble model
const bool g_bEnsembleModel = false;

