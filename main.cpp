#include <QApplication>

#include "MainWindow.h"


#include <iostream>
#include "INIReader.h"

double g_fThreshold = 273.16 - 15;
bool g_bGlobalArea = false;
enumMeteModel g_usedModel = T2_ECMWF;

int g_nWidth = 0;
int g_nHeight = 0;
int g_nFocusX = 0;
int g_nFocusY = 0;
int g_nFocusW = 0;
int g_nFocusH = 0;
int g_nWest = 0;
int g_nEast = 0;
int g_nNorth = 0;
int g_nSouth = 0;
char g_strFileName[50];

char g_strObsFileName[50];
int g_nBiasY=0;
int g_nBiasD=0;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-72.ini";
//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-96.ini";
//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-120.ini";
//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-144.ini";
//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-168.ini";
//	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-192.ini";
	const char* strConfigFileName = "../../data/config/t2-mod-ecmwf-20160105-00-216.ini";
//	const char* strConfigFileName = "../../data/config/pre_20171016_0-360-by-6.ini";
	// read config file
	INIReader reader(strConfigFileName);

	if (reader.ParseError() < 0) {
		std::cout << "Can't load '"<< strConfigFileName <<"'\n";
		return 1;
	}
	g_bGlobalArea = reader.GetBoolean("region", "global",true);
	g_usedModel = (enumMeteModel)reader.GetInteger("model", "model", -1);
	g_fThreshold = reader.GetReal("value", "threshold", -1);



	g_nWidth = reader.GetInteger("region", "width", 0);
	g_nHeight = reader.GetInteger("region", "height", 0);
	g_nFocusX = reader.GetInteger("region", "focusx", 0);
	g_nFocusY = reader.GetInteger("region", "focusy", 0);
	g_nFocusW = reader.GetInteger("region", "focusw", 0);
	g_nFocusH = reader.GetInteger("region", "focush", 0);
	g_nWest = reader.GetInteger("region", "west", 0);
	g_nEast = reader.GetInteger("region", "east", 0);
	g_nNorth = reader.GetInteger("region", "north", 0);
	g_nSouth = reader.GetInteger("region", "south", 0);
	strcpy(g_strFileName, reader.Get("file", "name", "").c_str());


	strcpy(g_strObsFileName, reader.Get("obs", "file", "").c_str());
	g_nBiasY = reader.GetInteger("obs", "biasy", 0);
	g_nBiasD = reader.GetInteger("obs", "biasd", 0);

	app.setApplicationName(app.translate("main", "EnsembleVis"));
	app.setOrganizationName("CG&CAD");
	app.setOrganizationDomain("cg&cad.tsinghua.edu.cn");

	g_fThreshold= 273.16+15;

	MainWindow mainwindow;
	mainwindow.show();
    return app.exec();
}
