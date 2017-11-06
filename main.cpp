#include <QApplication>

#include "MainWindow.h"


#include <iostream>
#include "INIReader.h"

double g_fThreshold = 273.16 - 15;
bool g_bGlobalArea = false;
enumMeteModel g_usedModel = T2_ECMWF;
double g_arrIsoValues[5] = { 273.16 - 20,273.16 - 10,273.16,273.16 + 10,273.16 + 20 };

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

	const char* strConfigFileName = "../../data/config/1.ini";
	//const char* strConfigFileName = "../../data/config/2.ini";
	// read config file
	INIReader reader(strConfigFileName);

	if (reader.ParseError() < 0) {
		std::cout << "Can't load '"<< strConfigFileName <<"'\n";
		return 1;
	}
	g_bGlobalArea = reader.GetBoolean("region", "global",true);
	g_usedModel = (enumMeteModel)reader.GetInteger("model", "model", -1);
	g_fThreshold = reader.GetReal("value", "threshold", -1);
	g_arrIsoValues[0] = reader.GetReal("value", "i0", 0);
	g_arrIsoValues[1] = reader.GetReal("value", "i1", 0);
	g_arrIsoValues[2] = reader.GetReal("value", "i2", 0);
	g_arrIsoValues[3] = reader.GetReal("value", "i3", 0);
	g_arrIsoValues[4] = reader.GetReal("value", "i4", 0);

	app.setApplicationName(app.translate("main", "EnsembleVis"));
	app.setOrganizationName("CG&CAD");
	app.setOrganizationDomain("cg&cad.tsinghua.edu.cn");

	MainWindow mainwindow;
	mainwindow.show();
    return app.exec();
}
