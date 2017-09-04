#include <QApplication>

#include "MainWindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

	app.setApplicationName(app.translate("main", "EnsembleVis"));
	app.setOrganizationName("CG&CAD");
	app.setOrganizationDomain("cg&cad.tsinghua.edu.cn");

	MainWindow mainwindow;
	mainwindow.show();
    return app.exec();
}
