#include "piston.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Piston w;
	w.show();
	return a.exec();
}
