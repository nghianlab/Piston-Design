#ifndef PISTON_H
#define PISTON_H

#include <QtGui/QMainWindow>
#include "ui_piston.h"

class MainGLWidget;
class Piston : public QMainWindow
{
	Q_OBJECT

public:
	Piston(QWidget *parent = 0, Qt::WFlags flags = 0);
	~Piston();

	void setupUi();
private slots: 
	void PistonLengthChange(double val);
	void PistonDiameterChange(double val);
	void PistonThicknessChange(double val);
	void PistonNumRingChange(int val);
	void makePiston();
	void exportPiston();
private:
	Ui::PistonClass ui;
	MainGLWidget *_glWidget;
	QVBoxLayout* _docomLayout;
};

#endif // PISTON_H
