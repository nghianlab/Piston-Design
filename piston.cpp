#include "piston.h"
#include "MainGLWidget.h"

Piston::Piston(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	this->setupUi();
}

Piston::~Piston()
{
	delete _glWidget;
}

void Piston::setupUi()
{
	ui.setupUi(this);
	_glWidget = new MainGLWidget(nullptr,this->centralWidget());
	_docomLayout = new QVBoxLayout(this->centralWidget());
	_docomLayout->setContentsMargins(0, 0, 0, 0);
	_docomLayout->addWidget(_glWidget);
	ui._glWidget->setLayout(_docomLayout);
	
	connect(ui._dsbLength, SIGNAL(valueChanged(double)), this, SLOT(PistonLengthChange(double)));
	connect(ui._dsbDiameter, SIGNAL(valueChanged(double)), this, SLOT(PistonDiameterChange(double)));
	connect(ui._dsbThickness, SIGNAL(valueChanged(double)), this, SLOT(PistonThicknessChange(double)));
	connect(ui._spbNumberRing, SIGNAL(valueChanged(int)), this, SLOT(PistonNumRingChange(int)));
	connect(ui._bntMake, SIGNAL(clicked()), this, SLOT(makePiston()));
	connect(ui._bntExport, SIGNAL(clicked()), this, SLOT(exportPiston()));

	QImage image((":/Piston/Resources/ReferenceImage.PNG"));
	if (!image.isNull()) {
		
		ui._referenceImage->setPixmap(QPixmap::fromImage(image));
		
	}
	//imageLabel->setPixmap(QPixmap::fromImage(image));
}

void Piston::PistonLengthChange(double val)
{
	_glWidget->setPistonLength(val);
}

void Piston::PistonDiameterChange(double val)
{
	_glWidget->setPistonDiameter(val);
}

void Piston::PistonThicknessChange(double val)
{
	_glWidget->setPistonThickness(val);
}

void Piston::PistonNumRingChange(int val)
{
	_glWidget->setPistonNumRing(val);
}

void Piston::makePiston()
{
	_glWidget->createPiston();
	_glWidget->updateGL();
}

void Piston::exportPiston()
{
	_glWidget->savePiston();
}
