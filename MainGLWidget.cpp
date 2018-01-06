#include "MainGLWidget.h"
#include <QMouseEvent>
#include <QtGui/QtGui>
#include <QDebug>
#include <iostream>
#include <ostream>
#include <strstream>
#include <fstream>
#include "Point3Dd.h"
#include <QMutexLocker>
#define PI_VAL 3.14159265358979
#define MAX_GL_ZOOM 1000000
#define MIN_GL_ZOOM 0.0001

#define MOTION_NUM_TRI_LIMITED 500000
//const unsigned int N_ORIENTATION_AREA = 26 * 4;
using namespace std;


MainGLWidget::MainGLWidget(QGLContext * context, QWidget * parent, const QGLWidget * shareWidget, Qt::WindowFlags f)
                           :QGLWidget(context, parent, shareWidget, f)
{
    setAcceptDrops(true);

    // Preview
    this->setContextMenuPolicy(Qt::CustomContextMenu);

    _renderingCenter[0] = 0.0;
    _renderingCenter[1] = 0.0;
    _renderingCenter[2] = 0.0;
    _renderingScale = 1.0;
	_cx = 0;
	_cy = 0;

	_pisLength = 150.0;
	_pisDiameter = 190.0;
	_pisThickness = 9.0;
	_pisNumRing = 3;
	createPiston();
}

MainGLWidget::~MainGLWidget() {
	clearPistion();
}

QSize MainGLWidget::minimumSizeHint() const {
    return QSize(50, 50);
}

QSize MainGLWidget::sizeHint() const {
    return QSize(400, 400);
}

void MainGLWidget::setDefaultGLView()
{
    updateGL();
}

void MainGLWidget::resizeGL(int width, int height)
{

    if (_cx > 0 && _cy > 0){
        float sx = 0.0f;
        if (static_cast<float>(abs(width - _cx)) / static_cast<float>(_cx) >
            static_cast<float>(abs(height - _cy)) / static_cast<float>(_cy)){
            sx = static_cast<float>(height) / static_cast<float>(_cy);
        }
        else{
            sx = static_cast<float>(width) / static_cast<float>(_cx);
        }
        _renderingScale = (_renderingScale *sx); // [Virtual]
    }
    _cx = width;
    _cy = height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //this->SetupViewFrustum(); // [Virtual]
    SetupViewFrustum();
    glMatrixMode(GL_MODELVIEW);
	resetView();

}

void MainGLWidget::initializeGL()
{
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_POLYGON_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	setupMaterial();
}


int MainGLWidget::setupMaterial()
{
    float   MatAmbient[] = { 0.15f, 0.15f, 0.15f, 1.0f };
    //  float   MatDiffuse[]  = {0.929524f, 0.796542f, 0.178823f,1.0f};//{0.829,0.829,  0.829,  0.922};
    float   MatSpecular[] = { 0.296648f, 0.296648f, 0.296648f, 1.0f };//{0.296648,   0.296648,0.296648,0.922};
    float   MatShininess = 128.0f;
    //float MatEmission[]  = {0.0f, 0.0f, 0.0f, 0.0f};
    float redcolor[] = { 1.0f, 0.0f, 0.0f };
    //float MatAmbientBack[]  = {0.0f, 0.5f, 1.0f, 0.0f}; // green material behind
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, MatAmbient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, redcolor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, MatSpecular);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, MatShininess);

    return 0;
}


void MainGLWidget::paintGL()
{
    makeCurrent();
 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    renderBackground();

    // coordinate system
	glDisable(GL_LIGHTING);
	glColor3f(1.0f,0.0f,0.0f);
	glLineWidth(2.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3f(50.0f/_renderingScale,0.0f,0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f,50.0f/_renderingScale,0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f,0.0f,50.0f/_renderingScale);
	glEnd();
	glEnable(GL_LIGHTING);

	double dd[16];
	glGetDoublev(GL_MODELVIEW_MATRIX,dd);

	QMutexLocker autolock(&_mutex);
	glColor3f(0.7f,0.7f,0.7f);
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < _vtt.size(); ++i){
		glNormal3dv(_vtt[i]->_normal.v);
		glVertex3dv(_vtt[i]->_vpp[0]->v);
		glVertex3dv(_vtt[i]->_vpp[1]->v);
		glVertex3dv(_vtt[i]->_vpp[2]->v);
	}

	glEnd();
    
}


void MainGLWidget::renderBackground() {
    glClearColor(0.85,0.9,0.9,0);
     return;


}


void MainGLWidget::mousePressEvent(QMouseEvent *event)
{

    _lastPos = event->pos();
    _mouseDownPos = _lastPos;
	switch (event->button()) {
	case Qt::LeftButton:{
        
		BeginTracking(_lastPos, GL_TRACKING_ROTATE);
	
		break;
						}
	case Qt::RightButton: {
        
		BeginTracking(_lastPos, GL_TRACKING_PAN);

		break;
						  }
	case Qt::MidButton: {
		
		break;
						}
	default:
		break;
	}
}

void MainGLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    
    QPoint p = event->pos();
    switch (event->button()) {
    case Qt::MidButton:
        //EndTracking();
        break;
    case Qt::RightButton:
    {
        EndTracking();
        break;
    }
    case Qt::LeftButton:{
        
        EndTracking();
        break;
    }
    default:
        break;
    }
    updateGL();
}

void MainGLWidget::mouseMoveEvent(QMouseEvent *event)
{

    _currentMousePos = event->pos();

    if (event->buttons() & Qt::LeftButton ||
        event->buttons() & Qt::RightButton) {
        DoTracking(_currentMousePos);
        updateGL();
    }
    _lastPos = event->pos();

    

}


void MainGLWidget::wheelEvent(QWheelEvent *event)
{
   

    QPoint p = event->pos();
    // Local zoom
    _lastPoint = p;
    QPoint screenCenterPoint(_cx/2, _cy/2);
    MoveScene(screenCenterPoint);

    // Do zooming
    int delta = event->delta();
    BeginTracking(p, GL_TRACKING_ZOOM);
    QPoint dp(p.x() - delta/3, p.y() - delta/3);
    DoTracking(dp);
    EndTracking();

    // Local zoom
    _lastPoint = screenCenterPoint;
    MoveScene(p);
    Sleep(75);

    updateGL();
}

void MainGLWidget::resetView()
{

	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//this->SetupViewFrustum(); // [Virtual]
	SetupViewFrustum();
	glMatrixMode(GL_MODELVIEW);
	double modelview[16] = {0.78582614660263062,-0.32773992419242859,0.52446490526199341,0.00000000000000000,0.61844176054000854,0.41289731860160828,-0.66861474514007568,0.00000000000000000,0.0025815307162702084,0.84976631402969360,0.52715325355529785,0.00000000000000000,	-3.4936046600341797,-53.493606567382812,-8.2887709140777588e-008,1.0000000000000000};//{0.86602540378443837, -0.171010071662, -0.469846310392, 0, 0.5, 0.29619813272602, 0.81379768134937, 0, 0, 0.9396926207, 0.342020143325, 0,0.0, 0.0, 0.0, 1.0};
	glLoadMatrixd(modelview);
	updateGL();
}

void MainGLWidget::BeginTracking(const QPoint& p, const GL_TRACKING_MODE tm) {
    _trackingMode = tm;
    if (GL_TRACKING_NONE != tm)
        _lastPoint = p;
}

void MainGLWidget::DoTracking(const QPoint& p) {
    switch (_trackingMode) {
    case GL_TRACKING_ROTATE: {
        RotateScene(p);
        _lastPoint = p;

        break;
    }
    case GL_TRACKING_ZOOM: {
        ZoomScene(p);
        _lastPoint = p;
        break;
    }
    case GL_TRACKING_PAN: {
        MoveScene(p);
        _lastPoint = p;
        break;
    }
    default:
        break;
    }
}

void MainGLWidget::EndTracking() {
	_trackingMode = GL_TRACKING_NONE;
}

#include <gl/GLU.h>
int MainGLWidget::convertScreenToOpenGLCoord(QPoint point, GLdouble outpos[3],
	bool onNear, GLdouble znear)
{
	makeCurrent();
	glPushMatrix();
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLdouble winX, winY, winZ;

	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);

	glGetIntegerv(GL_VIEWPORT, viewport);
	GLfloat zdep = 0;
	winX = (GLdouble)point.x();
	winY = (GLdouble)viewport[3] - (GLdouble)point.y();

	if (onNear){
		gluUnProject(winX, winY, znear, modelview, projection, viewport, &outpos[0], &outpos[1], &outpos[2]);
	}
	else{
		glReadPixels(point.x(), (viewport[3] - (GLint)point.y()), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zdep);
		winZ = zdep;
		gluUnProject(winX, winY, winZ, modelview, projection, viewport, &outpos[0], &outpos[1], &outpos[2]);
	}
	glPopMatrix();
	return 0;
}


int MainGLWidget::calcScreenNorm(GLdouble norm[3]) {
	glPushMatrix();
	GLdouble modelView[4][4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelView[0]);
	GLdouble iden[3] = { 0.0, 0.0, 1.0 };
	norm[0] = (iden[0] * modelView[0][0] +
		iden[1] * modelView[0][1] +
		iden[2] * modelView[0][2]);
	norm[1] = (iden[0] * modelView[1][0] +
		iden[1] * modelView[1][1] +
		iden[2] * modelView[1][2]);
	norm[2] = (iden[0] * modelView[2][0] +
		iden[1] * modelView[2][1] +
		iden[2] * modelView[2][2]);
	glPopMatrix();
	return 0;
}

void MainGLWidget::RotateScene(const QPoint& p) {
	Point3Dd sp,ep;
	convertScreenToOpenGLCoord(_lastPoint, sp.v);
	convertScreenToOpenGLCoord(p,ep.v);
	Point3Dd norm;
	calcScreenNorm(norm.v);
	Point3Dd org(0.0,0.0,0.0);
	Point3Dd sp1,ep1;
	sp.orthoToPlane(org,norm,sp1);
	ep.orthoToPlane(org,norm,ep1);
	Point3Dd dv1(ep1);
	dv1 -= sp1;

	Point3Dd rotAxis = norm*dv1;
	rotAxis.unit();

	int dx = p.x() - _lastPoint.x();
	int dy = p.y() - _lastPoint.y();
	double leng = sqrt(double(dx*dx + dy*dy));
	double bigleng = sqrt(double(_cx*_cx + _cy*_cy));
	double d_a = leng * 855.0 / bigleng;

	if (d_a > 360.0) d_a -= 360.0;
	if (d_a < -360.0) d_a += 360.0;
	double r_a = d_a*PI_VAL/180.0;
	glMatrixMode(GL_MODELVIEW);
	double m[16];
	double rotCenter[3] = {0.0,0.0,0.0};
	calcRotMatrix(r_a,rotCenter/*center.v*/,rotAxis.v,m);
	glMultMatrixd(m);
}

void MainGLWidget::calcRotMatrix(double angle, double *p,
	double *dir, double *m)
{
	double c_ = cos(angle);
	double s_ = sin(angle);
	double a = p[0];
	double b = p[1];
	double c = p[2];
	double u = dir[0];
	double v = dir[1];
	double w = dir[2];
	double u2 = u*u;
	double v2 = v*v;
	double w2 = w*w;

	m[0] = u2 + (v2 + w2)*c_;
	m[4] = u*v*(1 - c_) - w*s_;
	m[8] = u*w*(1 - c_) + v*s_;
	m[12] = (a*(v2 + w2) - u*(b*v + c*w))*(1 - c_) + (b*w -c*v)*s_;
	m[1] = u*v*(1 - c_) + w*s_;
	m[5] = v2 + (u2 + w2)*c_;
	m[9] = v*w*(1 - c_) - u*s_;
	m[13] = (b*(u2 + w2) - v*(a*u + c*w))*(1 - c_) + (c*u - a*w)*s_;
	m[2] = u*w*(1 - c_) - v*s_;
	m[6] = v*w*(1 - c_) + u*s_;
	m[10] = w2 + (u2 + v2)*c_;
	m[14] = (c*(u2 + v2) - w*(a*u + b*v))*(1 - c_) + (a*v - b*u)*s_;
	m[3] = 0.0;
	m[7] = 0.0;
	m[11] = 0.0;
	m[15] = 1.0;
}

void MainGLWidget::ZoomScene(const QPoint& p) {
    _renderingScale = _renderingScale*(float)pow(2.0, (_lastPoint.y() - p.y()) * 0.005);

    if (_renderingScale < MIN_GL_ZOOM){
        _renderingScale = MIN_GL_ZOOM;
		//return;
	}

    if (_renderingScale > MAX_GL_ZOOM){
        _renderingScale = MAX_GL_ZOOM;
		//return;
	}

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	SetupViewFrustum();

	glMatrixMode(GL_MODELVIEW);
}


void MainGLWidget::MoveScene(const QPoint& p,bool iscalcCenter) {
   
	 Point3Dd sp,ep;
	convertScreenToOpenGLCoord(_lastPoint, sp.v);
	convertScreenToOpenGLCoord(p,ep.v);
	Point3Dd norm;
	calcScreenNorm(norm.v);
	Point3Dd org(0.0,0.0,0.0);
	Point3Dd sp1,ep1;
	sp.orthoToPlane(org,norm,sp1);
	ep.orthoToPlane(org,norm,ep1);
	Point3Dd dv1(ep1);
	dv1 -= sp1;
	glTranslatef((float)dv1.v[0],(float)dv1.v[1],(float)dv1.v[2]);
   
}


void MainGLWidget::SetupViewFrustum() {
    

    double sizeBox = 5000.0;
    double minCoord[3] = { -sizeBox, -sizeBox, -sizeBox }, maxCoord[3] = { sizeBox, sizeBox, sizeBox };
    


    double left = -(double)_cx * 0.5 / _renderingScale; // [Virtual]
    double right = (double)_cx * 0.5 / _renderingScale;
    double top = (double)_cy * 0.5 / _renderingScale;
    double bottom = -(double)_cy * 0.5 / _renderingScale;


    minCoord[0] += _renderingCenter[0];
    minCoord[1] += _renderingCenter[1];
    minCoord[2] += _renderingCenter[2];
    maxCoord[0] += _renderingCenter[0];
    maxCoord[1] += _renderingCenter[1];
    maxCoord[2] += _renderingCenter[2];

    double viewRadius = max(right - left, top - bottom);
    double zNear, zFar;

    zNear = min(min(minCoord[0], minCoord[1]), minCoord[2]);
    zFar = max(max(maxCoord[0], maxCoord[1]), maxCoord[2]);
    double center = (zNear + zFar) / 2;
    zNear = 2 * zNear - center;
    zFar = 2 * zFar - center;

    if ((zFar - zNear) < viewRadius*2.0)
    {
        zNear = center - viewRadius;
        zFar = center + viewRadius;
    }

    
    zFar *= 2.0;
    
    zNear = -zFar*2.0;

    glOrtho(left, right, bottom, top, zNear, zFar);
    return;
}

void MainGLWidget::createPiston()
{
	QMutexLocker autolock(&_mutex);
	clearPistion();

	double pisRingSpace = 6.0;
	double pisRingWidth = 6.0;
	double pisRingDepth = _pisThickness*0.6;
	double pisDistanceToFirstRing = _pisLength - 20.0;
	double pisHoleRadius = 25.0;
	double pisDistanceToHole = pisDistanceToFirstRing - pisHoleRadius - _pisNumRing*(pisRingSpace + pisRingWidth)-5.0;


	int resol = 90;
	double stepAngle = 2.0*PI_VAL/(double)resol;
	double outterRad = _pisDiameter*0.5;
	double innerRad = outterRad - _pisThickness;
	double ringDepRad = outterRad - pisRingDepth;
	double innerLength = _pisLength - _pisThickness;
	std::vector<Point3Dd*> vppBottomOut;
	std::vector<Point3Dd*> vppTopOut;
	std::vector<Point3Dd*> vppBottomInner;
	std::vector<Point3Dd*> vppTopInner;
	std::vector<Point3Dd*> vppRingInDepth;
	vppBottomOut.resize(resol);
	vppTopOut.resize(resol);
	vppTopInner.resize(resol);
	vppBottomInner.resize(resol);
	vppRingInDepth.resize(resol);
	double angle = 0.0;
	double cAngle,sAngle;
	for(int i = 0; i < resol; ++i){
		cAngle = cos(angle);
		sAngle = sin(angle);
		vppBottomOut[i] = new Point3Dd(outterRad*cAngle, outterRad*sAngle,0.0);
		vppTopOut[i] = new Point3Dd(outterRad*cAngle, outterRad*sAngle,_pisLength);
		vppBottomInner[i] = new Point3Dd(innerRad*cAngle, innerRad*sAngle,0.0);
		vppTopInner[i] = new Point3Dd(innerRad*cAngle, innerRad*sAngle,innerLength);
		vppRingInDepth[i] = new Point3Dd(ringDepRad*cAngle, ringDepRad*sAngle,pisDistanceToFirstRing);
		angle += stepAngle;
	}
	_vpp.insert(_vpp.end(),vppBottomOut.begin(),vppBottomOut.end());
	_vpp.insert(_vpp.end(),vppTopOut.begin(),vppTopOut.end());
	_vpp.insert(_vpp.end(),vppBottomInner.begin(),vppBottomInner.end());
	_vpp.insert(_vpp.end(),vppTopInner.begin(),vppTopInner.end());
	// topCap
	Point3Dd *topCenter = new Point3Dd(0,0,_pisLength);
	Point3Dd *topCenterIn = new Point3Dd(0.0,0.0,innerLength);
	_vpp.push_back(topCenterIn);
	_vpp.push_back(topCenter);
	Triangle *trig;
	for(int i = 0; i < resol; ++i){
		trig = new Triangle(topCenter,vppTopOut[i], vppTopOut[(i+1)%resol]);
		_vtt.push_back(trig);
		trig = new Triangle(topCenterIn,vppTopInner[(i+1)%resol], vppTopInner[i]);
		_vtt.push_back(trig);
	}

	// bottom annular
	for(int i = 0; i < resol; ++i){
		trig = new Triangle(vppBottomInner[i], vppBottomOut[(i+1)%resol],vppBottomOut[i]);
		_vtt.push_back(trig);
		trig = new Triangle(vppBottomInner[i], vppBottomInner[(i+1)%resol],vppBottomOut[(i+1)%resol]);
		_vtt.push_back(trig);
	}

	//Rings
	std::vector<Point3Dd*> vpp1,vpp2,vpp3,vpp4,vpp5;
	vpp2.resize(resol);
	vpp3.resize(resol);
	vpp4.resize(resol);
	vpp5.resize(resol);
	vpp1.assign(vppTopOut.begin(),vppTopOut.end());
	double h0 = pisDistanceToFirstRing;
	for(int j = 0; j < _pisNumRing; ++j){
		for(int i = 0; i < resol; ++i){
			vpp2[i] = new Point3Dd(vpp1[i]->v[0],vpp1[i]->v[1],h0);
			vpp3[i] = new Point3Dd(vppRingInDepth[i]->v[0],vppRingInDepth[i]->v[1],h0);
			vpp4[i] = new Point3Dd(vppRingInDepth[i]->v[0],vppRingInDepth[i]->v[1],h0-pisRingWidth);
			vpp5[i] = new Point3Dd(vpp1[i]->v[0],vpp1[i]->v[1],h0-pisRingWidth);
		}
		_vpp.insert(_vpp.end(),vpp2.begin(),vpp2.end());
		_vpp.insert(_vpp.end(),vpp3.begin(),vpp3.end());
		_vpp.insert(_vpp.end(),vpp4.begin(),vpp4.end());
		_vpp.insert(_vpp.end(),vpp5.begin(),vpp5.end());

		for(int i = 0; i < resol; ++i){
			trig = new Triangle(vpp1[i], vpp2[i],vpp2[(i+1)%resol]);
			_vtt.push_back(trig);
			trig = new Triangle(vpp1[i], vpp2[(i+1)%resol],vpp1[(i+1)%resol]);
			_vtt.push_back(trig);

			trig = new Triangle(vpp2[i], vpp3[i],vpp3[(i+1)%resol]);
			_vtt.push_back(trig);
			trig = new Triangle(vpp2[i], vpp3[(i+1)%resol],vpp2[(i+1)%resol]);
			_vtt.push_back(trig);

			trig = new Triangle(vpp3[i], vpp4[i],vpp4[(i+1)%resol]);
			_vtt.push_back(trig);
			trig = new Triangle(vpp3[i], vpp4[(i+1)%resol],vpp3[(i+1)%resol]);
			_vtt.push_back(trig);

			trig = new Triangle(vpp4[i], vpp5[i],vpp5[(i+1)%resol]);
			_vtt.push_back(trig);
			trig = new Triangle(vpp4[i], vpp5[(i+1)%resol],vpp4[(i+1)%resol]);
			_vtt.push_back(trig);

			
			
		}
		h0 -= (pisRingWidth + pisRingSpace);
		vpp1.assign(vpp5.begin(),vpp5.end());
	}

	// Circle on ZX
	double xcoord, ycoord,zcoord;
	angle = 0.0;
	for(int i = 0; i < resol; ++i){
		xcoord = pisHoleRadius*cos(angle);
		zcoord = pisDistanceToHole + pisHoleRadius*sin(angle);
		ycoord = sqrt(outterRad*outterRad - xcoord*xcoord);
		vpp1[i] = new Point3Dd(xcoord,-ycoord,zcoord);
		vpp4[i] = new Point3Dd(xcoord,ycoord,zcoord);
		ycoord = sqrt(innerRad*innerRad - xcoord*xcoord);
		vpp2[i] = new Point3Dd(xcoord,-ycoord,zcoord);
		vpp3[i] = new Point3Dd(xcoord,ycoord,zcoord);
		angle += stepAngle;
	}

	_vpp.insert(_vpp.end(),vpp1.begin(),vpp1.end());
	_vpp.insert(_vpp.end(),vpp2.begin(),vpp2.end());
	_vpp.insert(_vpp.end(),vpp3.begin(),vpp3.end());
	_vpp.insert(_vpp.end(),vpp4.begin(),vpp4.end());

	// hole wall
	for(int i = 0; i < resol; ++i){
		
		trig = new Triangle(vpp1[i],vpp2[(i+1)%resol], vpp2[i]);
		_vtt.push_back(trig);
		trig = new Triangle(vpp1[i],vpp1[(i+1)%resol], vpp2[(i+1)%resol]);
		_vtt.push_back(trig);


		trig = new Triangle(vpp3[i],vpp4[(i+1)%resol], vpp4[i]);
		_vtt.push_back(trig);
		trig = new Triangle(vpp3[i],vpp3[(i+1)%resol], vpp4[(i+1)%resol]);
		_vtt.push_back(trig);
	}

    // out wall
	int sInd1 = -1,  sInd2 = -1;
	for(int i = 0; i < resol; ++i){
		if((vpp5[i]->v[0] > pisHoleRadius && vpp5[(i+1)%resol]->v[0] > pisHoleRadius) ||
			(vpp5[i]->v[0] < -pisHoleRadius && vpp5[(i+1)%resol]->v[0] < -pisHoleRadius)){
			trig = new Triangle(vpp5[i], vppBottomOut[i],vppBottomOut[(i+1)%resol]);
			_vtt.push_back(trig);
			trig = new Triangle(vpp5[i], vppBottomOut[(i+1)%resol],vpp5[(i+1)%resol]);
			_vtt.push_back(trig);
		}
		else if(-1 == sInd1){
			sInd1 = i;
		}

		if((vppBottomInner[i]->v[0] > pisHoleRadius && vppBottomInner[(i+1)%resol]->v[0] > pisHoleRadius) ||
			(vppBottomInner[i]->v[0] < -pisHoleRadius && vppBottomInner[(i+1)%resol]->v[0] < -pisHoleRadius)){
				trig = new Triangle(vppBottomInner[i], vppTopInner[i],vppTopInner[(i+1)%resol]);
				_vtt.push_back(trig);
				trig = new Triangle(vppBottomInner[i], vppTopInner[(i+1)%resol],vppBottomInner[(i+1)%resol]);
				_vtt.push_back(trig);
		}
		else if(sInd2 == -1){
			sInd2 = i;
		}
	}

	// around hole
	Point3Dd *Pbo_1 = vpp5[sInd1], *Pbo_2 = vppBottomOut[sInd1],*ph0,*ph1,*pu0,*pu1;
	Point3Dd *Pbi_1 = vppTopInner[sInd2], *Pbi_2 = vppBottomInner[sInd2],*vh0,*vh1,*vu0,*vu1;

	trig = new Triangle(Pbo_1,Pbo_2, vpp4[0]);
	_vtt.push_back(trig);
	trig = new Triangle(Pbi_1,Pbi_2, vpp3[0]);
	_vtt.push_back(trig);
	int resol2 = resol/2;
	int i1 = sInd1;
	int i2 = sInd2;
	ph0 = vpp4[0];
	pu0 = ph0;
	vh0 = vpp3[0];
	vu0 = vh0;
	for(int i = 0; i < resol2; ++i){
		ph1 = vpp4[i+1];
		trig = new Triangle(Pbo_1,ph0,ph1);
		_vtt.push_back(trig);
		ph0 = ph1;

		pu1 = vpp4[resol - i -1];
		trig = new Triangle(Pbo_2,pu1,pu0);
		_vtt.push_back(trig);
		pu0 = pu1;

		// inner
		vh1 = vpp3[i+1];
		trig = new Triangle(Pbi_1,vh0,vh1);
		_vtt.push_back(trig);
		vh0 = vh1;

		vu1 = vpp3[resol - i -1];
		trig = new Triangle(Pbi_2,vu1,vu0);
		_vtt.push_back(trig);
		vu0 = vu1;

		if(i == (resol2 -1)){
			if(fabs(Pbo_1->v[0]) < pisHoleRadius){
				trig = new Triangle(Pbo_1,ph1,vpp5[i1+1]);
				_vtt.push_back(trig);

				trig = new Triangle(Pbo_2,vppBottomOut[i1+1],pu1);
				_vtt.push_back(trig);
				i1++;
				Pbo_1 = vpp5[i1];
				Pbo_2 = vppBottomOut[i1];
			}
			trig = new Triangle(Pbo_2,Pbo_1,ph1);
			_vtt.push_back(trig);

			//inner
			if(fabs(Pbi_1->v[0]) < pisHoleRadius){
				trig = new Triangle(Pbi_1,vh1,vppTopInner[i2+1]);
				_vtt.push_back(trig);

				trig = new Triangle(Pbi_2,vppBottomInner[i2+1],vu1);
				_vtt.push_back(trig);
				i2++;
				Pbi_1 = vppTopInner[i2];
				Pbi_2 = vppBottomInner[i2];
			}
			trig = new Triangle(Pbi_2,Pbi_1,vh1);
			_vtt.push_back(trig);
		}
		
		if(abs(ph1->v[0] - Pbo_1->v[0]) > abs(ph1->v[0] - vpp5[i1+1]->v[0])){
			trig = new Triangle(Pbo_1,ph1,vpp5[i1+1]);
			_vtt.push_back(trig);

			trig = new Triangle(Pbo_2,vppBottomOut[i1+1],pu1);
			_vtt.push_back(trig);

			i1++;
			Pbo_1 = vpp5[i1];
			Pbo_2 = vppBottomOut[i1];
		}
		//inner
		if(abs(vh1->v[0] - Pbi_1->v[0]) > abs(vh1->v[0] - vppTopInner[i2+1]->v[0])){
			trig = new Triangle(Pbi_1,vh1,vppTopInner[i2+1]);
			_vtt.push_back(trig);

			trig = new Triangle(Pbi_2,vppBottomInner[i2+1],vu1);
			_vtt.push_back(trig);

			i2++;
			Pbi_1 = vppTopInner[i2];
			Pbi_2 = vppBottomInner[i2];
		}
	}

	sInd1 = resol - sInd1;
	sInd2 = resol - sInd2;
	Pbo_1 = vpp5[sInd1];
	Pbo_2 = vppBottomOut[sInd1];
	Pbi_1 = vppTopInner[sInd2];
	Pbi_2 = vppBottomInner[sInd2];
	i1 = sInd1;
	i2 = sInd2;
	ph0 = vpp1[0];
	pu0 = ph0;
	vh0 = vpp2[0];
	vu0 = vh0;

	trig = new Triangle(Pbo_2,Pbo_1,ph0);
	_vtt.push_back(trig);
	trig = new Triangle(Pbi_2,Pbi_1, vh0);
	_vtt.push_back(trig);
	
	for(int i = 0; i < resol2; ++i){
		ph1 = vpp1[i+1];
		trig = new Triangle(Pbo_1,ph1,ph0);
		_vtt.push_back(trig);
		ph0 = ph1;

		pu1 = vpp1[resol - i -1];
		trig = new Triangle(Pbo_2,pu0,pu1);
		_vtt.push_back(trig);
		pu0 = pu1;

		// inner
		vh1 = vpp2[i+1];
		trig = new Triangle(Pbi_1,vh1,vh0);
		_vtt.push_back(trig);
		vh0 = vh1;

		vu1 = vpp2[resol - i -1];
		trig = new Triangle(Pbi_2,vu0,vu1);
		_vtt.push_back(trig);
		vu0 = vu1;

		if(i == (resol2 -1)){
			if(fabs(Pbo_1->v[0]) < pisHoleRadius){
				trig = new Triangle(Pbo_1,vpp5[i1-1],ph1);
				_vtt.push_back(trig);

				trig = new Triangle(Pbo_2,pu1,vppBottomOut[i1-1]);
				_vtt.push_back(trig);
				i1--;
				Pbo_1 = vpp5[i1];
				Pbo_2 = vppBottomOut[i1];
			}
			trig = new Triangle(Pbo_2,ph1,Pbo_1);
			_vtt.push_back(trig);

			//inner
			if(fabs(Pbi_1->v[0]) < pisHoleRadius){
				trig = new Triangle(Pbi_1,vppTopInner[i2-1],vh1);
				_vtt.push_back(trig);

				trig = new Triangle(Pbi_2,vu1,vppBottomInner[i2-1]);
				_vtt.push_back(trig);
				i2--;
				Pbi_1 = vppTopInner[i2];
				Pbi_2 = vppBottomInner[i2];
			}
			trig = new Triangle(Pbi_2,vh1,Pbi_1);
			_vtt.push_back(trig);
		}

		if(abs(ph1->v[0] - Pbo_1->v[0]) > abs(ph1->v[0] - vpp5[i1-1]->v[0])){
			trig = new Triangle(Pbo_1,vpp5[i1-1],ph1);
			_vtt.push_back(trig);

			trig = new Triangle(Pbo_2,pu1,vppBottomOut[i1-1]);
			_vtt.push_back(trig);

			i1--;
			Pbo_1 = vpp5[i1];
			Pbo_2 = vppBottomOut[i1];
		}
		//inner
		if(abs(vh1->v[0] - Pbi_1->v[0]) > abs(vh1->v[0] - vppTopInner[i2-1]->v[0])){
			trig = new Triangle(Pbi_1,vppTopInner[i2-1],vh1);
			_vtt.push_back(trig);

			trig = new Triangle(Pbi_2,vu1,vppBottomInner[i2-1]);
			_vtt.push_back(trig);

			i2--;
			Pbi_1 = vppTopInner[i2];
			Pbi_2 = vppBottomInner[i2];
		}
	}
	//resetView();
}

void MainGLWidget::clearPistion()
{
	for(auto it = _vpp.begin(); it != _vpp.end(); ++it){
		delete (*it);
	}
	_vpp.clear();

	for(auto it = _vtt.begin(); it != _vtt.end(); ++it){
		delete (*it);
	}
	_vtt.clear();
}

void MainGLWidget::savePiston()
{
	if(_vtt.empty())
		return;

	QString filterStlBin = tr("STL-STereoLithography File Format, Binary (*.stl)");
	QString filterSel;
	QFileDialog saveDialog;

	QString fileNamePath = saveDialog.getSaveFileName(nullptr, tr("Save Object"), "", filterStlBin, &filterSel);
	// ------------------------------------------------------------------------------------------------------------
	if (fileNamePath.isEmpty())
	{
		return;
	}

	ofstream file(fileNamePath.toStdString(), ios::binary | ios::out );

	if (!file.is_open())
		return ;
	char headerFile[88] = "";
	file.write(headerFile, 80);
	int nts = _vtt.size();
	file.write(reinterpret_cast<char *>(&nts), sizeof(nts));
	Triangle *trig;
	for (int i = 0; i < nts; ++i) {
		trig = _vtt[i];
		writeBytesFromFloat(file, (float)(trig->_normal.v[0]));
		writeBytesFromFloat(file, (float)(trig->_normal.v[1]));
		writeBytesFromFloat(file, (float)(trig->_normal.v[2]));

		for (int j = 0; j < 3; j ++) {
			writeBytesFromFloat(file, (float)(trig->_vpp[j]->v[0]));
			writeBytesFromFloat(file, (float)(trig->_vpp[j]->v[1]));
			writeBytesFromFloat(file, (float)(trig->_vpp[j]->v[2]));
		}
		file.put(0);
		file.put(0);
	}
}

void writeBytesFromFloat(::std::ofstream &file, float valueIn)
{
	union {
		float floatValue;
		char charValue[4];
	} value;
	value.floatValue = valueIn;
	int newValue  = value.charValue[0] & 0xFF;
	newValue |= (value.charValue[1] & 0xFF) << 0x08;
	newValue |= (value.charValue[2] & 0xFF) << 0x10;
	newValue |= (value.charValue[3] & 0xFF) << 0x18;
	file.write(reinterpret_cast<char *>(&newValue), sizeof(newValue));
}