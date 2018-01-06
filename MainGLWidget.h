
#include <QGLWidget>
#include <vector>
#include <QMutex>

class Point3Dd;
class Triangle;

class MainGLWidget : public QGLWidget
{
    Q_OBJECT

public:


	enum GL_TRACKING_MODE{
		GL_TRACKING_NONE,
		GL_TRACKING_ROTATE,
		GL_TRACKING_ZOOM,
		GL_TRACKING_PAN
	};



public:
    // MainWindows' MainWidget's GL Widget:
    MainGLWidget(QGLContext * context, QWidget * parent, const QGLWidget * shareWidget = 0, Qt::WindowFlags f = 0); 
    ~MainGLWidget();

public:
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void setDefaultGLView();

protected:
    void paintGL();
    void renderBackground();
    void resizeGL(int width, int height);
    void initializeGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
	void resetView();
public:

	void BeginTracking(const QPoint& p, const GL_TRACKING_MODE tm);
	void DoTracking(const QPoint& p);
	void EndTracking();
	int convertScreenToOpenGLCoord(QPoint point, GLdouble outpos[3],
		bool onNear = true, GLdouble znear = 0.001);
	int calcScreenNorm(GLdouble norm[3]);
	void RotateScene(const QPoint& p);
	void ZoomScene(const QPoint& p);
	void MoveScene(const QPoint& p,bool iscalcCenter = true);
	void calcRotMatrix(double angle, double *p,
		double *dir, double *m);
    void SetupViewFrustum();

    int setupMaterial();

	void setPistonLength(double length){ _pisLength = length;}
	void setPistonDiameter(double diameter){ _pisDiameter = diameter;}
	void setPistonThickness(double thickness){ _pisThickness = thickness;}
	void setPistonNumRing(int n){ _pisNumRing = n;}
	void createPiston();
	void clearPistion();
	void savePiston();
signals:
    void updateMainWidget();

private:

	GL_TRACKING_MODE _trackingMode;  // mouse tracking mode

	QMutex _mutex;
    int _cx, _cy;
    double _renderingScale;
    double _renderingCenter[3];
    QPoint _lastPos;
    QPoint _lastPoint;
    QPoint _currentMousePos;
    QPoint _mouseDownPos;

	///--------------------
	double _pisLength;
	double _pisDiameter;
	double _pisThickness;
	double _pisNumRing;
	std::vector<Point3Dd*> _vpp;
	std::vector<Triangle*> _vtt;
	
};

void writeBytesFromFloat(::std::ofstream &file, float valueIn);