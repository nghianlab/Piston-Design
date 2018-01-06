#pragma once
class Point3Dd
{
public:
	double v[3];
public:
	Point3Dd(void);
	Point3Dd(const double *mv);
	Point3Dd(const double& mx, const double& my, const double& mz);
	Point3Dd(const Point3Dd& other) ;
	~Point3Dd(void);

	Point3Dd& operator = (const Point3Dd& other);
	Point3Dd& operator = (const double other[3]);
	bool operator == (const Point3Dd& v1) const;
	bool operator != (const Point3Dd& v1) const;
	Point3Dd operator +(const Point3Dd& v1) const;
	const Point3Dd& operator +=(const Point3Dd& v1);
	Point3Dd operator +(double d) const;
	const Point3Dd& operator +=(double d);
	Point3Dd operator -() const;
	Point3Dd operator -(const Point3Dd& v1) const;
	const Point3Dd& operator -=(const Point3Dd& v1);
	Point3Dd operator -(double d) const;
	const Point3Dd& operator -=(double d);
	Point3Dd operator *(double d) const;
	const Point3Dd& operator *=(double d);
	Point3Dd operator *(const Point3Dd& v1) const;
	Point3Dd operator /(double d) const;
	const Point3Dd& operator /=(double d);

	void set(const double& mx, const double& my, const double& mz);
	void set(const double& val);
	double unit();
	double scalarProduct(const Point3Dd &v1) const;
	int orthoToPlane(const Point3Dd &pv, const Point3Dd &nv,
		Point3Dd &orth) const;
};


class Triangle{
public:
	Point3Dd* _vpp[3];
	Point3Dd _normal;
	Triangle(){
		_vpp[0] = nullptr;	
		_vpp[1] = nullptr;	
		_vpp[2] = nullptr;	
	}
	Triangle(Point3Dd *p1, Point3Dd *p2, Point3Dd *p3){
		_vpp[0] = p1;
		_vpp[1] = p2;
		_vpp[2] = p3;
		calcNormal();
	}
	~Triangle(){};

	void calcNormal(){
		Point3Dd dv1(*_vpp[1]);
		dv1 -= *_vpp[0];
		Point3Dd dv2(*_vpp[2]);
		dv2 -= *_vpp[0];
		
		_normal = dv1*dv2;

		_normal.unit();
	}
};