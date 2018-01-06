#include "Point3Dd.h"
#include <math.h>


Point3Dd::Point3Dd()
{
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

Point3Dd::Point3Dd(const double *mv)
{
	v[0] = mv[0];
	v[1] = mv[1];
	v[2] = mv[2];
}

Point3Dd::Point3Dd(const double& mx, const double& my, const double& mz)
{
	v[0] = mx;
	v[1] = my;
	v[2] = mz;
}

void Point3Dd::set(const double& mx, const double& my, const double& mz)
{
	v[0] = mx;
	v[1] = my;
	v[2] = mz;
}

void Point3Dd::set(const double& val)
{
	v[0] = val;
	v[1] = val;
	v[2] = val;
}


Point3Dd::Point3Dd(const Point3Dd& other)
{
	v[0] = other.v[0];
	v[1] = other.v[1];
	v[2] = other.v[2];
}


Point3Dd::~Point3Dd(void)
{
}


Point3Dd& Point3Dd::operator = (const Point3Dd& other) {
	if (this != (&other)) {
		v[0] = other.v[0], v[1] = other.v[1], v[2] = other.v[2];
	}
	return *this;
}

Point3Dd& Point3Dd::operator = (const double other[3])
{

	v[0] = other[0], v[1] = other[1], v[2] = other[2];

	return *this;
}

bool Point3Dd::operator == (const Point3Dd& v1) const
{
	if( (fabs(v[0] - v1.v[0]) <= 1e-7) && (fabs(v[1] - v1.v[1]) <= 1e-7) && (fabs(v[2] - v1.v[2]) <= 1e-7) ) {
		return true;
	} else {
		double l = sqrt(((v[0] - v1.v[0])*(v[0] - v1.v[0])) + ((v[1] - v1.v[1])*(v[1] - v1.v[1])) + ((v[2] - v1.v[2])*(v[2] - v1.v[2])));
		if( l <= 1e-7 )
			return true;
	}
	return false;
}

bool Point3Dd::operator != (const Point3Dd& v1) const
{
	if( (*this) == v1)
		return false;
	else
		return true;
}

Point3Dd Point3Dd::operator +(const Point3Dd& v1) const {
	Point3Dd temp;
	temp.v[0] = v[0] + v1.v[0];
	temp.v[1] = v[1] + v1.v[1];
	temp.v[2] = v[2] + v1.v[2];
	return temp;
}

const Point3Dd& Point3Dd::operator +=(const Point3Dd& v1) {
	v[0] += v1.v[0];
	v[1] += v1.v[1];
	v[2] += v1.v[2];
	return *this;
}

Point3Dd Point3Dd::operator +(double d) const {
	Point3Dd temp;
	temp.v[0] = v[0] + d;
	temp.v[1] = v[1] + d;
	temp.v[2] = v[2] + d;
	return temp;
}

const Point3Dd& Point3Dd::operator +=(double d) {
	v[0] += d;
	v[1] += d;
	v[2] += d;
	return *this;
}

Point3Dd Point3Dd::operator -() const {
	Point3Dd temp;
	temp.v[0] = -v[0];
	temp.v[1] = -v[1] ;
	temp.v[2] = -v[2] ;
	return temp;
}

Point3Dd Point3Dd::operator -(const Point3Dd& v1) const {
	Point3Dd temp;
	temp.v[0] = v[0] - v1.v[0];
	temp.v[1] = v[1] - v1.v[1];
	temp.v[2] = v[2] - v1.v[2];
	return temp;
}

const Point3Dd& Point3Dd::operator -=(const Point3Dd& v1){
	v[0] -= v1.v[0];
	v[1] -= v1.v[1];
	v[2] -= v1.v[2];
	return *this;
}

Point3Dd Point3Dd::operator -(double d) const {
	Point3Dd temp;
	temp.v[0] = v[0] - d;
	temp.v[1] = v[1] - d;
	temp.v[2] = v[2] - d;
	return temp;
}

const Point3Dd& Point3Dd::operator -=(double d){
	v[0] -= d;
	v[1] -= d;
	v[2] -= d;
	return *this;
}


Point3Dd Point3Dd::operator *(double d) const {
	Point3Dd temp;
	temp.v[0] = v[0] * d;
	temp.v[1] = v[1] * d;
	temp.v[2] = v[2] * d;
	return temp;
}


const Point3Dd& Point3Dd::operator *=(double d) {
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
	return *this;
}

Point3Dd Point3Dd::operator *(const Point3Dd& v1) const {
	Point3Dd temp;
	temp.v[0] = v[1] * v1.v[2] - v[2] * v1.v[1];
	temp.v[1] = v[2] * v1.v[0] - v[0] * v1.v[2];
	temp.v[2] = v[0] * v1.v[1] - v[1] * v1.v[0];
	return temp;
}



Point3Dd Point3Dd::operator /(double d) const {
	Point3Dd temp;
	temp.v[0] = v[0] / d;
	temp.v[1] = v[1] / d;
	temp.v[2] = v[2] / d;
	return temp;
}

const Point3Dd& Point3Dd::operator /=(double d){
	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
	return *this;
}

double Point3Dd::unit()
{
	double dd = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if(dd < 1e-15){
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
		return 0.0;
	} else if (fabs(dd - 1.0) < 1e-12) {
		return 1.0;
	}
	double d = sqrt(dd);
	v[0] /= d;
	v[1] /= d;
	v[2] /= d;
	return d;
}

double Point3Dd::scalarProduct(const Point3Dd &v1) const
{
	double d = v[0] * v1.v[0] + v[1] * v1.v[1] + v[2] * v1.v[2];
	return d;
}

int Point3Dd::orthoToPlane(const Point3Dd &pv, const Point3Dd &nv,
	Point3Dd &orth) const
{
	Point3Dd v = *this;
	v -= pv;
	double dis2 = v.scalarProduct(v);
	if (dis2 < 1e-10)
	{
		orth = pv;
		return 0;
	}

	if (fabs(v.scalarProduct(nv)) < 1e-7)
	{
		orth = *this;
		return 0;
	}

	double num = v.scalarProduct(nv);
	double den = nv.scalarProduct(nv);
	//num = pv.scalarProduct(nv)-num;  ko can boi ben tren da co v-=pv
	double t = -num/den;
	orth = *this + nv*t;

	return 0;
}