#include "Point.h"
#include<cmath>

void Point::set_xy(double x, double y)
{
	this->x = x;
	this->y = y;
}

double Point::get_x()
{
	return x;
}

double Point::get_y()
{
	return y;
}

double Point::distance_from_point(Point pt)
{

	return sqrt(pow((x - pt.x), 2) + pow((y - pt.y), 2));
}
