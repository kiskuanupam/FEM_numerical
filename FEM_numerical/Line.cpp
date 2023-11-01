#include "Line.h"
#include<math.h>

double Line::calculate_slope()
{
	double slope;
	slope = atan(this->end_point.get_y() - this->start_point.get_y()) / (this->end_point.get_x() - this->start_point.get_x());
	return slope;
}

double Line::calculate_length()
{

	return this->start_point.distance_from_point(this->end_point);

}

Line::Line(Point a, Point b)

{
	this->start_point = a;
	this->end_point = b;

}

