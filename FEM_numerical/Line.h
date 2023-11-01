#pragma once
#include"Point.h"


class Line
{
private:
	Point start_point, end_point;
	double slope{}, length{};

public:

	double calculate_slope();
	double calculate_length();

	Line(Point, Point);
};
