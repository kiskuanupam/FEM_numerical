#pragma once
class Point
{
protected:
	double x{}, y{};
public:
	void set_xy(double, double);
	double get_x();
	double get_y();

	double distance_from_point(Point);
};

