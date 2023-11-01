// FEM_numerical.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include<iostream>
#include"Point.h"
int main()
{
	Point p1, p2;
	p1.set_xy(0, 0);
	p2.set_xy(3, 4);
	std::cout << p1.distance_from_point(p2);
}
