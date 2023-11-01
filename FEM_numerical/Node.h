#pragma once
#include"Point.h"


class Node :
	public Point
{
protected:
	int id{};
public:
	Node(Point);
	void set_ID(int);
	int get_ID();

	Node();
	~Node();
};

