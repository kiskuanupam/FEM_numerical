#include "Node.h"
Node::Node(Point pt)
{
	this->x = pt.get_x();
	this->y = pt.get_y();
}

void Node::set_ID(int)
{
	this->id = id;

}

int Node::get_ID()
{
	return id;
}

Node::Node()
{
}

Node::~Node()
{
}

