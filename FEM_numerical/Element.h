#pragma once
#include "Node.h"
#include"Material.h"
#include<vector>
#include<Eigen/Dense>
#include"Line.h"


class Element
{
public:
	int id{}, ele_type{};
	double length{};
	std::vector<Node*>connectivity;
	Material* mat;

	Eigen::MatrixXd T, B, K;

	void calculate_stiffness();
	void apply_BC(int n_BC, std::vector<Boundary_Condition> BC_list);
};
