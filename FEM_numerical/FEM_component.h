#pragma once
#include"Node.h"
#include"Element.h"
#include"Force.h"
#include"Material.h"
#include"boundary_condition.h"
#include<vector>
#include<fstream>
#include<Eigen/Dense>


class FEM_component
{
public:
	int no_of_nodes, no_of_elements, no_of_force, no_of_materials, no_of_BC;


	std::vector<Node>node_list;
	std::vector<Material>material_list;
	std::vector<Element>ele_list;
	std::vector<Force>force_list;
	std::vector<Boundary_Condition>bc_list;
	Eigen::MatrixXd RHS;

	void read_node_data();
	void read_element_data();
	void read_material_data();
	void read_force_data();
	void read_BC();
	void cal_RHS();


};


