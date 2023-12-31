#pragma once
#include"Node.h"
#include"Element.h"
#include"Force.h"
#include"Material.h"
#include"boundary_condition.h"
#include<vector>
#include<fstream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include"Linear_solver.h"


class FEM_component
{
public:
	int no_of_nodes, no_of_elements, no_of_force, no_of_materials, no_of_BC, global_total_DOF, total_ele_dofs;


	std::vector<Node>node_list;
	std::vector<Material>material_list;
	std::vector<Element>ele_list;
	std::vector<Force>force_list;
	std::vector<Boundary_Condition>bc_list;

	Eigen::MatrixXd RHS;
	Eigen::MatrixXd lhs;
	Eigen::MatrixXd displacement_vector{};
	Eigen::MatrixXd reaction_vector{};

	void read_data();
	void read_node_data();
	void read_element_data();
	void read_material_data();
	void read_force_data();
	void read_BC();

	void cal_RHS();
	void cal_stiffness();
	void apply_boundary_con();

	void solve_MUMPS();
	void get_secondary_vars();
	int determine_size_of_Global_A(int symmetry);
	void Create_MUMPS_parameters(int* eltptr, int* eltvar, double* a_elt, double* rhs_);


};


