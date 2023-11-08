#include "FEM_component.h"

void FEM_component::read_node_data()
{
	std::ifstream file;
	file.open("Node.txt");

	file >> no_of_nodes;
	global_total_DOF = 2 * no_of_nodes;

	int id{};
	double x, y;

	Node temp;

	for (size_t i = 0; i < no_of_nodes; i++)
	{
		file >> id >> x >> y;

		temp.set_ID(id);
		temp.set_xy(x, y);

		node_list.push_back(temp);
	}
	file.close();
}

void FEM_component::read_element_data()
{
	std::ifstream file;
	file.open("Element.txt");

	file >> no_of_elements;

	int id, n, mat_type;

	for (size_t i = 0; i < no_of_elements; i++) //no of element loop
	{
		Element temp;
		file >> id ;
		temp.id = id;
		temp.ele_type = 2;

		for (size_t j = 0; j < temp.ele_type; j++) //element type
		{
			file >> n;
			for (size_t k = 0; k < no_of_nodes; k++)  //no of nodes
			{
				if (node_list[k].get_ID() == n)
				{
					temp.connectivity.push_back(&node_list[k]);
				}
				break;
			}
		}

		file >> mat_type;
		for (size_t l = 0; l < no_of_materials; l++)
		{
			if (mat_type == material_list[i].id)
			{
				temp.mat = &material_list[l];
			}
			break;
		}
		ele_list.push_back(temp);
	}
	file.close();

}

void FEM_component::read_material_data()
{
	std::ifstream file;
	file.open("Material_Property.txt");

	int id;
	double Youngs_mod, area;

	file >> no_of_materials;

	Material temp;

	for (size_t i = 0; i < no_of_materials; i++)
	{
		file >> id >> Youngs_mod >> area;
		temp.id = id;
		temp.youngs_mod = Youngs_mod;
		temp.area = area;

		material_list.push_back(temp);

	}
	file.close();


}

void FEM_component::read_force_data()
{
	std::ifstream file;
	file.open("Force.txt");

	file >> no_of_force;

	int id, dof_id;
	double force_value;

	Force temp;

	for (size_t i = 0; i < no_of_force; i++)
	{
		file >> id >> dof_id >> force_value;
		temp.id = id;
		temp.dof_id = dof_id;
		temp.force_value = force_value;

		force_list.push_back(temp);

	}
	file.close();
}

void FEM_component::read_BC()
{

	std::ifstream file;
	file.open("Boundary_Condition.txt");

	file >> no_of_BC;

	int id, dof_id;
	double bc_value;

	Boundary_Condition temp;

	for (size_t i = 0; i < no_of_BC; i++)
	{
		file >> id >> dof_id >> bc_value;
		temp.node_ID = id;
		temp.node_dof = dof_id;
		temp.bc_value = bc_value;

		bc_list.push_back(temp);
	}
	file.close();
}

void FEM_component::cal_RHS()
{
	Eigen::MatrixXd temp(2 * no_of_elements, 1);
	temp.setZero();

	for (size_t i = 0; i < no_of_force; i++)
	{
		for (size_t j = 0; j < no_of_nodes; j++)
		{
			if (force_list[i].id == node_list[j].get_ID())
			{
				temp(2 * node_list[j].get_ID() + force_list[i].dof_id, 0) = force_list[i].force_value;
			}
		}
	}
	RHS = temp;
}

void FEM_component::cal_stiffness()
{
	for (size_t i = 0; i < no_of_elements; i++)
	{
		ele_list[i].calculate_stiffness();
	}
}

void FEM_component::apply_boundary_con()
{
	std::cout << "" << std::endl;
	std::cout << "Elemental k matrix after BC-" << std::endl;
	for (size_t i = 0; i < no_of_elements; i++)
	{
		ele_list[i].apply_BC(no_of_BC, bc_list);
	}
}
/*


*/
int FEM_component::determine_size_of_Global_A(int symmetry)
{
	int glCnt = 0;
	for (int i = 0; i < this->no_of_elements; i++) {
		int n = 4; // no of dof of each element
		if (symmetry == 1) glCnt += n * (n + 1) / 2;
		if (symmetry == 0) glCnt += n * n;
	}
	return glCnt;
}

void FEM_component::Create_MUMPS_parameters(int* eltptr, int* eltvar, double* a_elt, double* rhs_)
{
	long sum = 0;
	for (int i = 0; i < this->no_of_elements; i++) {
		int nn = 4;
		for (int i_ = 0; i_ < nn; i_++) {
			for (int j_ = i_; j_ < nn; j_++) {
				a_elt[sum++] = this->ele_list[i].K(j_, i_);
			}
		}
	}
	sum = 0;
	eltptr[0] = 1;
	for (int i = 0; i < this->no_of_elements; i++) {
		if (i > 0) eltptr[i] = eltptr[i - 1] + 4;
		for (int m = 0; m < 2; m++) {
			eltvar[sum++] = 2 * ele_list[i].connectivity[m]->get_ID() + 1;
			eltvar[sum++] = 2 * ele_list[i].connectivity[m]->get_ID() + 2;
		}
	}
	eltptr[this->no_of_elements] = eltptr[this->no_of_elements - 1] + 4;

	for (int i = 0; i < this->global_total_DOF; i++) {
		rhs_[i] = this->RHS(i, 0);
	}

}

void FEM_component::solve_MUMPS()
{
	int N{};
	int NELT{};
	int NVAR{};
	int NVAL{};
	int* ELTPTR = nullptr;
	int* ELTVAR = nullptr;
	double* A_ELT = nullptr;
	double* RHS = nullptr;

	N = this->global_total_DOF; // total no of variable
	NELT = this->no_of_elements;
	NVAR = total_ele_dofs; // no of dof for each element

	ELTPTR = new int[NELT + 1];
	ELTVAR = new int[NVAR];
	RHS = new double[N];
	NVAL = this->determine_size_of_Global_A(1);
	A_ELT = new double[NVAL];

	this->Create_MUMPS_parameters(ELTPTR, ELTVAR, A_ELT, RHS);
	MUMPS_paramters param;
	param.N = N;
	param.NELT = NELT;
	param.NVAR = NVAR;
	param.ELTPTR = ELTPTR;
	param.ELTVAR = ELTVAR;
	param.A_ELT = A_ELT;
	param.rhs = RHS;

	Linear_solver solve;
	param.rhs = solve.Evaluate(param);

	Matrix<double> temp(this->global_total_DOF, 1);
	this->displacement_vector = temp;

	for (int i = 0; i < this->global_total_DOF; i++) this->displacement_vector(i, 0) = param.rhs[i];

	std::cout << "Printing displacement............\n";
	std::cout << this->displacement_vector << std::endl;
	std::cin.get();

	delete[] param.ELTPTR;
	delete[] param.ELTVAR;
	delete[] param.A_ELT;
	delete[] param.rhs;

}

void FEM_component::get_secondary_vars()
{
	Matrix<double> ele_disp(4, 1);
	Matrix<double> rec_vect(global_total_DOF, 1);
	for (int i = 0; i < this->no_of_elements; i++) {
		int nn = ele_list[i].ele_type;
		for (int j = 0; j < nn; j++) { // 2 nodes
			for (int m = 0; m < 2; m++) { // 2 dofs
				ele_disp(2 * j + m, 0) = displacement_vector(2 * ele_list[i].connectivity[j]->get_ID() + m, 0);
			}
		}
		ele_list[i].cal_strain(ele_disp);
		ele_list[i].cal_stress();
		Matrix<double> ele_rec = ele_list[i].get_reaction();
		for (int j = 0; j < nn; j++) { // 2 nodes
			for (int m = 0; m < 2; m++) { // 2 dofs
				rec_vect(2 * ele_list[i].connectivity[j]->get_ID() + m, 0) += ele_rec(2 * j + m, 0);
			}
		}
	}
	reaction_vector = rec_vect;
	std::cout << "Printing Reaction Vector............\n";
	std::cout << reaction_vector << std::endl;
	std::cin.get();
}

