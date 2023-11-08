#include "FEM_component.h"

void FEM_component::read_node_data()
{
	std::ifstream file;
	file.open("Node.txt");

	file >> no_of_nodes;

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
	for (size_t i = 0; i < no_of_elements; i++)
	{
		ele_list[i].apply_BC(no_of_BC, bc_list);
	}
}

