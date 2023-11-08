#include "Element.h"

void Element::calculate_stiffness()
{
	// T matrix
	Line ele_line(*connectivity[0], *connectivity[1]);

	length = ele_line.calculate_length();
	double angle = ele_line.calculate_slope();

	Eigen::MatrixXd Temp_T(2 * ele_type, 2 * ele_type);
	Temp_T.setZero();

	for (size_t i = 0; i < 2 * ele_type; i++)
	{
		Temp_T(i, i) = cos(angle);
	}
	Temp_T(0, 1) = Temp_T(2, 3) = -sin(angle);
	Temp_T(1, 0) = Temp_T(3, 2) = sin(angle);

	T = Temp_T;


	//B matrix
	Eigen::MatrixXd Temp_B(1, 2 * ele_type);
	Temp_B.setZero();
	for (size_t i = 0; i < 2 * ele_type; i++)
	{
		Temp_B(0, 0) = -1 / length;
		Temp_B(0, 2) = 1 / length;
	}
	B = Temp_B;


	// K global
	double mult = (mat->area * mat->youngs_mod) * length;

	Eigen::MatrixXd K_local = B.transpose() * B * mult;

	Eigen::MatrixXd K_global = T.transpose() * K_local * T;

	K = K_global;

}

void Element::apply_BC(int n_BC, std::vector<Boundary_Condition> BC_list)
{
	for (size_t i = 0; i < ele_type; i++)
	{
		for (size_t j = 0; j < n_BC; j++)
		{
			if (connectivity[i]->get_ID() == BC_list[j].node_ID)
			{
				int n;
				n = 2 * i + BC_list[j].node_dof;
				K(n, n) += alpha;
			}
		}
	}
}

void Element::cal_strain(Eigen::MatrixXd disp)
{
	strain = B * T * disp;
}

void Element::cal_stress()
{
	stress = strain * (mat->youngs_mod);
}

Eigen::MatrixXd Element::get_reaction()
{
	Eigen::MatrixXd temp;
	temp = T.transpose() * (B.transpose() * stress) * (mat->area * length);
	return temp;
}


