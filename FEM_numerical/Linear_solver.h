#pragma once
#include <mpi.h>
#include "dmumps_c.h"
# include<fstream>
#include <string.h>
# define JOB_INIT			-1
# define JOB_END			-2
# define JOB_SOLVE			6
# define USE_COMM_WORLD		-987654
struct MUMPS_paramters
{
	int		N{};
	int		NELT{};
	int		NVAR{};
	int*	ELTPTR = nullptr;
	int*	ELTVAR = nullptr;
	double*	A_ELT = nullptr;
	double* rhs = nullptr;
};
class Linear_solver
{
private:
	DMUMPS_STRUC_C		mumps_par; //mumps_par is an object of type DMUMPS_STRUC_C structure in dmumps_c.h
public:
	double*				Evaluate(MUMPS_paramters);
	void				InitializeMUMPSData(MUMPS_paramters);
	double*				MUMPSolution(MUMPS_paramters);
	void				EndMUMPS();

	Linear_solver();
	~Linear_solver();
};