#include "stdafx.h"
#include "Linear_solver.h"


double * Linear_solver::Evaluate(MUMPS_paramters param)
{
	this->InitializeMUMPSData(param);
	param.rhs = this->MUMPSolution(param);
	this->EndMUMPS();
	return param.rhs;
}

void Linear_solver::InitializeMUMPSData(MUMPS_paramters param)
{
	int ierr, myid;
	int argc; char ** argv;
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	/**** Initialize a MUMPS instance. Use MPI_COMM_WORLD *****/
	/**** Define a F90_communicator for the package****/
	mumps_par.comm_fortran = USE_COMM_WORLD;
	/**** Ask for symetric code****/
	mumps_par.sym = 1;
	//	id.cntl(0);
	/**** Host working****/
	mumps_par.par = 1;
	/**** Initialize instance of the package****/
	mumps_par.job = JOB_INIT;
	dmumps_c(&mumps_par);

	if (myid == 0) {
		mumps_par.n = param.N; // number of equation
		mumps_par.nelt = param.NELT; // element number
		mumps_par.eltptr = param.ELTPTR;
		mumps_par.eltvar = param.ELTVAR;
	}
}

double * Linear_solver::MUMPSolution(MUMPS_paramters param)
{
	int ierr, myid;
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0) {
		mumps_par.a_elt = param.A_ELT;
		mumps_par.rhs = param.rhs;
	}
#define ICNTL(I) icntl[(I)-1] /* macro so that indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro so that indices match documentation */
	/*** Specify element entry ***/
	mumps_par.ICNTL(5) = 1;
	/*** Call package for solution ***/
	mumps_par.CNTL(3) = 1.0E-05;
	mumps_par.CNTL(5) = 1.0E+20;
	mumps_par.ICNTL(1) = -1;
	mumps_par.ICNTL(2) = -1;
	mumps_par.ICNTL(3) = -1;
	mumps_par.ICNTL(13) = 1;
	mumps_par.ICNTL(24) = 1;
	mumps_par.ICNTL(25) = 0;
	//id.ICNTL(4)=2;
	mumps_par.job = JOB_SOLVE;
	//	system("PAUSE");
	dmumps_c(&mumps_par);
	//	for(int m=0;m<eqNum;m++) {printf("%le\n ",id.rhs[m]);} system("PAUSE");;
	return mumps_par.rhs;
}

void Linear_solver::EndMUMPS()
{
	mumps_par.job = JOB_END;
	dmumps_c(&mumps_par);
	MPI_Finalize();
}

Linear_solver::Linear_solver()
{
}


Linear_solver::~Linear_solver()
{
}