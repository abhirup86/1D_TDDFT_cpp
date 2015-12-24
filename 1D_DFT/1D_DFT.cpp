// 1D_DFT.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mesh.h"
#include "system.h"
#include "potential.h"
#include "scfsolver.h"
#include "propagator.h"
#include "linear_response.h"
#include "GreensFunction.h"
#include "propagatorsolver.h"
#include "PHM.h"
#include "TDM.h"
#include "smooth.h"


int _tmain(int argc, _TCHAR* argv[])
{
	grid<double> box(200, 0.05);
	//smooth<double> myfun(box);
	//myfun.LoadFunction("222.txt");
	//myfun.Variance();
	//return 0;
	//trimple wells
	//only_Hartree<double> sys(8, box);
	//double w[] = {  0.200,   0.450,   0.550,  0.800,   1.000};
	//double d[] = {-14.000,   0.000,  -4.000,  0.000, -12.000};
	//vector<double> vw(w, w+5);
	//vector<double> vd(d, d+5);
	//1D molecule
	only_Hartree<double> sys(14, box);
	double w[] = {0.050,  0.200,  0.300,  0.375,  0.400,  0.475,  0.625,  0.900,  1.000};
	double d[] = {0.000, -6.000,  0.000,  2.000,  0.000, -2.000,  0.000, -3.000,  0.000};
	vector<double> vw(w, w+9);
	vector<double> vd(d, d+9);
	sys.SetBackground(vw, vd);
	sys.PrintVextra();
	cout << "here" << endl;
	SCFSolver<double> mySCFSolver1(10000, sys, 1e-12);
	//mySCFSolver1.Solve(sys);
	//mySCFSolver1.Output(sys);
	mySCFSolver1.Load(sys, "HartreeOnlygroundwavefunction.txt", "HartreeOnlygroundenergy.txt");
	Hartree<double> KernelHolder(box);
	//LinearResponseTammDancoff<double> lr(box, 50);
	LinearResponseCasida<double> lr(box, 50);


	//lr.Solve(sys, KernelHolder);
	//lr.Output(sys);
	//exit(0);



	//stationary
	//GreensFunction<double> mygf(box, lr, 20);
	//mygf.CalOrbResponse(sys);
	//psolver<double> mypsolver(box, lr, sys, 10);
	//mypsolver.Solve();
	//PHM<double> myphm(box, mygf.maxresponselvl);
	//for(int i = 0; i < mygf.maxresponselvl; ++i) myphm.Solve2(mygf.OrbResponse[i], sys.Basis);
	//myphm.Output();
	//TDM<double> mytdm(box, mygf.maxresponselvl);
	//for(int i = 0; i < mygf.maxresponselvl; ++i) mytdm.Solve(lr, sys.Basis);
	//mytdm.Output();
	//time evolution
	//2nd excitation
	Laser<double> mylaser(box, 0.001, 1.34087, 140.577);
	PropagatorCN<double> mypropagator(box, mylaser, 50000, 0.004f);
	//3rd excitation
	//LaserS2<double> mylaser(box, 0.005, 1.81353, 103.935);
	//PropagatorCN<double> mypropagator(box, mylaser, 70000, 0.002f);
	//Laser<double> mylaser(box, 0.0000003, 0.0143941, 7000.0);
	mypropagator.Propagate(sys);
	return 0;
}

