#include "stdafx.h"
#include "system.h"
#include "linearalgebra.h"

template <typename T>
class SCFSolver
{
	int		Maxiter;
	int		Counter;
	double	ConvergeDensity;
	double	ChangeinDensity;
	T		*OldDensity;
	bool	isConverged(const general_system<T> &sys);
	void	CoutConvReport();
public:
			SCFSolver(const int &Max, const general_system<T> &sys, const double &delta_density = -1.0);
	virtual ~SCFSolver();
	void	Solve(general_system<T> &sys);
	void	Load(general_system<T> &sys, const string &filename1, const string &filename2);
	void	Output(general_system<T> &sys);
};

template <typename T>
SCFSolver<T>::SCFSolver(const int &Max, const general_system<T> &sys, const double &delta_density = -1.0) : Maxiter(Max), ConvergeDensity(delta_density) {
	ChangeinDensity = HUGE_VAL;
	OldDensity = new T[sys.nHamiltonian]();
	memcpy(OldDensity, sys.density, sizeof(T)*sys.nHamiltonian);
	ChangeinDensity = sys.Ne;
	Counter = 0;
}

template <typename T>
SCFSolver<T>::~SCFSolver() {
	delete [] OldDensity;
}

template <typename T>
bool SCFSolver<T>::isConverged(const general_system<T> &sys) {
	ChangeinDensity = 0.0;
	for(int i = 0; i < sys.nHamiltonian; ++i)
		ChangeinDensity += abs(sys.density[i] - OldDensity[i]);
	++Counter;
	memcpy(OldDensity, sys.density, sizeof(T)*sys.nHamiltonian);
	return (ChangeinDensity <= ConvergeDensity || Counter >= Maxiter);
}

template <typename T>
void SCFSolver<T>::CoutConvReport() {
	cout << "iter =" << Counter << " < Max =" << Maxiter << "   dens change =" << ChangeinDensity << endl;
}

template <typename T>
void SCFSolver<T>::Solve(general_system<T> &sys) {
	//memcpy(sys.EigenVectors, sys.HamiltonianFix, sys.nHamiltonian*sys.nHamiltonian*sizeof(T));
	EigenSolverLAPACK(sys.EigenVectors, sys.nHamiltonian, sys.EigenValues);
	sys.Normalization();
	sys.CalDensity(1.0);
	do {
		sys.Update();
		EigenSolverLAPACK(sys.EigenVectors, sys.nHamiltonian, sys.EigenValues);
		sys.Normalization();
		sys.CalDensity(0.002f);
		bool Converged = isConverged(sys);
		CoutConvReport();
		if(Converged) break;
	} while(1);
	memcpy(sys.gsdensity, sys.density, sys.nHamiltonian*sizeof(T));
	sys.Normalization(true);
	memcpy(sys.EigenVectors0, sys.EigenVectors, sys.nHamiltonian*sys.nHamiltonian*sizeof(T));
	return;
}

template <typename T>
void SCFSolver<T>::Load(general_system<T> &sys, const string &filename1, const string &filename2) {
	string buff;
	ifstream wfinp(filename1);
	for(int i = 0; i < sys.nHamiltonian; ++i) {
		for(int j = 0; j < sys.nHamiltonian; ++j) {
			getline(wfinp, buff);
			sys.EigenVectors[i*sys.nHamiltonian+j] = stod(buff);
		}
		getline(wfinp, buff);
	}
	sys.CalDensity(1.0);
	wfinp.close();
	ifstream wfenergy(filename2);
	for(int i = 0; i < sys.nHamiltonian; ++i) {
		getline(wfenergy, buff);
		sys.EigenValues[i] = stod(buff);
	}
	wfenergy.close();
	memcpy(sys.gsdensity, sys.density, sys.nHamiltonian*sizeof(T));
	memcpy(sys.EigenVectors0, sys.EigenVectors, sys.nHamiltonian*sys.nHamiltonian*sizeof(T));
	cout << "loaded." << endl;
	return;
}

template <typename T>
void SCFSolver<T>::Output(general_system<T> &sys) {
	ofstream densityfile(sys.systemname + "grounddensity.txt");
	densityfile.precision(std::numeric_limits<T>::digits10);
	densityfile.setf(std::ios::floatfield, std::ios::scientific);
	for(int i = 0; i < sys.nHamiltonian; ++i) densityfile << sys.density[i] << endl;
	densityfile.close();
	ofstream energyfile(sys.systemname + "groundenergy.txt");
	energyfile.precision(std::numeric_limits<T>::digits10);
	energyfile.setf(std::ios::floatfield, std::ios::scientific);
	for(int i = 0; i < sys.nHamiltonian; ++i) energyfile << sys.EigenValues[i] << endl;
	energyfile.close();
	ofstream wffile(sys.systemname + "groundwavefunction.txt");
	wffile.precision(std::numeric_limits<T>::digits10);
	wffile.setf(std::ios::floatfield, std::ios::scientific);
	for(int i = 0; i < sys.nHamiltonian; ++i) {
		for(int j = 0; j < sys.nHamiltonian; ++j) 
			wffile << sys.EigenVectors[i*sys.nHamiltonian + j] << endl;
		wffile << endl;
	}
	wffile.close();
	return;
}
