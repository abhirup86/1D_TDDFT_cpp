#include "stdafx.h"
#include "mesh.h"
#include "potential.h"
#include "linearalgebra.h"
#include "debug.h"
#include "PHM.h"

#ifndef PROPAGATOR_H
#define PROPAGATOR_H


template <typename T>
class Propagator
{
	potential<T>    *Field;
	string			Psistartfile;
	ofstream		laserfile;
	ofstream		tddensityfile;
	ofstream		tddensdifffile;
	ofstream		tdKSdensdifffile;
	ofstream		tdKSwffile;
	ofstream		tdKSdensdifffileProjectOut;
	ofstream		tdKSdensdifffileProjectIn;
	bool			newfile;

protected:
	float					dt;
	int						nsteps;
	const grid<T>			*Space;
	T						*TDHamiltonian;
	vector<complex<T>* >	Psi_t;
	complex<T>				*H_Psi;
	vector<T*>				OrbDens0;
	vector<T*>				OrbDensRespInPlane;
	vector<T*>				OrbDensRespOutPlane;
	PHM<T>					InPlanePHM;
	PHM<T>					OutPlanePHM;

public:
					Propagator(const grid<T> &g, potential<T> &TDPotential, const int &ntsteps_ = 1000, const float &dt_ = 0.02);
	virtual			~Propagator();
	virtual void	SolveTDSE() = 0;
	void			Propagate(general_system<T> &sys);
	void			CalTDDensity(general_system<T> &sys);
	void			Output(general_system<T> &sys);
};

template <typename T>
Propagator<T>::Propagator(const grid<T> &g, potential<T> &TDPotential, const int &nsteps_, const float &dt_) : 
	Space(&g), Field(&TDPotential), nsteps(nsteps_), dt(dt_), InPlanePHM(g), OutPlanePHM(g) {
	TDHamiltonian = new T[Space->ngrid*Space->ngrid]();
	H_Psi = new complex<T>[Space->ngrid]();
	InPlanePHM.SetFileName("inplane");
	OutPlanePHM.SetFileName("outplane");
}

template <typename T>
Propagator<T>::~Propagator() {
	TDHamiltonian = nullptr;
	delete [] H_Psi;
	for(vector<complex<T>* >::iterator it = Psi_t.begin(); it < Psi_t.end(); ++it) delete [] *it;
	Space = nullptr;
	Field = nullptr;
}

template <typename T>
void Propagator<T>::Propagate(general_system<T> &sys)
{
	Psi_t.resize(sys.Ne/2 + sys.Ne%2);
	OrbDens0.resize(sys.Ne/2 + sys.Ne%2);
	for(unsigned i = 0; i < Psi_t.size(); ++i) {
		OrbDens0[i] = new T[Space->ngrid]();
	}
	for(int i = 0; i < Space->ngrid; ++i) {
		for(unsigned j = 0; j < Psi_t.size(); ++j){
			OrbDens0[j][i] = 2*sys.Basis[j][i]*sys.Basis[j][i];
		}
	}
	OrbDensRespInPlane.resize((Psi_t.size()));
	for(unsigned i = 0; i < Psi_t.size(); ++i) {
		OrbDensRespInPlane[i] = new T[Space->ngrid]();
	}
	OrbDensRespOutPlane.resize((Psi_t.size()));
	for(unsigned i = 0; i < Psi_t.size(); ++i) {
		OrbDensRespOutPlane[i] = new T[Space->ngrid]();
	}
	if(Psistartfile.empty()) {
		for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) {
			Psi_t[i] = new complex<T>[Space->ngrid]();
			for(int j = 0; j < Space->ngrid; ++j) Psi_t[i][j] = sys.Psi[i][j];
		}
	} else {
		for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) {
			Psi_t[i] = new complex<T>[Space->ngrid]();
		}
		ifstream Pfile(Psistartfile.c_str());
		for(int j = 0; j < Space->ngrid; ++j) {
			int counter(0);
			while (counter < sys.Ne/2 + sys.Ne%2)
			{
				T tmp;
				Pfile >> tmp;
				Psi_t[counter][j].real(tmp);
				//cout<< Psi_t[counter][j] << endl;
				++counter;
			}
			Pfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		Pfile.close();
	}
	float time;
	int iter;
	for(iter = 1; iter < nsteps; ++iter) {
		time = (float)(iter - 0.5)*dt;
		if(iter%100 == 0)cout << iter*dt << endl;
		CalTDDensity(sys);
		sys.Update();
		TDHamiltonian = sys.EigenVectors;
		Field->Update(time);
		Field->FillHamiltonianMatrix(TDHamiltonian);
		SolveTDSE();
		if(iter%50 == 0) {
			Output(sys);
			InPlanePHM.Solve(OrbDensRespInPlane, OrbDens0);
			InPlanePHM.Output(iter);
			OutPlanePHM.Solve(OrbDensRespOutPlane, OrbDens0);
			OutPlanePHM.Output(iter);
		}
		if(iter == 49000) {
			ofstream f01("t49000.txt");
			for(int i = 0; i < Space->ngrid; ++i) {
				for(unsigned j = 0; j < Psi_t.size(); ++j){
					f01 << std::real(Psi_t[j][i]) << " " << std::imag(Psi_t[j][i]) << " ";
				}
				f01 << endl;
			}
			f01.close();
		}
	}
	return;
};

template <typename T>
void Propagator<T>::CalTDDensity(general_system<T> &sys) {
	for(int i = 0; i < Space->ngrid; ++i) sys.density[i] = 0.0;
	for(int i = 0; i < sys.Ne/2; ++i) {
		for(int j = 0; j < Space->ngrid; ++j) {
			sys.density[j] += 2*std::norm(Psi_t[i][j]);
		}
	}
	if(sys.Ne%2 == 1) {
		for(int j = 0; j < Space->ngrid; ++j) sys.density[j] += std::norm(Psi_t[sys.Ne/2][j]);
	}
}

template <typename T>
void Propagator<T>::Output(general_system<T> &sys) {
	if(newfile) {
		laserfile.open(sys.systemname + "tdlaser.txt");
		tddensityfile.open(sys.systemname + "tddensity.txt");
		tddensdifffile.open(sys.systemname + "tddensdiff.txt");
		tdKSdensdifffile.open(sys.systemname + "tdKSdensdiff.txt");
		tdKSdensdifffileProjectOut.open(sys.systemname + "tdKSdensdiffProjectOut.txt");
		tdKSdensdifffileProjectIn.open(sys.systemname + "tdKSdensdiffProjectIn.txt");
		newfile = false;
	}
	for(int i = 0; i < Space->ngrid; ++i) laserfile << Field->Data[i] << endl;
	laserfile << endl;
	for(int i = 0; i < Space->ngrid; ++i) tddensityfile << sys.density[i] << endl;
	tddensityfile << endl;
	for(int i = 0; i < Space->ngrid; ++i) tddensdifffile << sys.density[i] - sys.gsdensity[i] << endl;
	tddensdifffile << endl;
	for(int i = 0; i < Space->ngrid; ++i) {
		for(unsigned j = 0; j < Psi_t.size(); ++j){
			tdKSdensdifffile << 2*std::norm(Psi_t[j][i]) - 2*sys.Basis[j][i]*sys.Basis[j][i] << " ";
		}
		tdKSdensdifffile << endl;
	}
	tdKSdensdifffile << endl;
	////////////////////////////////////
	////compute orbital response
	//vector<T*> WfRes(Psi_t.size());
	//T* orbresp1 = new T[Space->ngrid]();
	//for(unsigned i = 0; i < Psi_t.size(); ++i) {
	//	WfRes[i] = new T[Space->ngrid]();
	//	for(int grd = 0; grd < Space->ngrid; ++grd) {
	//		orbresp1[grd] = std::norm(Psi_t[i][grd])/sys.Basis[i][grd] - sys.Basis[i][grd];
	//	}
	//	memcpy(WfRes[i], orbresp1, Space->ngrid*sizeof(T));
	//}
	////project out
	//for(unsigned j = 0; j < Psi_t.size(); ++j){
	//	for(unsigned pj = 0; pj < Psi_t.size(); ++ pj) {
	//		T prod = std::inner_product(sys.Basis[pj], sys.Basis[pj] + Space->ngrid, WfRes[j], 0.0)*Space->dx;
	//		for(int ix = 0; ix < Space->ngrid; ++ix) WfRes[j][ix] -= prod*sys.Basis[pj][ix];
	//	}
	//}
	////write
	//for(int i = 0; i < Space->ngrid; ++i) {
	//	for(unsigned j = 0; j < Psi_t.size(); ++j){
	//		tdKSdensdifffileProjectOut << WfRes[j][i]*sys.Basis[j][i] << " ";
	//	}
	//	tdKSdensdifffileProjectOut << endl;
	//}
	//tdKSdensdifffileProjectOut << endl;
	////////////////////////////////////
	//project out
	complex<T> C0;
	for(unsigned j = 0; j < Psi_t.size(); ++j){
		complex<T> prod;
		complex<T>* WFOutPlane = new complex<T>[Space->ngrid]();
		memcpy(WFOutPlane, Psi_t[j], Space->ngrid*sizeof(complex<T>));
		for(unsigned pj = 0; pj < Psi_t.size(); ++ pj) {
			prod = complex<T>(0.0, 0.0);
			for(int grd = 0; grd < Space->ngrid; ++grd) {
				prod += sys.Basis[pj][grd]*Psi_t[j][grd];
			}
			prod *=Space->dx;
			if(j == pj) C0 = prod;
			for(int grd = 0; grd < Space->ngrid; ++grd) {
				WFOutPlane[grd] -= prod*sys.Basis[pj][grd];
			}
		}
		for(int grd = 0; grd < Space->ngrid; ++grd) {
			OrbDensRespOutPlane[j][grd] = 2.0*2.0*std::real(std::conj(C0)*WFOutPlane[grd])*sys.Basis[j][grd];//first 2 is doubly occupied
		}
		delete [] WFOutPlane;
	}
	for(int i = 0; i < Space->ngrid; ++i) {
		for(unsigned j = 0; j < Psi_t.size(); ++j){
			OrbDensRespInPlane[j][i] = 2*std::norm(Psi_t[j][i]) - OrbDens0[j][i] - OrbDensRespOutPlane[j][i];
		}
	}
	//write
	for(int i = 0; i < Space->ngrid; ++i) {
		for(unsigned j = 0; j < Psi_t.size(); ++j){
			tdKSdensdifffileProjectOut << OrbDensRespOutPlane[j][i] << " ";
			tdKSdensdifffileProjectIn << OrbDensRespInPlane[j][i] << " ";
		}
		tdKSdensdifffileProjectOut << endl;
		tdKSdensdifffileProjectIn << endl;
	}
	tdKSdensdifffileProjectOut << endl;
	tdKSdensdifffileProjectIn << endl;
	////////////////////////////////////
}




//////////////////////////////
// Crank-Nicholson propagator
template <typename T>
class PropagatorCN : public Propagator<T>
{
	void			CNsetup();
	void			ParamterSetup();
	complex<T>		*CNplus;
	complex<T>		*CNminus;
	bool			CNready;
public:
					PropagatorCN(const grid<T> &g, potential<T> &TDPotential, const int &ntsteps_ = 1000, const float &dt_ = 0.02);
					PropagatorCN(const grid<T> &g, string orbfiletostart, const int &ntsteps_ = 1000, const float &dt_ = 0.02);
	virtual			~PropagatorCN();
	virtual void	SolveTDSE();
};

template <typename T>
PropagatorCN<T>::PropagatorCN(const grid<T> &g, potential<T> &TDPotential, const int &nsteps_, const float &dt_) : Propagator(g, TDPotential, nsteps_, dt_) {
	ParamterSetup();
}

template <typename T>
PropagatorCN<T>::PropagatorCN(const grid<T> &g, string orbfiletostart, const int &nsteps_, const float &dt_) : Propagator(g, orbfiletostart, nsteps_, dt_) {
	ParamterSetup();
}

template <typename T>
PropagatorCN<T>::~PropagatorCN() {
	delete [] CNplus;
	delete [] CNminus;
}

template <typename T>
void PropagatorCN<T>::ParamterSetup() {
	CNplus = new complex<T>[Space->ngrid*Space->ngrid]();
	CNminus = new complex<T>[Space->ngrid*Space->ngrid]();
	CNready = false;
}

template <typename T>
void PropagatorCN<T>::CNsetup() {
	for(int i = 0; i < Space->ngrid*Space->ngrid; ++i) {
		if(abs(TDHamiltonian[i]) < 1e-6) continue;
		if(i%Space->ngrid == i/Space->ngrid) continue;
		CNplus[i] = complex<T>(0.0, 0.5*dt*TDHamiltonian[i]);
		CNminus[i] = complex<T>(0.0, -0.5*dt*TDHamiltonian[i]);
	}
	CNready = true;
}

template <typename T>
void PropagatorCN<T>::SolveTDSE() {
	if(!CNready) CNsetup();
	for(int i = 0; i < Space->ngrid; ++i) {
		CNplus[i+i*Space->ngrid] = complex<T>(1.0, 0.5*dt*TDHamiltonian[i+i*Space->ngrid]);
		CNminus[i+i*Space->ngrid] = complex<T>(1.0, -0.5*dt*TDHamiltonian[i+i*Space->ngrid]);
	}
	complex<T> *Htmp = new complex<T>[Space->ngrid*Space->ngrid]();
	for(unsigned orbital = 0; orbital < Psi_t.size(); ++orbital) {
		memcpy(Htmp, CNplus, Space->ngrid*Space->ngrid*sizeof(complex<T>));
		//Psi_t = CNplus*Psi_t;
		LinearAlgebraMul(CNminus, Psi_t[orbital], Space->ngrid);
		//call ZGESV(PntNum,NRHS,CNminus,PntNum,IPIV,temp,LDB,info)
		LinearSolverLAPACK(Htmp, Space->ngrid, Psi_t[orbital]);
	}
	delete [] Htmp;
}




//////////////////////////////
// Exponential propagator
template <typename T>
class PropagatorExp : public Propagator<T>
{
	int				order;
	complex<T>		*TaylorExp;
	void			ParamterSetup();
public:
					PropagatorExp(const grid<T> &g, potential<T> &TDPotential, const int &ntsteps_ = 1000, const float &dt_ = 0.02);
					PropagatorExp(const grid<T> &g, string orbfiletostart, const int &ntsteps_ = 1000, const float &dt_ = 0.02);
	virtual			~PropagatorExp();
	virtual void	SolveTDSE();
};

template <typename T>
PropagatorExp<T>::PropagatorExp(const grid<T> &g, potential<T> &TDPotential, const int &nsteps_, const float &dt_) : Propagator(g, TDPotential, nsteps_, dt_) {
	ParamterSetup();
}

template <typename T>
PropagatorExp<T>::PropagatorExp(const grid<T> &g, string orbfiletostart, const int &nsteps_, const float &dt_) : Propagator(g, orbfiletostart, nsteps_, dt_) {
	ParamterSetup();
}

template <typename T>
PropagatorExp<T>::~PropagatorExp() {
	delete [] TaylorExp;
}

template <typename T>
void PropagatorExp<T>::ParamterSetup() {
	order = 4;
	TaylorExp = new complex<T>[order + 1]();
	TaylorExp[0] = complex<T>(1.0, 0.0);
	complex<T> a(1.0, 0.0);
	for(int i = 0; i < order; ++i) {
		float tmp = -dt/(i+1);
		a = complex<T>(0.0, 1.0)*a*((T)tmp);
		TaylorExp[i+1] = a;
	}
}

template <typename T>
void PropagatorExp<T>::SolveTDSE() {
	complex<T> *temp = new complex<T>[Space->ngrid]();
	for(unsigned orbital = 0; orbital < Psi_t.size(); ++orbital) {
		memcpy(H_Psi, Psi_t[orbital], Space->ngrid*sizeof(complex<T>));
		for(int i = 1; i <= order; ++i) {
			//H_Psi = Hamiltonian*H_Psi;
			LinearAlgebraMul(TDHamiltonian, H_Psi, Space->ngrid);
			//*it += TaylorExp*H_Psi;
			LinearAlgebraMul(TaylorExp[i], H_Psi, Space->ngrid, temp);
			for(int it = 0; it < Space->ngrid; ++it) Psi_t[orbital][it] += temp[it];
		}
	}
	delete [] temp;
}

#endif //PROPAGATOR_H
