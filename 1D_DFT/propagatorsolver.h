#include "stdafx.h"
#include "potential.h"
#include "linear_response.h"
#include "system.h"

template <typename T>
class psolver
{
	const LinearResponse<T>		*LR;
	const grid<T>				*Space;
	Hartree<T>					HartreeFix;
	Hartree<T>					HartreeResponse;
	vector<complex<T>* >		Psi_t;
	vector<T* >					Psi_0;
	complex<T>					*CNplusUP;
	complex<T>					*CNminusUP;
	complex<T>					*CNplusDOWN;
	complex<T>					*CNminusDOWN;
	T							dt;
	int							maxdn;
public:
								psolver(const grid<T> &g, const LinearResponse<T> &LinResp, general_system<T> &sys, const int nresponse = 1);
	virtual						~psolver();
	void						Solve();
};

template <typename T>
psolver<T>::psolver(const grid<T> &g, const LinearResponse<T> &LinResp, general_system<T> &sys, const int nresponse) : Space(&g), LR(&LinResp), HartreeResponse(g), HartreeFix(g) {
	Psi_t.resize(sys.Ne/2 + sys.Ne%2);
	Psi_0.resize(sys.Ne/2 + sys.Ne%2);
	for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) {
		Psi_t[i] = new complex<T>[Space->ngrid]();
		Psi_0[i] = new T[Space->ngrid]();
		for(int j = 0; j < Space->ngrid; ++j) {
			Psi_t[i][j] = sys.Psi[i][j];
			Psi_0[i][j] = sys.Psi[i][j];
		}
	}
	dt = 0.02;
	maxdn = nresponse;
	CNplusUP = new complex<T>[Space->ngrid*Space->ngrid]();
	CNplusDOWN = new complex<T>[Space->ngrid*Space->ngrid]();
	CNminusUP = new complex<T>[Space->ngrid*Space->ngrid]();
	CNminusDOWN = new complex<T>[Space->ngrid*Space->ngrid]();
	HartreeFix.Update(sys.density, Space->ngrid);
	for(int i = 0; i < Space->ngrid*Space->ngrid; ++i) {
		if(abs(sys.HamiltonianFix[i]) < 1e-6) continue;
		if(i%Space->ngrid == i/Space->ngrid) continue;
		CNplusUP[i] = complex<T>(0.0, 0.5*dt*sys.HamiltonianFix[i]);
		CNminusUP[i] = complex<T>(0.0, -0.5*dt*sys.HamiltonianFix[i]);
		CNplusDOWN[i] = complex<T>(0.0, 0.5*dt*sys.HamiltonianFix[i]);
		CNminusDOWN[i] = complex<T>(0.0, -0.5*dt*sys.HamiltonianFix[i]);
	}
}

template <typename T>
psolver<T>::~psolver() {
	for(auto i : Psi_t) delete [] i;
	delete [] CNplusUP;
	delete [] CNminusUP;
	delete [] CNplusDOWN;
	delete [] CNminusDOWN;
}

template <typename T>
void psolver<T>::Solve() {
	complex<T> innerpd(0.0, 0.0), innerpd0(0.0, 0.0);
	T norm;
	ofstream f1("HartreeResponse.txt");
	ofstream f2("testpropagator.txt");
	T *DensityAndResponse = new T[Space->ngrid]();
	for(int w = 0; w < 1; ++w) {
		HartreeResponse.Update(LR->dn[w+1], Space->ngrid);
		for(int i = 0; i < Space->ngrid; ++i) HartreeResponse.Data[i] = 2*HartreeResponse.Data[i];
		for(int i = 0; i < Space->ngrid; ++i) {
			f1 << HartreeFix.Data[i] << " " << HartreeResponse.Data[i] << endl;
			CNplusUP[i+i*Space->ngrid] = complex<T>(1.0, 0.5*dt*(HartreeFix.Data[i] - HartreeResponse.Data[i]));
			CNminusUP[i+i*Space->ngrid] = complex<T>(1.0, -0.5*dt*(HartreeFix.Data[i] - HartreeResponse.Data[i]));
			CNplusDOWN[i+i*Space->ngrid] = complex<T>(1.0, 0.5*dt*(HartreeFix.Data[i] + HartreeResponse.Data[i]));
			CNminusDOWN[i+i*Space->ngrid] = complex<T>(1.0, -0.5*dt*(HartreeFix.Data[i] + HartreeResponse.Data[i]));
		}
		f1 << endl;
		for(int it = 0; it < 500; ++it) {
			cout << it << endl;
			complex<T> *Htmp = new complex<T>[Space->ngrid*Space->ngrid]();
			for(unsigned orbital = 0; orbital < Psi_t.size(); ++orbital) {
				memcpy(Htmp, CNplusUP, Space->ngrid*Space->ngrid*sizeof(complex<T>));
				//Psi_t = CNminus*Psi_t;
				LinearAlgebraMul(CNminusUP, Psi_t[orbital], Space->ngrid);
				//call ZGESV(PntNum,NRHS,CNplus,PntNum,IPIV,temp,LDB,info)
				LinearSolverLAPACK(Htmp, Space->ngrid, Psi_t[orbital]);
				//normalization
				norm = 0.0;
				for(int i = 0; i < Space->ngrid; ++i) norm += std::norm(Psi_t[orbital][i]);
				norm = norm*Space->dx;
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = Psi_t[orbital][i]/sqrt(norm);
				innerpd = complex<T>(0.0, 0.0);
				for(int i = 0; i < Space->ngrid; ++i) innerpd += Psi_t[orbital][i]*Psi_0[orbital][i];
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = 2.0*std::real(innerpd)*Psi_0[orbital][i] - Psi_t[orbital][i];
				norm = 0.0;
				for(int i = 0; i < Space->ngrid; ++i) norm += std::norm(Psi_t[orbital][i]);
				norm = norm*Space->dx;
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = Psi_t[orbital][i]/sqrt(norm);

				memcpy(Htmp, CNminusDOWN, Space->ngrid*Space->ngrid*sizeof(complex<T>));
				//Psi_t = CNplus*Psi_t;
				LinearAlgebraMul(CNplusDOWN, Psi_t[orbital], Space->ngrid);
				//call ZGESV(PntNum,NRHS,CNminus,PntNum,IPIV,temp,LDB,info)
				LinearSolverLAPACK(Htmp, Space->ngrid, Psi_t[orbital]);
				//normalization
				norm = 0.0;
				for(int i = 0; i < Space->ngrid; ++i) norm += std::norm(Psi_t[orbital][i]);
				norm = norm*Space->dx;
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = Psi_t[orbital][i]/sqrt(norm);
				innerpd = complex<T>(0.0, 0.0);
				for(int i = 0; i < Space->ngrid; ++i) innerpd += Psi_t[orbital][i]*Psi_0[orbital][i];
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = 2.0*std::real(innerpd)*Psi_0[orbital][i] - Psi_t[orbital][i];
				norm = 0.0;
				for(int i = 0; i < Space->ngrid; ++i) norm += std::norm(Psi_t[orbital][i]);
				norm = norm*Space->dx;
				for(int i = 0; i < Space->ngrid; ++i) Psi_t[orbital][i] = Psi_t[orbital][i]/sqrt(norm);

				//innerpd0 = complex<T>(0.0, 0.0);
				//for(int i = 0; i < Space->ngrid; ++i) innerpd0 += Psi_0[orbital][i]*Psi_0[orbital][i];
				//cout << innerpd << " " << innerpd0 << endl;
				if(orbital != 6) continue;
			}
			delete [] Htmp;
		}
		for(int i = 0; i < Space->ngrid; ++i) {
			for(unsigned orbital = 0; orbital < Psi_t.size(); ++orbital) {
				//f2 << Psi_t[orbital][i] << endl;
				f2 << std::norm(Psi_t[orbital][i]) - Psi_0[orbital][i]*Psi_0[orbital][i] << " ";
			}
			f2 << endl;
		}
	}
}