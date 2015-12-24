#include "stdafx.h"
#include "mesh.h"
#include "potential.h"

#ifndef SYSTEM_H
#define SYSTEM_H

template<typename T>
class general_system
{
protected:
	const grid<T>	*Space;
	potential<T>	Vextra;
	bool			BkgrdInitialized;
public:
					general_system(const int &es, const grid<T> &g);
	virtual			~general_system();
	virtual void	Update() = 0;
	virtual T		Kernel(const T &x, const T &xp) const = 0;
	void			CalDensity(const float &mixin);
	void			PrintVextra();
	void			Normalization(bool NormalizeRest = false);
	void			ZeroDiagonal();
	void			SetBackground(const vector<T> &W, const vector<T> &D);
	string			systemname;
	Hartree<T>		VH;
	ExCorr<T>		Vxc;
	vector<T*>		Psi;
	vector<T*>		Basis;
	int				nHamiltonian;
	T				*EigenVectors;
	T				*EigenVectors0;
	T				*EigenValues;
	int				Ne;
	T				*density;
	T				*gsdensity;
	T				*HamiltonianFix;
};

template <typename T>
general_system<T>::general_system(const int &es, const grid<T> &g) : Ne(es), Space(&g), Vextra(g), VH(g), Vxc(g) {
	BkgrdInitialized = false;
	HamiltonianFix = new T[g.ngrid*g.ngrid]();
	EigenVectors = new T[g.ngrid*g.ngrid]();
	EigenVectors0 = new T[g.ngrid*g.ngrid]();
	EigenValues = new T[g.ngrid]();
	nHamiltonian = g.ngrid;
	for(int i = 0; i < Ne/2 + Ne%2; ++i) Psi.push_back(EigenVectors + i*g.ngrid);
	for(int i = 0; i < nHamiltonian; ++i) Basis.push_back(EigenVectors0 + i*g.ngrid);
	density = new T[g.ngrid]();
	gsdensity = new T[g.ngrid]();
	for(int i = 0; i < g.ngrid; ++i) {
		for(int j = i - g.order; j <= i + g.order; ++j) {
			if(j < 0 || j >= g.ngrid) continue;
			HamiltonianFix[j+i*g.ngrid] += -0.5*g.LapCoeff[j+g.order-i]/g.dx/g.dx;
		}
	}
}

template <typename T>
general_system<T>::~general_system() {
	delete [] HamiltonianFix;
	for(unsigned i = 0; i < Psi.size(); ++i) Psi[i] = nullptr;
	delete [] EigenVectors;
	delete [] EigenVectors0;
	delete [] EigenValues;
	delete [] density;
	delete [] gsdensity;
	Space = nullptr;
}

template <typename T>
void general_system<T>::CalDensity(const float &mixin) {
	T* tmp_density = new T[Space->ngrid]();
	for(int i = 0; i < Ne/2; ++i) {
		for(int j = 0; j < Space->ngrid; ++j) {
			tmp_density[j] += 2*std::norm(Psi[i][j]);
		}
	}
	if(Ne%2 == 1) {
		for(int j = 0; j < Space->ngrid; ++j) tmp_density[j] += std::norm(Psi[Ne/2][j]);
	}
	for(int i = 0; i < Space->ngrid; ++i) density[i] = density[i]*(1-mixin) + tmp_density[i]*mixin;
	delete [] tmp_density;
}

template <typename T>
void general_system<T>::ZeroDiagonal() {
	for(int i = 0; i < nHamiltonian; ++i)  {
		int j = i + i*nHamiltonian;
		EigenVectors[j] = 0.0;
	}
}

template <typename T>
void general_system<T>::PrintVextra() {
	ofstream f1("1.txt");
	for(int i = 0; i < Space->ngrid; ++i)  {
		f1 << Vextra.Data[i] << endl;
	}
	f1.close();
}

template <typename T>
void general_system<T>::Normalization(bool NormalizeRest) {
	if(NormalizeRest) {
		for(int i = Ne/2 + Ne%2; i < Space->ngrid; ++i)  {
			for(int j = 0; j < Space->ngrid; ++j) EigenVectors[i*Space->ngrid+j] = EigenVectors[i*Space->ngrid+j]/sqrt(Space->dx);
		}
	} else {
		for(int i = 0; i < Ne/2 + Ne%2; ++i)  {
			for(int j = 0; j < Space->ngrid; ++j) Psi[i][j] = Psi[i][j]/sqrt(Space->dx);
		}
	}
}

template <typename T>
void general_system<T>::SetBackground(const vector<T> &W, const vector<T> &D) {
	if(BkgrdInitialized) return;
	Vextra.Update(W, D);
	for(int i = 0; i < Space->ngrid; ++i) {
		HamiltonianFix[i+i*Space->ngrid] += Vextra.Data[i];
	}
	BkgrdInitialized = true;
	ofstream f1("1.txt");
	for(int i = 0; i < Space->ngrid; ++i)f1 << Vextra.Data[i] << endl;
	f1.close();
}




// non_interacting system
template <typename T>
class non_interacting : public general_system<T>
{
public:
					non_interacting(const int &es, const grid<T> &g);
	virtual void	Update();
	virtual T		Kernel(const T &x, const T &xp) const;
};

template <typename T>
non_interacting<T>::non_interacting(const int &es, const grid<T> &g) : general_system(es, g) {
	systemname = "non-interacting";
	memcpy(EigenVectors, HamiltonianFix, g.ngrid*g.ngrid*sizeof(T));
}

template <typename T>
void non_interacting<T>::Update() {
	memcpy(EigenVectors, HamiltonianFix, nHamiltonian*nHamiltonian*sizeof(T));
}

template <typename T>
T non_interacting<T>::Kernel(const T &x, const T &xp) const {
	return (T) 0.0;
}





// Hartree only
template <typename T>
class only_Hartree : public general_system<T>
{
public:
					only_Hartree(const int &es, const grid<T> &g);
	virtual void	Update();
	virtual T		Kernel(const T &x, const T &xp) const;
};

template <typename T>
only_Hartree<T>::only_Hartree(const int &es, const grid<T> &g) : general_system(es, g) {
	systemname = "HartreeOnly";
	memcpy(EigenVectors, HamiltonianFix, nHamiltonian*nHamiltonian*sizeof(T));
}

template <typename T>
void only_Hartree<T>::Update() {
	VH.Update(density, Space->ngrid);
	memcpy(EigenVectors, HamiltonianFix, nHamiltonian*nHamiltonian*sizeof(T));
	VH.FillHamiltonianMatrix(EigenVectors);
}

template <typename T>
T only_Hartree<T>::Kernel(const T &x, const T &xp) const{
	return VH.Kernel(x, xp);
}




#endif SYSTEM_H