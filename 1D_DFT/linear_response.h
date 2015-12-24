#include "stdafx.h"
#include "mesh.h"
#include "system.h"
#include "linearalgebra.h"

#ifndef LINEAR_RESPONSE_H
#define LINEAR_RESPONSE_H

template <typename T>
class LinearResponse
{
	void					MakePair(const general_system<T> &sys);
protected:
	const grid<T>			*Space;
	int						Maxunocc;
	T						*ediff;
public:
							LinearResponse(const grid<T> &g, const int &unocc = -1);
	virtual void			BuildMatrix(const general_system<T> &sys, const potential<T> &V) = 0;
	virtual void			PostProcess() = 0;
	void					Solve(const general_system<T> &sys, const potential<T> &V);
	void					Output(const general_system<T> &sys);
	string					lrtheory;
	vector<vector<int>>		Pairia;
	T						*mat;
	int						LRMatDim;
	vector<T*>				dn;
	T						*Energy;
};

template <typename T>
LinearResponse<T>::LinearResponse(const grid<T> &g, const int &unocc) : Space(&g) {
	if(unocc < 0) Maxunocc = -1;
	else Maxunocc = unocc;
}

template <typename T>
void LinearResponse<T>::MakePair(const general_system<T> &sys) {
	int nocc = sys.Ne/2 + sys.Ne%2;
	if(Maxunocc < nocc || Maxunocc > sys.nHamiltonian) Maxunocc = sys.nHamiltonian;
	int nall = Maxunocc;
	for(int i = 0; i < nocc; ++i) {
		for(int a = nocc; a < nall; ++a) {
			int pair[] = {i, a};
			Pairia.push_back(vector<int>(pair, pair+2));
			dn.push_back(new T[Space->ngrid]());
		}
	}
}

template <typename T>
void LinearResponse<T>::Solve(const general_system<T> &sys, const potential<T> &V) {
	MakePair(sys);
	LRMatDim = Pairia.size();
	mat = new T[LRMatDim*LRMatDim]();
	Energy = new T[LRMatDim]();
	ediff = new T[LRMatDim]();
	for(int i = 0; i < LRMatDim; ++i) ediff[i] = sys.EigenValues[Pairia[i][1]] - sys.EigenValues[Pairia[i][0]];
	BuildMatrix(sys, V);
	EigenSolverLAPACK(mat, LRMatDim, Energy);
	PostProcess();
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = 0; j < LRMatDim; ++j) {
			T  Xia = mat[i*LRMatDim + j];
			for(int g = 0; g < sys.nHamiltonian; ++g) dn[i][g] += Xia*sys.Basis[Pairia[j][1]][g]*sys.Basis[Pairia[j][0]][g];
		}
	}
}

template <typename T>
void LinearResponse<T>::Output(const general_system<T> &sys) {
	ofstream lrenergyfile(sys.systemname + lrtheory+ "lrenergy.txt");
	std::sort(ediff, ediff+LRMatDim);
	for(int i = 0; i < LRMatDim; ++i) lrenergyfile << ediff[i] << " " << Energy[i] << endl;
	ofstream lrdensityfile(sys.systemname + lrtheory+ "lrdensdiff.txt");
	for(int i = 0; i < LRMatDim; ++i) {
		for(int g = 0; g < sys.nHamiltonian; ++g)
			lrdensityfile << dn[i][g] << endl;
		lrdensityfile << endl;
	}
	ofstream lrcoeffsfile(sys.systemname + lrtheory+ "lrcoeffs.txt");
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = 0; j < LRMatDim; ++j)
			lrcoeffsfile << Pairia[j][0]+1 << " " << Pairia[j][1]+1 << " " << mat[i*LRMatDim + j] << endl;
		lrcoeffsfile << endl;
	}
	ofstream lrKSdensdiff1file(sys.systemname + lrtheory+ "lrKSdensdiff1.txt");
	vector<vector<T*> > tmp(LRMatDim);
	for(int j = 0; j < LRMatDim; ++j) {
		for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i)
			tmp[j].push_back(new T[sys.nHamiltonian]());
	}
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = 0; j < LRMatDim; ++j) {
			T  Xia = mat[i*LRMatDim + j];
			for(int g = 0; g < sys.nHamiltonian; ++g) tmp[i][Pairia[j][0]][g] += Xia*sys.Basis[Pairia[j][1]][g]*sys.Basis[Pairia[j][0]][g];
		}
	}
	for(int i = 0; i < LRMatDim; ++i) {
		for(int g = 0; g < sys.nHamiltonian; ++g) {
			for(int j = 0; j < sys.Ne/2 + sys.Ne%2; ++j) lrKSdensdiff1file << tmp[i][j][g] << " ";
			lrKSdensdiff1file << endl;
		}
		lrKSdensdiff1file << endl;
	}
}





template <typename T>
class LinearResponseTammDancoff : public LinearResponse<T>
{
public:
			LinearResponseTammDancoff(const grid<T> &g, const int &unocc = -1);
	void	BuildMatrix(const general_system<T> &sys, const potential<T> &V);
	void	PostProcess() {};
};

template <typename T>
LinearResponseTammDancoff<T>::LinearResponseTammDancoff(const grid<T> &g, const int &unocc) : LinearResponse(g, unocc) {
	lrtheory = "TammDancoff";
}

template <typename T>
void LinearResponseTammDancoff<T>::BuildMatrix(const general_system<T> &sys, const potential<T> &V) {
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = i; j < LRMatDim; ++j) {
			T integral(0.0);
			//cout << Pairia[j][0] << " " << Pairia[j][1] << endl;
			for(int gi = 0; gi < Space->ngrid; ++gi) {
				T x = (gi + 0.5)*Space->dx;
				for(int gj = 0; gj < Space->ngrid; ++gj) {
					T xp = (gj + 0.5)*Space->dx;
					integral += sys.Basis[Pairia[j][0]][gj]*sys.Basis[Pairia[j][1]][gj]*(V.Kernel(x, xp))*sys.Basis[Pairia[i][0]][gi]*sys.Basis[Pairia[i][1]][gi];
				}
			}
			integral = integral*Space->dx*Space->dx;
			//occupation is dropped since I think it will be handled by the equation itself.
			mat[i + j*LRMatDim] = integral;
			mat[j + i*LRMatDim] = integral;
		}
	}
	for(int i = 0; i < LRMatDim; ++i) mat[i + i*LRMatDim] += sys.EigenValues[Pairia[i][1]] - sys.EigenValues[Pairia[i][0]];
}





//"Time-Dependent Density-Functional Theory Concepts and Applications" Carsten A. Ullrich eq. 7.138~7.141
template <typename T>
class LinearResponseCasida : public LinearResponse<T>
{
	const general_system<T> *psys;
	const potential<T> *pV;
public:
			LinearResponseCasida(const grid<T> &g, const int &unocc);
	void	BuildMatrix(const general_system<T> &sys, const potential<T> &V);
	void	PostProcess();
};

template <typename T>
LinearResponseCasida<T>::LinearResponseCasida(const grid<T> &g, const int &unocc) : LinearResponse(g, unocc) {
	lrtheory = "FullCasida";
}

template <typename T>
void LinearResponseCasida<T>::BuildMatrix(const general_system<T> &sys, const potential<T> &V) {
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = i; j < LRMatDim; ++j) {
			T integral(0.0);
			for(int gi = 0; gi < Space->ngrid; ++gi) {
				T x = (gi + 0.5)*Space->dx;
				for(int gj = 0; gj < Space->ngrid; ++gj) {
					T xp = (gj + 0.5)*Space->dx;
					integral += sys.Basis[Pairia[j][0]][gj]*sys.Basis[Pairia[j][1]][gj]*(V.Kernel(x, xp))*sys.Basis[Pairia[i][0]][gi]*sys.Basis[Pairia[i][1]][gi];
				}
			}
			integral = integral*Space->dx*Space->dx;
			//occupation is dropped since I think it will be handled by the equation itself.
			mat[i + j*LRMatDim] = 2.0*integral*sqrt(ediff[i]*ediff[j]);
			mat[j + i*LRMatDim] = mat[i + j*LRMatDim];
		}
	}
	for(int i = 0; i < LRMatDim; ++i) mat[i + i*LRMatDim] += ediff[i]*ediff[i];
	psys = &sys;
	pV = &V;
}

template <typename T>
void LinearResponseCasida<T>::PostProcess() {
	for(int i = 0; i < LRMatDim; ++i) Energy[i] = sqrt(Energy[i]);
	//solve X-Y from Z
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = 0; j < LRMatDim; ++j) {
			mat[i*LRMatDim + j] = mat[i*LRMatDim + j]/sqrt(ediff[j]);
		}
	}
	//solve X+Y since (A+B)(X+Y) = -Omega*(X-Y)
	T *ApB = new T[LRMatDim*LRMatDim]();
	T *ApBomega = new T[LRMatDim*LRMatDim]();
	for(int i = 0; i < LRMatDim; ++i) {
		for(int j = i; j < LRMatDim; ++j) {
			T integral(0.0);
			for(int gi = 0; gi < Space->ngrid; ++gi) {
				T x = (gi + 0.5)*Space->dx;
				for(int gj = 0; gj < Space->ngrid; ++gj) {
					T xp = (gj + 0.5)*Space->dx;
					integral += psys->Basis[Pairia[j][0]][gj]*psys->Basis[Pairia[j][1]][gj]*(pV->Kernel(x, xp))*psys->Basis[Pairia[i][0]][gi]*psys->Basis[Pairia[i][1]][gi];
				}
			}
			integral = integral*Space->dx*Space->dx;
			ApB[i + j*LRMatDim] = 2.0*integral;
			ApB[j + i*LRMatDim] = 2.0*integral;
		}
	}
	for(int i = 0; i < LRMatDim; ++i) ApB[i + i*LRMatDim] += ediff[i];
	for(int omega = 0; omega < LRMatDim; ++omega) {
		for(int i = 0; i < LRMatDim*LRMatDim; ++i) ApBomega[i] = -1.0/Energy[omega]*ApB[i];
		LinearSolverLAPACK(ApBomega, LRMatDim, mat+omega*LRMatDim);
		cout << "omega " << omega << endl;
	}
	delete [] ApB;
	delete [] ApBomega;
	//normalization
	for(int i = 0; i < LRMatDim; ++i) {
		T Norm = 0;
		for(int j = 0; j < LRMatDim; ++j) {
			Norm += mat[i*LRMatDim + j]*mat[i*LRMatDim + j];
		}
		Norm = sqrt(Norm);
		for(int j = 0; j < LRMatDim; ++j) {
			mat[i*LRMatDim + j] = mat[i*LRMatDim + j]/Norm;
		}
	}
}





#endif //LINEAR_RESPONSE_H
