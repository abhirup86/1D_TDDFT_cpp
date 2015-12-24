#include "stdafx.h"
#include "linear_response.h"
#include "mesh.h"
#include "potential.h"
#include "system.h"

#ifndef GREENSFUNCTION_H
#define GREENSFUNCTION_H

template <typename T>
class GreensFunction
{
	void						SetMat(const general_system<T> &sys, const int &freqind);
	const LinearResponse<T>		*LR;
	const grid<T>				*Space;
	T							*KernelMat;
	int							maxorb;
	Hartree<T>					HartreeResponse;
	ofstream					f2;
public:
								GreensFunction(const grid<T> &g, const LinearResponse<T> &LinResp, const int maxlvl = 1);
	virtual						~GreensFunction();
	void						CalOrbResponse(const general_system<T> &sys);
	int							maxresponselvl;
	vector<vector<T*>>			OrbResponse;
};

template <typename T>
GreensFunction<T>::GreensFunction(const grid<T> &g, const LinearResponse<T> &LinResp, const int maxlvl) : Space(&g), LR(&LinResp), HartreeResponse(g) {
	maxorb = -1;
	for(unsigned i = 0; i < LR->Pairia.size(); ++i) maxorb = std::max(maxorb, LR->Pairia[i][1]);
	++maxorb;
	KernelMat = new T[maxorb*maxorb]();
	maxresponselvl = maxlvl;
	f2.open("22.txt");
}

template <typename T>
GreensFunction<T>::~GreensFunction() {
	LR = nullptr;
	delete [] KernelMat;
	for(int w = 0; w < maxresponselvl; ++w) {
		for(unsigned i = 0; i < OrbResponse[0].size(); ++i) delete [] OrbResponse[w][i];
	}
}

template <typename T>
void GreensFunction<T>::SetMat(const general_system<T> &sys, const int &freqind) {
	HartreeResponse.Update(LR->dn[freqind], Space->ngrid);
	for(int i = 0; i < Space->ngrid; ++i) {
		HartreeResponse.Data[i] += 0.002*(i - Space->ngrid/2 + 0.5)*Space->dx;
		f2 << (i + 0.5)*Space->dx << " " << LR->dn[freqind][i] << " " << HartreeResponse.Data[i] << endl;
	}
	f2 << endl;
	for(int i = 0; i < maxorb; ++i) {
		for(int j = i; j < maxorb; ++j) {
			T tmp(0);
			for(int g = 0; g < Space->ngrid; ++g) tmp += sys.Basis[i][g]*sys.Basis[j][g]*HartreeResponse.Data[g];
			KernelMat[j+i*maxorb] = Space->dx*tmp;
			if(i != j)KernelMat[i+j*maxorb] = KernelMat[j+i*maxorb];
		}
	}
}

template <typename T>
void GreensFunction<T>::CalOrbResponse(const general_system<T> &sys) {
	OrbResponse.resize(maxresponselvl);
	for(int w = 0; w < maxresponselvl; ++w) OrbResponse[w].resize(sys.Ne/2 + sys.Ne%2);
	for(int w = 0; w < maxresponselvl; ++w) {
		for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) OrbResponse[w][i] = new T[Space->ngrid]();
	}
	ofstream ff1("orbitalresponse.txt");
	ofstream ff2("orbitaldensityresponse.txt");
	for(int w = 0; w < maxresponselvl; ++w) {
		SetMat(sys, w);
		for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) {
			//for(int j = sys.Ne/2 + sys.Ne%2; j < maxorb; ++j) {
			for(int j = 0; j < maxorb; ++j) {
				T factor;
				if(LR->lrtheory == "TammDancoff")
					factor = 1.0/(LR->Energy[w] - sys.EigenValues[j] + sys.EigenValues[i]);
				else if(LR->lrtheory == "FullCasida")
					factor = 2.0*(sys.EigenValues[j] - sys.EigenValues[i])/(LR->Energy[w]*LR->Energy[w] - pow(sys.EigenValues[j] - sys.EigenValues[i],2));
				//if(i == j) continue;
				for(int g = 0; g < Space->ngrid; ++g) OrbResponse[w][i][g] += factor*KernelMat[i+j*maxorb]*sys.Basis[j][g];
			}
		}
		for(int g = 0; g < Space->ngrid; ++g)  {
			for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) ff1 << OrbResponse[w][i][g] << " ";
			for(int i = 0; i < sys.Ne/2 + sys.Ne%2; ++i) ff2 << OrbResponse[w][i][g]*sys.Basis[i][g] << " ";
			ff1 << endl;
			ff2 << endl;
		}
		ff1 << endl;
		ff2 << endl;
	}
}



#endif //GREENSFUNCTION_H