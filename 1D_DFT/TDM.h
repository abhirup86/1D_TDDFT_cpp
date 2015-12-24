#include "stdafx.h"
#include "mesh.h"
#include "linear_response.h"

template <typename T>
class TDM
{
	grid<T>				*Space;
	unsigned			nmap;
	vector<T*>			map;
	unsigned			indmap;
public:
			TDM(grid<T> &g, const int &map_number);
			~TDM();
	void	Solve(const LinearResponse<T> &LR, const vector<T*> &Basis);
	void	Output();
};

template <typename T>
TDM<T>::TDM(grid<T> &g, const int &map_number) : Space(&g), nmap(map_number) {
	map.resize(nmap);
	for(unsigned i = 0; i < nmap; ++i) map[i] = new T[Space->ngrid*Space->ngrid]();
	indmap = 0;
}

template <typename T>
TDM<T>::~TDM() {
	for(unsigned i = 0; i < nmap; ++i) delete [] map[i];
}

template <typename T>
void TDM<T>::Solve(const LinearResponse<T> &LR, const vector<T*> &Basis) {
	for(;indmap < nmap; ++indmap) {
		for(int j = 0; j < LR.LRMatDim; ++j) {
			T  Xia = LR.mat[indmap*(LR.LRMatDim) + j];
			for(int y = 0; y < Space->ngrid; ++y) {
				for(int x = 0; x < Space->ngrid; ++x) {
					map[indmap][x + y*Space->ngrid] += Xia*Basis[LR.Pairia[j][1]][y]*Basis[LR.Pairia[j][0]][x];
				}
			}
		}
	}
	//take matrix square.
	//for(indmap = 0; indmap < nmap; ++indmap)
	//	for(int i = 0; i < Space->ngrid*Space->ngrid; ++i) map[indmap][i] = map[indmap][i]*map[indmap][i];
}

template <typename T>
void TDM<T>::Output() {
	ofstream f1("TDM.txt");
	for(unsigned i = 0; i < nmap; ++i) {
		for(int y = 0; y < Space->ngrid; ++y) {
			for(int x = 0; x < Space->ngrid; ++x) {
				f1 << map[i][x + y*Space->ngrid] << " ";
			}
			f1 << endl;
		}
		f1 << endl;
	}
	++indmap;
}