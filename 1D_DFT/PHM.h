#include "stdafx.h"
#include "mesh.h"

#ifndef PHM_H
#define PHM_H

template <typename T>
class PHM
{
	const grid<T>		*Space;
	T*					map;
	string				prefix;
public:
			PHM(const grid<T> &g);
			~PHM();
	void	Solve(const vector<T*> &dn, const vector<T*> &n0);	//use when you have dn and n0
	void	Solve2(const vector<T*> &dn, const vector<T*> &n0);	//use when you have delta psi and psi0
	void	Output(const int &idx);
	void	SetFileName(const string &name);
};

template <typename T>
PHM<T>::PHM(const grid<T> &g) : Space(&g) {
	map = new T[Space->ngrid*Space->ngrid]();
}

template<typename T>
void PHM<T>::SetFileName(const string &name){
	prefix = name;
	ofstream f1(prefix + "PHM.txt", std::ofstream::out | std::ofstream::trunc);
	f1.close();
}

template <typename T>
PHM<T>::~PHM() {
	delete [] map;
}

template <typename T>
void PHM<T>::Solve(const vector<T*> &dn, const vector<T*> &n0) {
	for(int i = 0; i < Space->ngrid*Space->ngrid; ++i) map[i] = 0.0;
	for(unsigned i = 0; i < dn.size(); ++i) {
		for(int y = 0; y < Space->ngrid; ++y) {
			for(int x = 0; x < Space->ngrid; ++x) {
				map[x + y*Space->ngrid] += n0[i][x]*dn[i][y];
			}
		}
	}
}

template <typename T>
void PHM<T>::Solve2(const vector<T*> &dPsi, const vector<T*> &Psi0) {
	for(int i = 0; i < Space->ngrid*Space->ngrid; ++i) map[i] = 0.0;
	for(unsigned i = 0; i < dPsi.size(); ++i) {
		for(int y = 0; y < Space->ngrid; ++y) {
			for(int x = 0; x < Space->ngrid; ++x) {
				map[x + y*Space->ngrid] += dPsi[i][y]*Psi0[i][y]*Psi0[i][x]*Psi0[i][x];
			}
		}
	}
}

template <typename T>
void PHM<T>::Output(const int &idx) {
	std::stringstream ss;
	ss << std::setw(7) << std::setfill('0') << idx;
	ofstream f1(prefix + ss.str() + "PHM.txt", std::ofstream::out | std::ofstream::app);
	for(int y = 0; y < Space->ngrid; ++y) {
		for(int x = 0; x < Space->ngrid; ++x) {
			f1 << map[x + y*Space->ngrid] << " ";
		}
		f1 << endl;
	}
	f1 << endl;
	f1.close();
}


#endif //PROPAGATOR_H
