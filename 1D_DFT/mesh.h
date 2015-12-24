#include "stdafx.h"
#ifndef MESH_H
#define MESH_H

enum BoundaryType
{
	Wall,
	Periodic,
	Absorption
};


template <typename T>
class grid
{
public:
	grid(const int &n, const T &deltax, const BoundaryType &boundary = Wall, const int &neworder = 1);
	~grid();
	int				ngrid;
	T				dx;
	BoundaryType	bc;
	T				ApplyLap(const T* x);
	int				order;			// finite difference order
	T				*NablaCoeff;	// first order derivative
	T				*LapCoeff;		// Second order derivative
};

template <typename T>
grid<T>::grid(const int &n, const T &deltax, const BoundaryType &boundary = Wall, const int &neworder = 1) : ngrid(n), dx(deltax), bc(boundary), order(neworder) {
	if(ngrid < order + 1) order = ngrid - 1;
	if(ngrid <= 1) throw std::logic_error("grid point number is wrong.");
	if(order > 4) order = 4;
	static T NablaCoeff1[] = {-1.0/2, 0.0, 1.0/2};
	static T NablaCoeff2[] = { 1.0/12, -2.0/3, 0.0, 2.0/3, -1.0/12};
	static T NablaCoeff3[] = {-1.0/60, 3.0/20, -3.0/4, 0.0, 3.0/4, -3.0/20, 1.0/60};
	static T NablaCoeff4[] = { 1.0/280, -4.0/105, 1.0/5, -4.0/5, 0.0, 4/5, -1.0/5, 4.0/105, -1.0/280};
	static T LapCoeff1[] = { 1.0, -2.0, 1.0};
	static T LapCoeff2[] = {-1.0/12, 4.0/3, -5.0/2, 4.0/3, -1.0/12};
	static T LapCoeff3[] = { 1.0/90, -3.0/20, 3.0/2, -49.0/18, 3.0/2, -3.0/20, 1.0/90};
	static T LapCoeff4[] = {-1.0/560, 8.0/315, -1.0/5, 8.0/5, -205.0/72.0, 8.0/5, -1.0/5, 8.0/315, -1.0/560};
	switch (order)
	{
	case(1):
		NablaCoeff = NablaCoeff1;
		LapCoeff = LapCoeff1;
		break;
	case(2):
		NablaCoeff = NablaCoeff2;
		LapCoeff = LapCoeff2;
		break;
	case(3):
		NablaCoeff = NablaCoeff3;
		LapCoeff = LapCoeff3;
		break;
	case(4):
		NablaCoeff = NablaCoeff4;
		LapCoeff = LapCoeff4;
		break;
	default:
		break;
	}
}

template <typename T>
grid<T>::~grid() {
	NablaCoeff = nullptr;
	LapCoeff = nullptr;
}

template <typename T>
T grid<T>::ApplyLap(const T* x) {
	T res(0);
	for(int i = 0; i < 2*order + 1; ++i) res += *LapCoeff * *(x - order + i);
	return res;
}


#endif //MESH_H
