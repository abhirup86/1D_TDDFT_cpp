#include "stdafx.h"
#include "mesh.h"
#ifndef POTENTIAL_H
#define POTENTIAL_H

// general potential
template <typename T>
class potential
{
public:
					potential(const grid<T> &g);
					potential(const grid<T> &g, const vector<T> &V);
					potential(const grid<T> &g, const T* pV, int n);
	virtual			~potential();
	virtual void	Update(const float &time) {};
	virtual void	Update(const vector<T> &W, const vector<T> &D);
	void			FillHamiltonianMatrix(T *Hamiltonian);
	virtual T		Kernel(const T &x, const T &xp) const;
	T				*Data;
protected:
	const grid<T>	*Space;
};

template <typename T>
potential<T>::potential(const grid<T> &g) : Space(&g) {
	Data = new T[Space->ngrid]();
}

template <typename T>
potential<T>::potential(const grid<T> &g, const vector<T> &V) : Space(&g) {
	if(V.size() != Space.ngrid) throw std::logic_error;
	Data = new T[Space.ngrid]();
	for(int i = 0; i < ngrid; ++i) Data[i] = V[i];
}

template <typename T>
potential<T>::potential(const grid<T> &g, const T* pV, int n) : Space(&g) {
	if(n != Space.ngrid) throw std::logic_error;
	Data = new T[Space.ngrid]();
	for(int i = 0; i < ngrid; ++i) Data[i] = pV[i];
}

template <typename T>
potential<T>::~potential() {
	delete [] Data;
	Space = nullptr;
}

template <typename T>
void potential<T>::Update(const vector<T> &W, const vector<T> &D) {
	for(vector<T>::const_iterator i = W.begin(); i < W.end(); ++i)
		if(*i > 1.0) throw std::logic_error("width must be a number between 0 and 1 (percentage)");
	if(W.size() != D.size()) throw std::logic_error("information mismatch");
	vector<int> Wgrid(W.size());
	T dxhalf = (Space->dx)/2.0;
	T L = Space->dx*Space->ngrid;
	unsigned p1(0), p2(0);
	int i;
	for(i = 0; i < Space->ngrid; ++i) {
		T x = i*Space->dx + dxhalf;
		T ratio = x/L;
		if(ratio > W[p1]) {
			++p1;
			if(p1 >= W.size()) break;
		}
		Data[i] = D[p1];
	}
	if(i < Space->ngrid) {
		for(i; i < Space->ngrid; ++i) Data[i] = D.back();
	}
}

template <typename T>
void potential<T>::FillHamiltonianMatrix(T *Hamiltonian) {
	int j = 0;
	for(T* i = Hamiltonian; i < Hamiltonian + Space->ngrid*Space->ngrid; i += Space->ngrid + 1) {
		*i += Data[j];
		++j;
	}
}

template <typename T>
T potential<T>::Kernel(const T &x, const T &xp) const {
	return (T)0.0;
}



/////////////////////////////////////////////////
// now Hartree
template <typename T>
class Hartree : public potential<T>
{
	T					AlphaConst;
	T					Afactor;
	T					SoftCoulomb(const T &x, const T &xp);
	void				SetPara();
public:
						Hartree(const grid<T> &g);
						Hartree(const grid<T> &g, const T *density, const int &n);
	void				Update(const T *density, const int &n);
	T					Kernel(const T &x, const T &xp) const;
};

template <typename T>
Hartree<T>::Hartree(const grid<T> &g, const T *density, const int &n) : potential(g) {
	SetPara();
	Update(density, n);
}

template <typename T>
Hartree<T>::Hartree(const grid<T> &g) : potential(g) {
	SetPara();
}

template <typename T>
void Hartree<T>::SetPara() {
	AlphaConst = 0.1;
	Afactor = 1.0;
}

template <typename T>
void Hartree<T>::Update(const T *density, const int &n) {
	T* kernel = new T[Space->ngrid]();
	T sum;
	for(int i = 0; i < Space->ngrid; ++i) {
		T x = (i + 0.5)*Space->dx;
		for(int j = 0; j < Space->ngrid; ++j) {
			T xp = (j + 0.5)*Space->dx;
			kernel[j] = Kernel(x, xp);
		}
		sum = 0;
		for(int it = 0; it < n; ++it) sum += kernel[it]*density[it];
		Data[i] = sum*Space->dx;
	}
	delete [] kernel;
}

template <typename T>
T Hartree<T>::Kernel(const T &x, const T &xp) const {
	double diff = x - xp;
	return Afactor/sqrt(diff*diff + AlphaConst*AlphaConst);
}




/////////////////////////////////////////////////
// xc later
template <typename T>
class ExCorr : public potential<T>
{
	string		Type;
public:
				ExCorr(const grid<T> &g);
				ExCorr(const grid<T> &g, const T *density, const int &n);
	void		Update(const T *density, int n);
};

template <typename T>
ExCorr<T>::ExCorr(const grid<T> &g) : potential(g) {
}

template <typename T>
ExCorr<T>::ExCorr(const grid<T> &g, const T *density, const int &n) : potential(g) {
}

template <typename T>
void ExCorr<T>::Update(const T *density, int n) {
}




/////////////////////////////////////////////////
// laser (time-dependent potential)
template <typename T>
class Laser : public potential<T>
{
	T		Amp;
	T		Freq;
	T		Width;
	T		Phase;
	T		*Shape;
public:
			Laser(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase = 0.0);
	void	Update(const float &time);
};

template <typename T>
Laser<T>::Laser(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase) : Amp(Amplitude), Freq(Frequncy), Width(PulseWidth), Phase(PulsePhase), potential(g) {
	Shape = new T[g.ngrid]();
	T	L = g.ngrid*g.dx;
	cout << L << endl;
	for(int i = 0; i < g.ngrid; ++i) {
		Shape[i] = -L/2.0 + (i+0.5)*g.dx;
	}
}

template <typename T>
void Laser<T>::Update(const float &time) {
	double Amp_d = Amp;
	double Freq_d = Freq;
	double Width_d = Width;
	double Phase_d = Phase;
	double Factor_d = Amp_d*pow(sin(3.14159265359/Width_d*time),2)*sin(Freq_d*time+Phase_d);
	T	Factor = Factor_d;
	if(time > Width_d) Factor = 0.0;
	for(int i = 0; i < Space->ngrid; ++i) Data[i] = Factor*Shape[i];
	return;
}




/////////////////////////////////////////////////
// laser style 1 (time-dependent potential)
template <typename T>
class LaserS1 : public potential<T>
{
	T		Amp;
	T		Freq;
	T		Width;
	T		Phase;
	T		*Shape;
public:
			LaserS1(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase = 0.0);
	void	Update(const float &time);
};

template <typename T>
LaserS1<T>::LaserS1(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase) : Amp(Amplitude), Freq(Frequncy), Width(PulseWidth), Phase(PulsePhase), potential(g) {
	Shape = new T[g.ngrid]();
	for(int i = 0; i < 20; ++i) {
		Shape[i] = 10.0;
		Shape[g.ngrid - i - 1] = -10.0;
	}
}

template <typename T>
void LaserS1<T>::Update(const float &time) {
	double Amp_d = Amp;
	double Freq_d = Freq;
	double Width_d = Width;
	double Phase_d = Phase;
	double Factor_d = Amp_d*pow(sin(3.14159265359/Width_d*time),2)*sin(Freq_d*time+Phase_d);
	T	Factor = Factor_d;
	if(time > Width_d) Factor = 0.0;
	for(int i = 0; i < Space->ngrid; ++i) Data[i] = Factor*Shape[i];
	return;
}




/////////////////////////////////////////////////
// laser style 2 (time-dependent potential)
template <typename T>
class LaserS2 : public potential<T>
{
	T		Amp;
	T		Freq;
	T		Width;
	T		Phase;
	T		*Shape;
public:
			LaserS2(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase = 0.0);
	void	Update(const float &time);
};

template <typename T>
LaserS2<T>::LaserS2(const grid<T> &g, const T &Amplitude, const T &Frequncy, const T &PulseWidth, const T &PulsePhase) : Amp(Amplitude), Freq(Frequncy), Width(PulseWidth), Phase(PulsePhase), potential(g) {
	Shape = new T[g.ngrid]();
	T	L = g.ngrid*g.dx;
	cout << L << endl;
	for(int i = 0; i < g.ngrid; ++i) {
		Shape[i] = 10.0*sin(-L/2.0 + (i+0.5)*g.dx);
	}
}

template <typename T>
void LaserS2<T>::Update(const float &time) {
	double Amp_d = Amp;
	double Freq_d = Freq;
	double Width_d = Width;
	double Phase_d = Phase;
	double Factor_d = Amp_d*pow(sin(3.14159265359/Width_d*time),2)*sin(Freq_d*time+Phase_d);
	T	Factor = Factor_d;
	if(time > Width_d) Factor = 0.0;
	for(int i = 0; i < Space->ngrid; ++i) Data[i] = Factor*Shape[i];
	return;
}


#endif //POTENTIAL_H