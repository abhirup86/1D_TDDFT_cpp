#include "stdafx.h"
#include "mesh.h"
#include "potential.h"
#include "system.h"
#include "debug.h"

//FIRE mixing is not finished yet
template <typename T>
class SIESTA_fireparameters
{
public:
	int n;
	int step;
	int npos;
	T alpha;
	T dt;
	int nmin;
	T dt_inc;
	T dt_dec;
	T alphamax;
	T alpha_factor;
	T dtmax;
	T *v;
			SIESTA_fireparameters(const int &n_, const T &dt_, const int &nmin_ = 5, const T &dt_inc_ = 1.1, const T &dt_dec_ = 0.5, const T &alphamax_ = 0.1, const T &alpha_factor_ = 0.99);
			~SIESTA_fireparameters() {delete [] v;};
	void	GetPara(T const * const f, T * const x, T const * const maxstep, const bool &normalize);
};

template <typename T>
SIESTA_fireparameters<T>::SIESTA_fireparameters(const int &n_, const T &dt_, const int &nmin_ = 5, const T &dt_inc_ = 1.1, const T &dt_dec_ = 0.5, const T &alphamax_ = 0.1, const T &alpha_factor_ = 0.99) :
	n(n_), dt(dt_), nmin(nmin_), dt_inc(dt_inc_), dt_dec(dt_dec_), alphamax(alphamax_), alpha_factor(alpha_factor_), step(0), npos(0), dtmax(10.0), v(nullptr) {
    alpha = alphamax;
	v = new double[n]();
}

// For a general problem, we should choose dt and the masses
// so that the first step does not travel too far compared to some
// "typical length". Say, delta_x = 0.1 bohr
template <typename T>
void SIESTA_fireparameters<T>::GetPara(T const * const f, T * const x, T const * const maxstep, const bool &normalize) {
	T norm_v, norm_f;
	T *dx = new T[n]();
    int *count = new int[n]();
	T dt1(dt);
	++step;
    T power = std::inner_product(f, f+n, v, 0);
	if (power < 0) {
		memset(v, 0, n*sizeof(T));
		dt *= dt_dec;
		alpha = alphamax;
		npos = 0;
	} else {
		norm_v = std::inner_product(v, v+n, v, 0);
		norm_v = sqrt(norm_v);
		norm_f = std::inner_product(f, f+n, f, 0);
		norm_f = sqrt(norm_f);
		for(int j = 0; j < n; ++j)
			v[j] = (1.0 - alpha)*v[j] + alpha*norm_v*f[j]/norm_f;
		if (npos > nmin) {
			dt = std::min(dt*dt_inc, dtmax);
			alpha = alpha*alpha_factor;
		}
		++npos;
	}
	// Euler with midpoint correction
    for(int j = 0; j < n; ++j) dx[j] = v[j]*dt1/2.0;		//half of the mid-point rule (v(0))
    for(int j = 0; j < n; ++j) v[j] = v[j] + f[j]*dt1;		//v at dt
    for(int j = 0; j < n; ++j) dx[j] = dx[j] + v[j]*dt1/2.0;//the other half  (v(t))
	if (normalize) {
		int counter(0);
		for(int j = 0; j < n; ++j) {
			if(abs(dx[j]) > maxstep[j]) {
				++counter;
				break;
			}
		}
		if (counter != 0) {
			T norm_x = sqrt(std::inner_product(dx, dx+n, dx, 0));
			for(int j = 0; j < n; ++j) {
				dx[j] = maxstep[j] * dx[j] / norm_x;
				x[j] += dx[j];
			}
		}
	} else {
		for(int j = 0; j < n; ++j) {
			if(abs(dx[j]) <= maxstep[j]) x[j] += dx[j];
			else {
				if(dx[j] > 0) x[j] += maxstep[j];
				else x[j] += -maxstep[j];
				++count;
			}
		}
	}
    delete [] dx;
	delete [] count;
}



template <typename T>
class inversion
{
	grid<T> *Space;
	T		para_a;
	T		para_b;
	SIESTA_fireparameters<T> *MixMethod;
	void	mixpotential(T * const potential, T const * const newpotential);
public:
			inversion(grid<T> &g);
	void	invert(general_system<T> &sys, T const * const density);
};

template <typename T>
inversion<T>::inversion(grid<T> &g) : Space(&g),  para_a(4.0), para_b(0.0) {
	MixMethod = new SIESTA_fireparameters<T>(Space->ngrid, 0.1);
}

template <typename T>
void inversion<T>::mixpotential(T * const potential, T const * const newpotential) {
	T *diffpotential = new T[Space->ngrid]();
	for(int i = 0; i < Space->ngrid; ++i) diffpotential[i] = newpotential[i] - potential[i];
	T dmax = *std::max_element(diffpotential, diffpotential + Space->ngrid);
	T *dmaxvec = new T[Space->ngrid]();
	for(int i = 0; i < Space->ngrid; ++i) dmaxvec[i] = dmax;
	MixMethod->GetPara(diffpotential, potential, dmaxvec, true);
	delete [] diffpotential;
	delete [] dmaxvec;
}

template <typename T>
void inversion<T>::invert(general_system<T> &sys, T const * const density) {
	potential<T> myV(*Space);
	T *myVnew = new T[Space->ngrid]();
	for(int i = 0; i < Space->ngrid; ++i) myV.Data[i] = 100.0;
	T error(0);
	for(int i = 0; i < Space->ngrid; ++i) error += abs(density[i] - sys.density[i]);
	int counter(0);
	while(error > 1e-9) {
		for(int i = 0; i < Space->ngrid; ++i) myVnew[i] = (sys.density[i]+para_a)/(density[i]+para_a);
		//mixpotential(myVnew, myV.Data);
		for(int i = 0; i < Space->ngrid; ++i) myV.Data[i] = myV.Data[i]*(0.999+0.001*(sys.density[i]+para_a)/(density[i]+para_a));
		sys.Update();
		sys.ZeroDiagonal();
		myV.FillHamiltonianMatrix(sys.EigenVectors);
		EigenSolverLAPACK(sys.EigenVectors, sys.nHamiltonian, sys.EigenValues);
		sys.Normalization();
		sys.CalDensity(1.0);
		error = 0.0;
		for(int i = 0; i < Space->ngrid; ++i) error += abs(density[i] - sys.density[i]);
		cout << "\r" << error;
		//cout << error << endl;
		++counter;
		if(counter > 0) break; //after 2000 iterations we cut it anyway.
	}
	cout << endl;
	ofstream invertVfile(sys.systemname + "invertV.txt");
	PrintVector(myV.Data, Space->ngrid, invertVfile);
	ofstream invertKSdensfile(sys.systemname + "invertKSdens.txt");
	for(int g = 0; g < Space->ngrid; ++g) {
		for(unsigned i = 0; i < sys.Psi.size(); ++i) invertKSdensfile << pow(sys.EigenVectors[g+i*Space->ngrid], 2) - pow(sys.EigenVectors0[g+i*Space->ngrid], 2) << " ";
		invertKSdensfile << endl;
	}
	ofstream invertwffile(sys.systemname + "invertwf.txt");
	PrintMatrix(sys.EigenVectors, sys.nHamiltonian, invertwffile, false);
		    /*do
       error=0.D0;Vmax=0.D0;Vmin=100000000.D0
       call SolveSE_single_particle(V_ks_Guess,tmp_Evalues,tmp_Evectors)
       call ground_Dens_recon(tmp_Evectors,Dens0)
       do i=1,GridNum
          V_ks(i)=V_ks_Guess(i)*(Dens0(i)+para_a)/(Dens(i)+para_a)
          error_tmp=dabs(Dens0(i)-Dens(i))
          if (error<error_tmp) error=error_tmp
          if (Vmax<dabs(V_ks(i))) Vmax=V_ks(i)           
          if (Vmin>dabs(V_ks(i)))Vmin=V_ks(i)
       enddo
       if (error<=tolerence) then
          write(unit=*,fmt='(A45I10A6A1)'),"static calculation done...  iteration goes ",imax," times",char(13)
          exit
       endif
       imax=imax+1;V_ks_Guess=V_ks
    enddo*/
}
