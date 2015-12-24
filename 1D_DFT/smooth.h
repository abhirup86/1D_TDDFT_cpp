#include "stdafx.h"
#include "mesh.h"

template<typename T>
class smooth
{
public:
	smooth(const grid<T> &g);
	~smooth();
	void				LoadFunction(const string &filename);
	void				Variance();		//order vs variance
private:
	int					maxorder;
	const grid<T>*		Space;
	T*					gridx;
	T*					y;
	T*					y_reg;
	T*					var;
};

template<typename T>
smooth<T>::smooth(const grid<T> &g) : Space(&g) {
	maxorder = 20;
	gridx = new T[Space->ngrid+1]();
	y = new T[Space->ngrid+1]();
	y_reg = new T[Space->ngrid+1]();
	var = new T[maxorder]();
	for(int i = 0; i < Space->ngrid; ++i) {
		gridx[i] = (i+0.5)*Space->dx;
	}
}

template<typename T>
smooth<T>::~smooth() {
	delete [] gridx;
	delete [] y;
	delete [] y_reg;
	delete [] var;
}

template<typename T>
void smooth<T>::LoadFunction(const string &filename) {
	string buff;
	ifstream f1(filename.c_str());
	for(int i = 0; i < Space->ngrid; ++i) {
		getline(f1, buff);
		y[i] = stod(buff);
	}
}

template<typename T>
void smooth<T>::Variance() {
	ofstream order_var("order_var.txt");
	ofstream order_reg("order_reg.txt");
	for(int order = 1; order <= maxorder; ++order) {
		T* x = new T[(order + 1)*Space->ngrid]();
		//fill x
		for(int i = 0; i < order; ++i) {
			for(int grd = 0; grd < Space->ngrid; ++grd) {
				x[grd + Space->ngrid*i] = pow(gridx[grd], i+1);
			}
		}
		for(int grd = 0; grd < Space->ngrid; ++grd) x[grd + Space->ngrid*order] = 1;
		int oi, oj;
		//compute x^T*x
		T* xTx = new T[(order + 1)*(order + 1)]();
		for(oi = 0; oi <= order; ++oi) {
			for(oj = oi; oj <= order; ++oj) {
				xTx[oi + oj*(order + 1)] = std::inner_product(x + oi*Space->ngrid, x + oi*Space->ngrid+ Space->ngrid, x + oj*Space->ngrid, 0.0);
				if(oi != oj) xTx[oj + oi*(order + 1)] = xTx[oi + oj*(order + 1)];
			}
		}
		//compute x^T*y
		T* xTy = new T[order + 1]();
		for(oi = 0; oi <= order; ++oi) {
			xTy[oi] = std::inner_product(x + oi*Space->ngrid, x + (oi + 1)*Space->ngrid, y, 0.0);
		}
		//solve coefficients and store in xTy
		LinearSolverLAPACK(xTx, order + 1, xTy);
		for(int i = 0; i <= order; ++i) cout << xTy[i] << " ";
		cout << endl;
		cout << endl;
		//compute y by linear regression
		for(int i = 0; i < Space->ngrid; ++i) y_reg[i] = 0.0;
		for(int i = 0; i < Space->ngrid; ++i) {
			for(int j = 0; j <= order; ++j) {
				y_reg[i] += xTy[j]*x[i + j*Space->ngrid];
			}
		}
		for(int i = 0; i < Space->ngrid; ++i) order_reg << y_reg[i] << endl;
		order_reg << endl;
		//compute the variance.
		T* diff = new T[Space->ngrid]();
		for(int i = 0; i < Space->ngrid; ++i) diff[i] = y[i] - y_reg[i];
		var[order - 1] += std::inner_product(diff, diff + Space->ngrid, diff, 0.0);
		order_var << order << " " << var[order - 1] << endl;
	}
}
