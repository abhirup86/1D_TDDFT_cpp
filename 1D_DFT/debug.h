#include "stdafx.h"

#ifndef DEBUG_H
#define DEBUG_H

ofstream dfile("debug.txt");

template <typename T>
void PrintMatrix(const T *mat, const int &dim, ofstream &file, const bool &noTrans = true) {
	if(noTrans) {
		for(int i = 0; i < dim; ++i) {
			for(int j = 0; j < dim; ++j)
				file << std::setw(20) << std::setprecision(12) << mat[j + i*dim] << " ";
			file << endl;
		}
		file << endl;
	} else {
		for(int i = 0; i < dim; ++i) {
			for(int j = 0; j < dim; ++j)
				file << std::setw(20) << std::setprecision(12) << mat[i + j*dim] << " ";
			file << endl;
		}
		file << endl;
	}
}

template <typename T>
void PrintVector(const T *vec, const int &dim, ofstream &file) {
	for(int i = 0; i < dim; ++i) {
		file << std::setw(20) << std::setprecision(12) << vec[i] << endl;
	}
	file << endl;
}


#endif //DEBUG_H