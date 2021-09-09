#pragma once

#include <complex>
#include <cublas_v2.h>
// forward declaration 

class gpu_matrix_mul
{
public:
	gpu_matrix_mul();
	~gpu_matrix_mul();
    int mul_matmat(std::complex<double>* mat_in1, std::complex<double>* mat_in2, std::complex<double>* mat_out, size_t row = 0, size_t col_row = 0, size_t col = 0);
    int mul_matmat(std::complex<float>* mat_in1, std::complex<float>* mat_in2, std::complex<float>* mat_out, size_t row = 0, size_t col_row = 0, size_t col = 0);
	   
private:
};
