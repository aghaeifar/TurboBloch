// gpu_matrix_mul.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"
#include "gpu_matrix_mul.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <cublas_v2.h>

#include <iostream>

//  mexcuda -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart ...

// TODO: This is an example of a library function
int gpu_matrix_mul::mul_matmat(std::complex<double> *mat_in1, std::complex<double> *mat_in2, std::complex<double> *mat_out, size_t row, size_t col_row, size_t col)
{
	//Allocate space for device copies in device memory
	cuDoubleComplex* cdc_mat_in1;
	cuDoubleComplex* cdc_mat_in2;
	cuDoubleComplex* cdc_mat_out;
	
	cudaMalloc(&cdc_mat_in1, row * col_row * sizeof(cuDoubleComplex));
	cudaMalloc(&cdc_mat_in2, col_row * col * sizeof(cuDoubleComplex));
	cudaMalloc(&cdc_mat_out, row * col * sizeof(cuDoubleComplex));

	if (cudaMemcpy(cdc_mat_in1, mat_in1, row * col_row * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) != cudaSuccess)
		std::cout << "copying from host to device failed" << std::endl;
	if (cudaMemcpy(cdc_mat_in2, mat_in2, col_row * col * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice) != cudaSuccess)
		std::cout << "copying from host to device failed" << std::endl;
	
	
	// https://docs.nvidia.com/cuda/cublas/index.html
	// https://stackoverflow.com/questions/43441573/using-cublas-with-complex-numbers-from-thrust
	cublasHandle_t handle;
	cublasStatus_t stat = cublasCreate(&handle);
	if (stat != CUBLAS_STATUS_SUCCESS) {
		std::cout<<"CUBLAS initialization failed"<<std::endl;
		return EXIT_FAILURE;
	}
	std::complex<double> alpha(1.0, 0.0);
	std::complex<double> beta(0.0, 0.0);

	cuDoubleComplex* _alpha = reinterpret_cast<cuDoubleComplex*>(&alpha);
	cuDoubleComplex* _beta = reinterpret_cast<cuDoubleComplex*>(&beta);

	if(cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, row, col, col_row,
				_alpha, cdc_mat_in1, row, cdc_mat_in2, col_row, _beta, cdc_mat_out, row) != CUBLAS_STATUS_SUCCESS)
		std::cout << "CUBLAS failed to multiply" << std::endl;

	if (cudaMemcpy(mat_out, cdc_mat_out, row * col * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost) != cudaSuccess)
		std::cout << "copying from device to host failed" << std::endl;

	cublasDestroy(handle);
	cudaFree(cdc_mat_in1);
	cudaFree(cdc_mat_in2);
	cudaFree(cdc_mat_out);
    
	cdc_mat_in1 = cdc_mat_in2 = cdc_mat_out = NULL;
    return EXIT_SUCCESS;
}

int gpu_matrix_mul::mul_matmat(std::complex<float> *mat_in1, std::complex<float> *mat_in2, std::complex<float> *mat_out, size_t row, size_t col_row, size_t col)
{
    cuFloatComplex* cdc_mat_in1;
    cuFloatComplex* cdc_mat_in2;
    cuFloatComplex* cdc_mat_out;

    cudaMalloc(&cdc_mat_in1, row * col_row * sizeof(cuFloatComplex));
    cudaMalloc(&cdc_mat_in2, col_row * col * sizeof(cuFloatComplex));
    cudaMalloc(&cdc_mat_out, row * col * sizeof(cuFloatComplex));

    if (cudaMemcpy(cdc_mat_in1, mat_in1, row * col_row * sizeof(cuFloatComplex), cudaMemcpyHostToDevice) != cudaSuccess)
        std::cout << "copying from host to device failed" << std::endl;
    if (cudaMemcpy(cdc_mat_in2, mat_in2, col_row * col * sizeof(cuFloatComplex), cudaMemcpyHostToDevice) != cudaSuccess)
        std::cout << "copying from host to device failed" << std::endl;


    // https://docs.nvidia.com/cuda/cublas/index.html
    // https://stackoverflow.com/questions/43441573/using-cublas-with-complex-numbers-from-thrust
    cublasHandle_t handle;
    cublasStatus_t stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        std::cout<<"CUBLAS initialization failed"<<std::endl;
        return EXIT_FAILURE;
    }
    std::complex<float> alpha(1.0, 0.0);
    std::complex<float> beta(0.0, 0.0);

    cuFloatComplex* _alpha = reinterpret_cast<cuFloatComplex*>(&alpha);
    cuFloatComplex* _beta = reinterpret_cast<cuFloatComplex*>(&beta);

    if(cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, row, col, col_row,
                _alpha, cdc_mat_in1, row, cdc_mat_in2, col_row, _beta, cdc_mat_out, row) != CUBLAS_STATUS_SUCCESS)
        std::cout << "CUBLAS failed to multiply" << std::endl;

    if (cudaMemcpy(mat_out, cdc_mat_out, row * col * sizeof(cuFloatComplex), cudaMemcpyDeviceToHost) != cudaSuccess)
        std::cout << "copying from device to host failed" << std::endl;

    cublasDestroy(handle);
    cudaFree(cdc_mat_in1);
    cudaFree(cdc_mat_in2);
    cudaFree(cdc_mat_out);

    cdc_mat_in1 = cdc_mat_in2 = cdc_mat_out = NULL;
    return EXIT_SUCCESS;
}

gpu_matrix_mul::gpu_matrix_mul()
{

}

gpu_matrix_mul::~gpu_matrix_mul()
{

}
