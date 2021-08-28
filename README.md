


# Bloch Equation Simulator (C++ implementation + GPU)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 
This work is inspired by Brian Hargreaves Bloch equation simulator [+](http://www-mrsrl.stanford.edu/~brian/blochsim/)

**Enabling GPU:**
 - install the most recent driver for your NVIDIA adapter 
 - install CUDA Toolkit : https://developer.nvidia.com/cuda-downloads

Run some CUDA samples and ensure they work fine. I tried this which uses cuBLAS:

https://github.com/NVIDIA/cuda-samples/tree/master/Samples/simpleCUBLAS

Please note some CUDA versions are compatible with specific MS Visual Studio. I tried cuda_11.4 with Visual Studio 2019.

**MATLAB Mex**
MATLAB mex wrapper is provided and should work with MATLAB R2018a and newer releases. Use command below to compile mex file:

     mexcuda -DUSE_GPU -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart bloch_sim_mex.cpp bloch_sim.cpp ./gpu_matrix_mul/gpu_matrix_mul.cu -R2018a
        

Program uses windows ppl library for multi-threading. This employs `parallel_for` loop to accelerate the calculation which limits the code to be compiled in windows OS [+](https://docs.microsoft.com/en-us/cpp/parallel/concrt/how-to-write-a-parallel-for-loop?view=msvc-160). 

Let me know if you find out other approaches to accelerate the program.

