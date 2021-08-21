

# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 
This work is inspired by Brian Hargreaves Bloch equation simulator [+](http://www-mrsrl.stanford.edu/~brian/blochsim/)

MATLAB mex wrapper is provided and should work with MATLAB R2018a and newer releases. Use command below to compile mex file:

    mex bloch_sim.cpp bloch_sim_mex.cpp -R2018a

I use one of following strategies to make code runs faster:

 Using windows ppl library for multi-threading. This employs `parallel_for` loop to accelerate the calculation which limits the code to be compiled in windows OS [+](https://docs.microsoft.com/en-us/cpp/parallel/concrt/how-to-write-a-parallel-for-loop?view=msvc-160). 

    mex COMPFLAGS='$COMPFLAGS /DUSE_PPL' bloch_sim_mex.cpp bloch_sim.cpp -R2018a

Using OpenMP interface for multi-threading

    mex COMPFLAGS='$COMPFLAGS /openmp' bloch_sim_mex.cpp bloch_sim.cpp -R2018a

Using built-in Intel® MKL optimizations and OpenMP interface, if you have MKL installed

    mex  COMPFLAGS='$COMPFLAGS /DEIGEN_USE_MKL_ALL /openmp /DMKL_ILP64' -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/' ...
         -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core...
         -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
         bloch_sim_mex.cpp bloch_sim.cpp -R2018a

The code relies on the [Eigen](https://eigen.tuxfamily.org) library which is incorporated in this repository, for convenience.
I used MSVC 2019 to compile mex file; however, following command should work if you are using MinGW64:

    mex CXXFLAGS='$CXXFLAGS -std=c++11 -fopenmp' LDFLAGS=-fopenmp bloch_sim_mex.cpp bloch_sim.cpp -R2018a

Let me know if you find out other approaches to accelerate the program.
