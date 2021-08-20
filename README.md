

# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 
This work is inspired by Brian Hargreaves Bloch equation simulator [+](http://www-mrsrl.stanford.edu/~brian/blochsim/)

MATLAB mex wrapper is provided and should work with MATLAB R2018a and newer releases. Use command below to compile mex file:

    mex bloch_sim.cpp bloch_sim_mex.cpp -R2018a

OpenMP interface allowed me run the code slightly faster. I used following command to compile mex

    mex COMPFLAGS='$COMPFLAGS /openmp' bloch_sim_mex.cpp bloch_sim.cpp -R2018a

The code relies on the [Eigen](https://eigen.tuxfamily.org) library which is incorporated in this repository, for convenience.
Eigen can benefit from built-in Intel® MKL optimizations, if you have it installed. I used following command to compile:

    mex  COMPFLAGS='$COMPFLAGS /DEIGEN_USE_MKL_ALL /openmp /DMKL_ILP64' -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/' ...
         -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core...
         -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
         bloch_sim_mex.cpp bloch_sim.cpp -R2018a

I am also using `parallel_for` loop to accelerate the calculation which limits the code to be compiled in windows OS [+](https://docs.microsoft.com/en-us/cpp/parallel/concrt/how-to-write-a-parallel-for-loop?view=msvc-160). One should be able to replace `parallel_for` with ordinary `for` loop if you wish to compile in other OS.
