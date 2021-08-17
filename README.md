
# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 
This work is inspired by Brian Hargreaves Bloch equation simulator [+](http://www-mrsrl.stanford.edu/~brian/blochsim/)

MATLAB mex wrapper is provided and should work with MATLAB R2018a and newer releases. Use command below to compile mex file:

    mex bloch_sim.cpp bloch_sim_mex.cpp -R2018a

The code relies on the [Eigen](https://eigen.tuxfamily.org) library which is incorporated in this repository, for convenience.

I am using `parallel_for` loop to accelerate the calculation which limits the code to be compiled in windows OS [+](https://docs.microsoft.com/en-us/cpp/parallel/concrt/how-to-write-a-parallel-for-loop?view=msvc-160). One should be able to replace `parallel_for` with ordinary `for` loop if you wish to compile in other OS.