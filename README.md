


# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 

Program needs Intel MKL library and uses Intel Threading Building Blocks (tbb) for multi threading. 

## MATLAB Mex

MATLAB mex wrapper is provided and should work with MATLAB R2018a and newer releases. Use command below to compile mex file:

     mex COMPFLAGS='$COMPFLAGS /DMKL_LP64' -I'C:/Program Files (x86)/Intel/oneAPI/tbb/2021.3.0/include/'...
    -L'C:/Program Files (x86)/Intel/oneAPI/tbb/2021.3.0/lib/intel64/vc14' -ltbb12 ...
    -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/'...
    -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ...
    -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
    bloch_sim_mex.cpp bloch_sim.cpp -R2018a
        

Let me know if you find out other approaches to accelerate the program.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
