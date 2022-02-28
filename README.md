



# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 

Program can use parallelism to speedup the simulation. The available methods are:

 - *parallel_for()* based on Parallel Patterns Library (PPL); this works if operating system is **Windows**.
 - *oneAPI Threading Building Blocks (oneTBB)*; a multi-platform library to add parallelism to applications. The latest release can be downloaded from its GitHub repository ([here](https://github.com/oneapi-src/oneTBB)). To use this method, add predefined macro *_TBB* to compiler.
 - *openMP*; it is another multi-platform method. A proper flag should be set for compiler. I found it is slower than other two methods in the current implementation. 

Sequential run will be used if none of above is available.
Program either can be compiled as shared library or be included and directly used in other applications. 

## Compiling
### Installing MKL
Intel Math Kernel Library (MKL), is a library of math routines optimized for science and engineering. Bloch simulator uses MKL for a matrix multiplication. You can download Intel oneAPI (collection of useful tools including MKL) from intel ([here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)). I chose online installer and unchecked all the tools except "*Intel OneAPI Math Kernel Library*"  to save space. The program only needs following libraries: *mkl_intel_lp64, mkl_intel_thread, mkl_core,* and *libiomp5md*.

### Create MATLAB mex
An interface is programmed to compile the simulator as mex file and use it in MATLAB.
I tried these commands in MATLAB to compile the program. Example is provided in *MATLAB_example* folder.

**using parallel_for()**

    mex -I'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\intel64' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win' ...
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -llibiomp5md ...
        MATLAB_example/bloch_mex.cpp bloch.cpp -R2018a

**using openMP**

    mex -I'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\intel64' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win' ...
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -llibiomp5md ...
        COMPFLAGS="$COMPFLAGS /openmp" ...
        MATLAB_example/bloch_mex.cpp bloch.cpp -R2018a

**using TBB library**

    mex -I'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\mkl\latest\lib\intel64' ...
        -L'C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\compiler\lib\intel64_win' ...
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -llibiomp5md ...
        -I'./oneapi-tbb/include' ...
        -L'./oneapi-tbb/lib/intel64/vc14' ...
        -ltbb -D_TBB ...    
        MATLAB_example/bloch_mex.cpp bloch.cpp -R2018a


## MATLAB Example

Please see the example code which executes a selective excitation using a sinc pulse. Pulse is designed for an RF coil with 8Tx (here all are homogeneous). Phase of pulse is designed to be 45 degree.


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

