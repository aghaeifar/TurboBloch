



# Bloch Equation Simulator (C++ implementation)
Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx). 

## Installing MKL
Intel Math Kernel Library (MKL), is a library of math routines optimized for science and engineering. Bloch simulator uses MKL for matrix multiplication and Intel Threading Building Blocks (tbb) for multi threading. You can download Intel oneAPI (collection of useful tools including MKL and tbb) from [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html). I chose online installer and unchecked all the tools except "*Intel OneAPI Math Kernel Library*" and "*Intel oneAPI Threading Building Blocks*" to save space.

## Compiling
we need first to run *setvars* script to set environment variables for use with the oneAPI. I used following lines in command-line to compile under Linux (Ubuntu)

    source /opt/intel/oneapi/setvars.sh 
    g++ bloch_sim.cpp ./CPU/bloch.cpp -o libbloch_sim.so \
        -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
        ${MKLROOT}/lib/intel64/libmkl_tbb_thread.a \
        ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group \
        -L${TBBROOT}/lib/intel64/gcc4.8 -ltbb -lstdc++ -lpthread -lm -ldl \
        -m64  -I"${MKLROOT}/include" 

I tried these commands in Windows command-line to compile a mex file, to be used in MATLAB:

    "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
    mex -I"%MKLROOT%\include" mkl_intel_lp64.lib mkl_tbb_thread.lib mkl_core.lib tbb12.lib bloch_sim.cpp CPU/bloch.cpp -R2018a
If you would like to compile the program directly in MATLAB command window, merge both commands and use *system()* function:

    system('"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" & mex -I"%MKLROOT%\include" mkl_intel_lp64.lib mkl_tbb_thread.lib mkl_core.lib tbb12.lib bloch_sim.cpp CPU/bloch.cpp -R2018a');


## MATLAB Example

Please see the example code which executes a selective excitation using a since pulse. Pulse is designed for an RF coil with 8Tx (here all are homogeneous). Phase of pulse is designed to be 45 degree.


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

