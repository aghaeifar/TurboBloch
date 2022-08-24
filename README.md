# Bloch Equation Simulator (C++ implementation)

Efficient and fast implementation of Bloch equation simulator for magnetic resonance imaging (MRI) sequences, supporting parallel transmission (pTx).

## Dependencies:

* Math Kernel Library (MKL)
* Git (http://git-scm.com/)
* Cmake build tool (http://www.cmake.org/)
* Parallelization is implemented based on STL algorithms. Using a compiler with c++17 support is required.


## Linux installation:

Installing dependencies:

```sh
sudo apt-get install g++ cmake git
sudo apt-get install libtbb-dev intel-mkl
```

Clone bloch simulator from repository:

```sh
git clone https://github.com/aghaeifar/bloch_simulator.git
```

Build and install Bloch simulator as a shared library:

```sh
$ cd bloch_simulator
$ mkdir build
$ cd build
$ cmake ..
$ make
$ sudo make install
```

## Windows installation:
MKL can be obtained by installing Intel oneAPI Toolkits ( click [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) to download). I chose online installer and unchecked all the tools except "*Intel OneAPI Math Kernel Library*" and "*Intel OneAPI Threading Building Blocks*" to save a lot of space. More information about configuring oneAPI can be found [here](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-base-linux/top/before-you-begin.html). Briefly, it is needed to set environmental variables by running *<oneAPI_install_dir>/setvars.bat*.

The cmake as stated above should also work in Windows. It will create a Visual Studio sln file. Open and build solution.
## Precompiled binaries
Shared library built for Linux and Windows can be downloaded in the repository releases in **release** mode.

## Macros:

One can define following macros to disable/enable some features in the program:
- ```__SEQUENTIAL__``` disables parallelization and the need for tbb library.
- ```__MEASURE_ELAPSED_TIME__ ``` enables printing elapsed time for calculations.
- ```__SINGLE_PRECISION__``` uses single precision floating-point format. Boost the speed at the cost of precision. All double inputs must be replaced with float.
- ```__EXPORT_CLASS_BLOCH__```, creates exports when building a shared library in Windows.
---

## Short manual:

### MATLAB Interface:

A MATLAB mex wrapper is provided in MATLAB folder with a few examples. Cmake should be able prepare mex file build too.


## Troubleshooting 

Getting this error when running mex file:
```
/sys/os/glnxa64/libstdc++.so.6: version CXXABI_1.3.8' not found
```
**Solution**: update softlink as:
```sh
sudo ln -vfns /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.30 libstdc++.so.6
```


---

## Contributing:

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.