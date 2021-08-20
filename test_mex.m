clc
mex -v COMPFLAGS='$COMPFLAGS /openmp' bloch_sim_mex.cpp bloch_sim.cpp -R2018a
% 
%%
clc
mex  COMPFLAGS='$COMPFLAGS /DEIGEN_USE_MKL_ALL /openmp /DMKL_ILP64' -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/' ...
     -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core...
     -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
     bloch_sim_mex.cpp bloch_sim.cpp -R2018a

% COPTIMFLAGS="-O3 -DNDEBUG"

%% test speed
gamma = 267522187.44;
flipangle = 20 * pi/180;
ntime = 800; % number of samples
dur = 2; % second
npos = 100000; % number of spatial positions

b1 = complex(ones(ntime,1) * flipangle/gamma/dur);
tp = ones(ntime,1) * dur/ntime;

gr = zeros(ntime, 3);
b0 = zeros(npos, 1);
pr = rand(npos, 3);

tic
bloch_sim_mex(b1, gr, tp, b0, pr);
toc

%% compare runtime
ntime = 800; % number of samples
npos = 100000; 
e_sens = rand(npos,1) + 1j+rand(npos,1);
e_b1 = rand(ntime,1)+ 1j+rand(ntime,1);
e_gr = rand(ntime, 3);
e_pr = rand(npos, 3);
e_b0 = rand(npos, 1);
e_tp = rand(ntime, 1);
GAMMA_T = 25656565;
e_m0 = rand(npos, 3);

ge_sens = gpuArray(e_sens);
ge_b1 = gpuArray(e_b1);


tic
e_b1comb = e_b1 * e_sens';
toc

tic 
    rotz = (e_gr * e_pr' + e_b0') .* e_tp * -1.0 * GAMMA_T;
    rotx = real(e_b1comb) .* e_tp * -1.0 * GAMMA_T;
    roty = imag(e_b1comb) .* e_tp * GAMMA_T; 
    e_m0t = e_m0';
toc

tic
fe_b1comb = ge_b1 * ge_sens';
toc
