%% compile with Microsoft Visual C++ 2019 
clc
mex bloch_sim_mex.cpp bloch_sim.cpp -R2018a

%% with ppl
clc
mex COMPFLAGS='$COMPFLAGS /DUSE_PPL' bloch_sim_mex.cpp bloch_sim.cpp -R2018a

%% with openblas
clc
mex  COMPFLAGS='$COMPFLAGS /DUSE_PPL /DEIGEN_USE_BLAS' -I'./openblas0.3.17/include/' ...
     -L'./openblas0.3.17/lib' -llibopenblas bloch_sim_mex.cpp bloch_sim.cpp -R2018a

 %% with openmp and mkl 
% DEIGEN_USE_MKL_ALL DUSE_PPL openmp
% include: C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/
% lib: C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/ & C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win
clc

mex  COMPFLAGS='$COMPFLAGS /DEIGEN_USE_MKL_ALL /DUSE_PPL /openmp /DMKL_LP64' -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/' ...
     -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core...
     -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
     bloch_sim_mex.cpp bloch_sim.cpp -R2018a
 
 %% with GPU - needs to install updated graphic driver
 % download CUDA from https://developer.nvidia.com/cuda-downloads
 % include: C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\include
 % lib: C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64
clc
cd('D:\OneDrive - University College London\Matlab\B1GUI\ipopt_lib\bloch_simulator')
mex  COMPFLAGS='$COMPFLAGS /DUSE_PPL /DUSE_GPU' -I'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\include' ...
     -L'./gpu_matrix_mul/x64/Release/' -lgpu_matrix_mul -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart_static ...
     bloch_sim_mex.cpp bloch_sim.cpp -R2018a
 
%% compile with MinGW64 Compiler (C++) and openmp
% see https://uk.mathworks.com/matlabcentral/answers/279171
clc
mex CXXFLAGS='$CXXFLAGS -std=c++11 -fno-math-errno -ffast-math -march=native -fopenmp' LDFLAGS=-fopenmp bloch_sim_mex.cpp bloch_sim.cpp -R2018a
%%
clc
cd('D:\OneDrive - University College London\Matlab\B1GUI\ipopt_lib\bloch_simulator')
% NVCCFLAGS='$NVCCFLAGS
mexcuda -DUSE_GPU -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart ...
       bloch_sim_mex.cpp bloch_sim.cpp ./gpu_matrix_mul/gpu_matrix_mul.cu -R2018a
   
%% test speed

gamma = 267522187.44;
flipangle = 90 * pi/180;
ntime = 375; % number of samples
dur = 2; % second
npos = 28000; % number of spatial positions
ncoil = 8;

b1 = complex(ones(ntime,ncoil) * flipangle/gamma/dur/ncoil) ;
sens = complex(ones(ncoil,npos));
tp = dur/ntime;

gr = zeros(ntime, 3);
b0 = zeros(npos, 1);
pr = rand(npos, 3);
m0 = [zeros(2,npos); ones(1,npos)];

tic
result = bloch_sim_mex(b1, gr, tp, b0, pr, 10000, 10000, sens, m0);
result';
toc
%%
x= rand(ntime, 8) + 1j*rand(ntime, 8);
y= rand(npos, 8) + 1j*rand(npos, 8);
tic
z = x * y';
toc
% tic
% c = mtimesx(x, y');
% toc
%% compare runtime

e_sens = rand(npos,1) + 1j*rand(npos,1);
e_b1 = rand(ntime,1)+ 1j*rand(ntime,1);
e_gr = rand(ntime, 3);
e_pr = rand(npos, 3);
e_b0 = rand(npos, 1);
e_tp = rand(ntime, 1);
GAMMA_T = 25656565;
e_m0 = rand(npos, 3);

tic
e_b1comb = e_b1 * e_sens';
toc

tic 
    rotz = (e_gr * e_pr' + e_b0') .* e_tp * -1.0 * GAMMA_T;
    rotx = real(e_b1comb) .* e_tp * -1.0 * GAMMA_T;
    roty = imag(e_b1comb) .* e_tp * GAMMA_T; 
    e_m0t = e_m0';
toc

%% with ppl
clc
mex COMPFLAGS='$COMPFLAGS /DUSE_PPL' blochCim.cpp
mex COMPFLAGS='$COMPFLAGS /openmp' blochCim.cpp
mex CXXFLAGS='$CXXFLAGS -std=c++11 -fno-math-errno -ffast-math -march=native' blochCim.cpp

tic
blochCim(e_b1, e_gr, 0.01, 1, 2, e_b0, e_pr, 4, e_sens');
toc

%% GPU
clear
gamma = 267522187.44;
flipangle = 20 * pi/180;
ntime = 389; % number of samples
dur = 2; % second
npos = 28500;

e_sens = rand(npos,1) + 1j*rand(npos,1);
e_b1 = rand(ntime,1)+ 1j*rand(ntime,1);
e_gr = rand(ntime, 3);
e_pr = rand(npos, 3);
e_b0 = rand(npos, 1);
e_tp = rand(ntime, 1);
GAMMA_T = 25656565;
e_m0 = rand(npos, 3);

ge_sens = gpuArray(e_sens);
ge_b1 = gpuArray(e_b1);
ge_gr = gpuArray(e_gr);
ge_pr = gpuArray(e_pr);
ge_b0 = gpuArray(e_b0);
ge_tp = gpuArray(e_tp);
gGAMMA_T = gpuArray(GAMMA_T);
ge_m0 = gpuArray(e_m0);
tic
ge_b1comb = ge_b1 * ge_sens';
% e_b1comb = gather(ge_b1comb);
toc

tic
    rotz = (ge_gr * ge_pr' + ge_b0') .* ge_tp * -1.0 * gGAMMA_T;
    rotx = real(ge_b1comb) .* ge_tp * -1.0 * gGAMMA_T;
    roty = imag(ge_b1comb) .* ge_tp * gGAMMA_T; 
    e_m0t = ge_m0';
    
%    e_m0t = gather(ge_m0t);
%     rotz = gather(grotz);
%     rotx = gather(grotx);
%     roty = gather(groty);
toc

%%
tic
phi = gpuArray(sqrt(rotx.*rotx + roty.*roty + rotz.*rotz));
hp = phi/2;
cp = cos(hp);
sp = sin(hp) ./ phi;
ar = cp;
ai = -rotz .* sp;
br =  roty .* sp;
bi = -rotx .* sp;

arar  = ar.*ar;
aiai  = ai.*ai;
arai2 = 2*ar.*ai;
brbr  = br.*br;
bibi  = bi.*bi;
brbi2 = 2*br.*bi;
arbi2 = 2*ar.*bi;
aibr2 = 2*ai.*br;
arbr2 = 2*ar.*br;
aibi2 = 2*ai.*bi;

toc
%%
for cpos=1:npos
    for ct=1:ntime
        if phi == 0.0
            rotmat = gpuArray(eye(3));
        else 
            
            rotmat(1,1) =  arar  -aiai -brbr +bibi;
            rotmat(1,2) = -arai2 -brbi2;
            rotmat(1,3) = -arbr2 +aibi2;
            rotmat(2,1) =  arai2 -brbi2;
            rotmat(2,2) =  arar  -aiai +brbr -bibi;
            rotmat(2,3) = -aibr2 -arbi2;
            rotmat(3,1) =  arbr2 +aibi2;
            rotmat(3,2) =  arbi2 -aibr2;
            rotmat(3,3) =  arar  +aiai -brbr -bibi;
        end

        m1 = rotmat * e_m0t(:,cpos);
        
        e_m0t(1,cpos) = m1(1);
        e_m0t(2,cpos) = m1(2) ;
        e_m0t(3,cpos) = m1(3) ;
    end
end
toc