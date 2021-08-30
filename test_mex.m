%% compile with Microsoft Visual C++ 2019 
clc
mex bloch_sim_mex.cpp bloch_sim.cpp -R2018a

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
 
%% compile with MinGW64 Compiler (C++) and openmp
% see https://uk.mathworks.com/matlabcentral/answers/279171
clc
mex CXXFLAGS='$CXXFLAGS -std=c++11 -fno-math-errno -ffast-math -march=native -fopenmp' LDFLAGS=-fopenmp bloch_sim_mex.cpp bloch_sim.cpp -R2018a
 %% with GPU - needs to install updated graphic driver
 % download CUDA from https://developer.nvidia.com/cuda-downloads
 % include: C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\include
 % lib: C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64
clc
cd('D:\OneDrive - University College London\Matlab\B1GUI\ipopt_lib\bloch_simulator')
% NVCCFLAGS='$NVCCFLAGS
mexcuda -DUSE_GPU -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart ...
        -DMKL_LP64 -I'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/'...
        -L'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64/' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ...
        -L'C:/Program Files (x86)/Intel/oneAPI/compiler/2021.3.0/windows/compiler/lib/intel64_win' -llibiomp5md...
        bloch_sim_mex.cpp bloch_sim.cpp ./gpu_matrix_mul/gpu_matrix_mul.cu -R2018a
    
 %%
 mexcuda -DUSE_GPU -L'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.4\lib\x64' -lcublas -lcudart ...
         bloch_sim_mex.cpp bloch_sim.cpp ./gpu_matrix_mul/gpu_matrix_mul.cu -R2018a
   
%% test speed
% Constants 
gamma   = 267522187.44;
gamma_hz= gamma/2/pi;
FA      = 45 * pi/180;
ntime   = 400; % number of samples
rf_len  = 4e-3; % second
td      = rf_len/ntime; %  dwell time 
ncoil   = 8;
voxel_sz= [3, 3, 3] * 1e-3;
fov     = [-60 45; -30 70; -90 90] * 1e-3;
rf_tbw  = 12;
is_sinc = true;
is_selective = true;


% positions
[prx, pry, prz] = ndgrid(fov(1,1):voxel_sz(1):fov(1,2), fov(2,1):voxel_sz(2):fov(2,2), fov(3,1):voxel_sz(3):fov(3,2));
sz      = size(prz);
pr      = transpose([prx(:), pry(:), prz(:)]);
npos    = size(pr, 2); % number of spatial positions

% off-resonance
b0 = zeros(npos, 1);

% sinc pulse
if is_sinc
    t = linspace(-rf_len/2, rf_len/2, ntime); % must be column
    BW = rf_tbw/rf_len;
    x = pi*t*BW + eps; % didn't use 2pi -> we are interested in full BW not only the positive
    snc = sin(x) ./ x; 
    hamming_window = 0.53836 + 0.46164*cos(2*pi * linspace(-0.5,0.5,ntime));
    rf = transpose(snc .* hamming_window);
    rf = repmat(rf / sum(rf), [1 ncoil]); % normalize
    b1 = complex(rf * FA/gamma/td/ncoil);
else
    b1 = complex(ones(ntime, ncoil) * FA/gamma/rf_len/ncoil) ;
end
sens = complex(ones(ncoil, npos));

% gradients
gr = zeros(3, ntime);
if is_selective
    BW = rf_tbw/rf_len; % [Hz]
    gr(3, :) = BW / voxel_sz(3) / gamma_hz;
end


m0 = [zeros(2, npos); ones(1, npos)];

% run
b1 = real(b1) + 1i*real(b1);
b1 = b1 / sqrt(2);
result = bloch_sim_mex(b1, gr, td, b0, pr, 10000, 10000, sens, m0);
result_xy = reshape(result(1,:), sz) + 1j*reshape(result(2,:), sz);
result_z = reshape(result(3,:), sz);

%vin(abs(result_xy))
close all
vin(reshape(result(1,:), sz))
vin(reshape(result(2,:), sz))
%%
tic
z = b1 * sens;
toc

tic 
ge_b1 = gpuArray(b1);
ge_gr = gpuArray(sens);
z  = ge_b1 * ge_gr;
z2 = gather(z);
toc

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