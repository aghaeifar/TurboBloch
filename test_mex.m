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

bloch_sim_mex(b1, gr, tp, b0, pr);


%% compare runtime
clear
ntime = 800; % number of samples
npos = 100000; 

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
e_b1comb = e_b1 * e_sens';
toc

tic 
    rotz = (e_gr * e_pr' + e_b0') .* e_tp * -1.0 * GAMMA_T;
    rotx = real(e_b1comb) .* e_tp * -1.0 * GAMMA_T;
    roty = imag(e_b1comb) .* e_tp * GAMMA_T; 
    e_m0t = e_m0';
toc

tic
ge_b1comb = ge_b1 * ge_sens';
toc

tic
    grotz = (ge_gr * ge_pr' + ge_b0') .* ge_tp * -1.0 * gGAMMA_T;
    grotx = real(ge_b1comb) .* ge_tp * -1.0 * gGAMMA_T;
    groty = imag(ge_b1comb) .* ge_tp * gGAMMA_T; 
    ge_m0t = ge_m0';
toc

%%
tic
for cpos=1:npos
    for ct=1:ntime
        phi = sqrt(rotx(ct,cpos)*rotx(ct,cpos) + roty(ct,cpos)*roty(ct,cpos) + rotz(ct,cpos)*rotz(ct,cpos));
        if phi == 0.0
            rotmat = eye(3);
        else
            hp = phi/2;
            cp = cos(hp);
            sp = sin(hp)/phi;
            ar = cp;
            ai = -rotz(ct,cpos)*sp;
            br = roty(ct,cpos)*sp;
            bi = -rotx(ct,cpos)*sp;
            
            arar  = ar*ar;
            aiai  = ai*ai;
            arai2 = 2*ar*ai;
            brbr  = br*br;
            bibi  = bi*bi;
            brbi2 = 2*br*bi;
            arbi2 = 2*ar*bi;
            aibr2 = 2*ai*br;
            arbr2 = 2*ar*br;
            aibi2 = 2*ai*bi;
            
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